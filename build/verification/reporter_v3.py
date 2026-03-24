import argparse
import sqlite3
from pathlib import Path

import polars as pl


SKIP_SUBFIELDS = {
    "ddt_w_adv_pc",
    "ddt_vn_apc_pc",
}


def load_data(db_path: str) -> pl.DataFrame:
    conn = sqlite3.connect(db_path)
    rows = conn.execute(
        "SELECT phys, field, sub_field, ss_a, ss_b, "
        "SNR_db, vSNR_db "
        "FROM comparisons ORDER BY phys"
    ).fetchall()
    conn.close()
    return pl.DataFrame(rows, schema=[
        "phys", "field", "sub_field", "ss_a", "ss_b",
        "SNR_db", "vSNR_db",
    ], orient="row")


def build_summary(df: pl.DataFrame) -> pl.DataFrame:
    enriched = df.with_columns(
        pl.when(pl.col("sub_field") == "-")
        .then(pl.col("field"))
        .otherwise(pl.col("field") + " % " + pl.col("sub_field"))
        .alias("name"),
        (pl.col("ss_a").cast(pl.Utf8) + "v" + pl.col("ss_b").cast(pl.Utf8))
        .alias("pair"),
    )

    worst_snr = (
        enriched.sort(["phys", "SNR_db"])
        .group_by("phys", maintain_order=True).first()
        .select([
            "phys",
            pl.col("SNR_db").alias("worst_SNR"),
            pl.when(pl.col("SNR_db").is_infinite())
            .then(pl.lit("-"))
            .otherwise(pl.col("name") + " @ " + pl.col("pair"))
            .alias("SNR_from"),
        ])
    )

    worst_vsnr = (
        enriched.sort(["phys", "vSNR_db"])
        .group_by("phys", maintain_order=True).first()
        .select([
            "phys",
            pl.col("vSNR_db").alias("worst_vSNR"),
            pl.when(pl.col("vSNR_db").is_infinite())
            .then(pl.lit("-"))
            .otherwise(pl.col("name") + " @ " + pl.col("pair"))
            .alias("vSNR_from"),
        ])
    )

    return worst_snr.join(worst_vsnr, on="phys").sort("phys")


def build_per_field_summary(df: pl.DataFrame) -> pl.DataFrame:
    enriched = df.with_columns(
        pl.when(pl.col("sub_field") == "-")
        .then(pl.col("field"))
        .otherwise(pl.col("field") + " % " + pl.col("sub_field"))
        .alias("name"),
        (pl.col("ss_a").cast(pl.Utf8) + "v" + pl.col("ss_b").cast(pl.Utf8))
        .alias("pair"),
    )

    worst_snr = (
        enriched.sort(["phys", "name", "SNR_db"])
        .group_by(["phys", "name"], maintain_order=True).first()
        .select([
            "phys", "name",
            pl.col("SNR_db").alias("worst_SNR"),
            pl.when(pl.col("SNR_db").is_infinite())
            .then(pl.lit("-"))
            .otherwise(pl.col("pair"))
            .alias("SNR_pair"),
        ])
    )

    worst_vsnr = (
        enriched.sort(["phys", "name", "vSNR_db"])
        .group_by(["phys", "name"], maintain_order=True).first()
        .select([
            "phys", "name",
            pl.col("vSNR_db").alias("worst_vSNR"),
            pl.when(pl.col("vSNR_db").is_infinite())
            .then(pl.lit("-"))
            .otherwise(pl.col("pair"))
            .alias("vSNR_pair"),
        ])
    )

    return worst_snr.join(worst_vsnr, on=["phys", "name"]).sort(["phys", "name"])


def df_to_md(df: pl.DataFrame, title: str) -> str:
    lines = [f"## {title}", ""]
    cols = df.columns
    lines.append("| " + " | ".join(cols) + " |")
    lines.append("| " + " | ".join("---" for _ in cols) + " |")
    for row in df.iter_rows():
        cells = []
        for v in row:
            if isinstance(v, float):
                if v == float("inf"):
                    cells.append("-")
                elif v == float("-inf"):
                    cells.append("-inf")
                else:
                    cells.append(f"{v:.2f}")
            else:
                cells.append(str(v))
        lines.append("| " + " | ".join(cells) + " |")
    lines.append("")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Report convergence from v3 SQLite DB.")
    parser.add_argument("db", help="Path to convergence.db")
    parser.add_argument("-o", "--output", default=None,
                        help="Output markdown report (default: <db>.report.md)")
    parser.add_argument("-p", "--phys", type=int, nargs="*", default=None)
    parser.add_argument("--per-field", action="store_true",
                        help="Show worst SNR per (phys, field) instead of aggregate")
    args = parser.parse_args()

    db_path = args.db
    out_path = args.output or db_path.replace(".db", ".report.md")

    df = load_data(db_path)
    if df.is_empty():
        print("No data in DB.")
        return

    if args.phys:
        df = df.filter(pl.col("phys").is_in(args.phys))

    df = df.filter(~pl.col("sub_field").is_in(SKIP_SUBFIELDS))

    names = df.with_columns(
        pl.when(pl.col("sub_field") == "-")
        .then(pl.col("field"))
        .otherwise(pl.col("field") + " % " + pl.col("sub_field"))
        .alias("name"),
    ).select("name").unique().sort("name")
    print(f"\nComparing {names.height} fields: {', '.join(names['name'].to_list())}")

    if args.per_field:
        table = build_per_field_summary(df)
        title = "Per-Field Worst-Case SNR"
    else:
        table = build_summary(df)
        title = "Worst-Case SNR per Physics Step"

    print()
    with pl.Config(tbl_rows=500, tbl_width_chars=220, tbl_cols=-1, float_precision=2):
        print(table)

    md = f"# Convergence Report\n\nDB: `{db_path}`\n\n" + df_to_md(table, title)
    Path(out_path).write_text(md)
    print(f"\nSaved to: {out_path}")


if __name__ == "__main__":
    main()

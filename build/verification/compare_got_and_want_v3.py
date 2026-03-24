import argparse
import filecmp
import itertools
import os
import re
import sqlite3
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import zstandard as zstd
from tqdm import tqdm

DEFAULT_WORKERS = int(
    os.environ.get(
        "SLURM_CPUS_PER_TASK", os.environ.get("SLURM_CPUS_ON_NODE", os.cpu_count() or 1)
    )
)

# ---------------------------------------------------------------------------
# SQLite schema
# ---------------------------------------------------------------------------

SCHEMA = """
CREATE TABLE IF NOT EXISTS comparisons (
    phys        INTEGER NOT NULL,
    field       TEXT    NOT NULL,
    sub_field   TEXT    NOT NULL,
    ss_a        INTEGER NOT NULL,
    ss_b        INTEGER NOT NULL,
    status      TEXT    NOT NULL,
    max_abs     REAL,
    max_rel     REAL,
    mae         REAL,
    rmse        REAL,
    SNR_db      REAL,
    vSNR_db     REAL,
    mismatches  INTEGER,
    total_elements INTEGER,
    PRIMARY KEY (phys, field, sub_field, ss_a, ss_b)
);
CREATE INDEX IF NOT EXISTS idx_phys_ss ON comparisons (phys, ss_a, ss_b);
"""


def init_db(db_path: Path) -> sqlite3.Connection:
    conn = sqlite3.connect(str(db_path))
    conn.executescript(SCHEMA)
    conn.execute("PRAGMA journal_mode=WAL")
    return conn


def pair_exists(conn: sqlite3.Connection, phys: int, field: str, ss_a: int, ss_b: int) -> bool:
    row = conn.execute(
        "SELECT 1 FROM comparisons WHERE phys=? AND field=? AND ss_a=? AND ss_b=? LIMIT 1",
        (phys, field, ss_a, ss_b),
    ).fetchone()
    return row is not None


def insert_rows(conn: sqlite3.Connection, rows: List[Dict[str, Any]]):
    if not rows:
        return
    conn.executemany(
        """INSERT OR REPLACE INTO comparisons
           (phys, field, sub_field, ss_a, ss_b, status,
            max_abs, max_rel, mae, rmse, SNR_db, vSNR_db,
            mismatches, total_elements)
           VALUES (:phys, :field, :sub_field, :ss_a, :ss_b, :status,
                   :max_abs, :max_rel, :mae, :rmse, :SNR_db, :vSNR_db,
                   :mismatches, :total_elements)""",
        rows,
    )
    conn.commit()


# ---------------------------------------------------------------------------
# File I/O
# ---------------------------------------------------------------------------

METADATA_TAGS = {"assoc", "alloc", "rank", "size", "lbound", "entries", "missing"}
_TAG_RE = re.compile(r"^#\s+(\S+)", re.MULTILINE)


def read_text(path: Path) -> str:
    if path.suffix == ".zst":
        with open(path, "rb") as f:
            return zstd.ZstdDecompressor().decompress(f.read()).decode("utf-8", errors="replace")
    return path.read_text(encoding="utf-8", errors="replace")


def parse_serde_file(path: Path) -> Dict[str, np.ndarray]:
    text = read_text(path)
    results: Dict[str, np.ndarray] = {}

    tags = []
    for m in _TAG_RE.finditer(text):
        eol = text.find("\n", m.end())
        data_start = eol + 1 if eol != -1 else len(text)
        tags.append((m.start(), m.group(1), data_start))

    if not tags:
        stem = path.stem.split(".")[0]
        try:
            arr = np.fromstring(text, sep="\n", dtype=np.float64)
            if arr.size > 0:
                results[stem] = arr
        except ValueError:
            pass
        return results

    current_field: Optional[str] = None
    for i, (_, tag_name, data_start) in enumerate(tags):
        data_end = tags[i + 1][0] if i + 1 < len(tags) else len(text)
        if tag_name == "entries":
            if current_field is not None:
                block = text[data_start:data_end]
                arr = np.fromstring(block, sep="\n", dtype=np.float64)
                if arr.size > 0:
                    if current_field in results:
                        results[current_field] = np.concatenate([results[current_field], arr])
                    else:
                        results[current_field] = arr
        elif tag_name not in METADATA_TAGS:
            current_field = tag_name

    return results


def count_entries(path: Path) -> int:
    text = read_text(path)
    count = 0
    in_entries = False
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("# "):
            tag = stripped[2:].strip().split()[0]
            in_entries = tag == "entries"
        elif in_entries and stripped:
            count += 1
    return count


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------


def compare_data(
    got: np.ndarray,
    want: np.ndarray,
    abs_tol: float = 1e-12,
    rel_tol: float = 1e-12,
) -> Dict[str, Any]:
    if got.size != want.size:
        return {
            "ok": False, "max_abs": -1.0, "max_rel": -1.0, "mae": -1.0,
            "rmse": -1.0, "SNR_db": -1.0, "vSNR_db": -1.0,
            "mismatches": abs(got.size - want.size),
            "count": max(got.size, want.size),
        }

    g = got.astype(np.float64)
    w = want.astype(np.float64)

    abs_diff = np.abs(g - w)
    scale = np.fmax(np.fmax(np.abs(g), np.abs(w)), abs_tol)
    with np.errstate(divide="ignore", invalid="ignore"):
        rel_diff = abs_diff / scale

    is_close = abs_diff <= np.maximum(
        rel_tol * np.maximum(np.abs(g), np.abs(w)), abs_tol
    )
    mismatches = np.count_nonzero(~is_close)

    norm_diff = np.linalg.norm(abs_diff)
    norm_want = np.linalg.norm(w)

    if norm_diff == 0:
        snr_db = float("inf")
    elif norm_want == 0:
        snr_db = -float("inf")
    else:
        snr_db = 20 * np.log10(norm_want / norm_diff)

    var_want = np.var(w)
    var_err = np.var(g - w)
    if var_err == 0:
        vsnr_db = float("inf")
    elif var_want == 0:
        vsnr_db = -float("inf")
    else:
        vsnr_db = 10 * np.log10(var_want / var_err)

    return {
        "ok": mismatches == 0,
        "max_abs": float(np.max(abs_diff)),
        "max_rel": float(np.max(rel_diff)),
        "mae": float(np.mean(abs_diff)),
        "rmse": float(norm_diff / np.sqrt(g.size)) if g.size > 0 else 0.0,
        "SNR_db": float(snr_db),
        "vSNR_db": float(vsnr_db),
        "mismatches": int(mismatches),
        "count": int(g.size),
    }


def compare_one_pair(
    field: str,
    path_a: Path,
    path_b: Path,
    atol: float,
    rtol: float,
) -> List[Dict[str, Any]]:
    # Fast path: byte-identical files
    if filecmp.cmp(str(path_a), str(path_b), shallow=False):
        total = count_entries(path_a)
        return [{
            "sub_field": "-",
            "status": "OK",
            "max_abs": 0.0, "max_rel": 0.0, "mae": 0.0, "rmse": 0.0,
            "SNR_db": float("inf"), "vSNR_db": float("inf"),
            "mismatches": 0, "total_elements": total,
        }]

    a_fields = parse_serde_file(path_a)
    b_fields = parse_serde_file(path_b)

    rows = []
    all_vars = sorted(set(a_fields) | set(b_fields))
    has_subfields = len(all_vars) > 1

    for var in all_vars:
        sf = var if has_subfields else "-"
        if var not in a_fields or var not in b_fields:
            rows.append({
                "sub_field": sf, "status": "MISSING",
                "max_abs": -1.0, "max_rel": -1.0, "mae": -1.0, "rmse": -1.0,
                "SNR_db": -1.0, "vSNR_db": -1.0,
                "mismatches": -1, "total_elements": -1,
            })
            continue

        stats = compare_data(a_fields[var], b_fields[var], abs_tol=atol, rel_tol=rtol)
        rows.append({
            "sub_field": sf,
            "status": "OK" if stats["ok"] else "DIFF",
            "max_abs": stats["max_abs"], "max_rel": stats["max_rel"],
            "mae": stats["mae"], "rmse": stats["rmse"],
            "SNR_db": stats["SNR_db"], "vSNR_db": stats["vSNR_db"],
            "mismatches": stats["mismatches"], "total_elements": stats["count"],
        })

    return rows


# ---------------------------------------------------------------------------
# File discovery
# ---------------------------------------------------------------------------

SKIP_FIELDS = {"dtime", "istep", "ntnd", "lvn_only", "ldeepatmo", "dt_linintp_ubc",
               "global_data", "p_patch", "p_int",
               "z_kin_hor_e", "z_vt_ie", "z_w_concorr_me"}

FILE_PAT = re.compile(
    r"^(?P<field>.+?)\.t0\.p(?P<phys>\d+)\.d(?P<d>\d+)\.vt(?P<vt>\d+)\.ss(?P<ss>\d+)\.data(?:\.zst)?$"
)


def discover_files(d: Path) -> List[dict]:
    records = []
    for f in d.iterdir():
        m = FILE_PAT.match(f.name)
        if m:
            field = m.group("field")
            if field.split(".")[0] in SKIP_FIELDS:
                continue
            records.append({
                "field": field,
                "phys": int(m.group("phys")),
                "ss": int(m.group("ss")),
                "path": f,
            })
    return records


# ---------------------------------------------------------------------------
# Worker: compare all pairs for one (phys, field) group
#   Reads each ss file once, then computes all pairwise comparisons in-memory.
# ---------------------------------------------------------------------------


def group_worker(
    task: Tuple[int, str, Dict[int, Path], List[Tuple[int, int]]],
    atol: float,
    rtol: float,
) -> List[Dict[str, Any]]:
    phys, field, ss_paths, pairs_needed = task

    # Read and parse each ss file once
    parsed: Dict[int, Dict[str, np.ndarray]] = {}
    raw_bytes: Dict[int, bytes] = {}
    for ss, path in ss_paths.items():
        with open(path, "rb") as f:
            raw_bytes[ss] = f.read()
        if path.suffix == ".zst":
            text = zstd.ZstdDecompressor().decompress(raw_bytes[ss]).decode("utf-8", errors="replace")
        else:
            text = raw_bytes[ss].decode("utf-8", errors="replace")
        parsed[ss] = _parse_text(text)

    results = []
    for ss_a, ss_b in pairs_needed:
        # Fast path: byte-identical
        if raw_bytes[ss_a] == raw_bytes[ss_b]:
            total = sum(a.size for a in parsed[ss_a].values())
            results.append({
                "phys": phys, "field": field, "sub_field": "-",
                "ss_a": ss_a, "ss_b": ss_b, "status": "OK",
                "max_abs": 0.0, "max_rel": 0.0, "mae": 0.0, "rmse": 0.0,
                "SNR_db": float("inf"), "vSNR_db": float("inf"),
                "mismatches": 0, "total_elements": total,
            })
            continue

        a_fields = parsed[ss_a]
        b_fields = parsed[ss_b]
        all_vars = sorted(set(a_fields) | set(b_fields))
        has_subfields = len(all_vars) > 1

        for var in all_vars:
            sf = var if has_subfields else "-"
            if var not in a_fields or var not in b_fields:
                results.append({
                    "phys": phys, "field": field, "sub_field": sf,
                    "ss_a": ss_a, "ss_b": ss_b, "status": "MISSING",
                    "max_abs": -1.0, "max_rel": -1.0, "mae": -1.0, "rmse": -1.0,
                    "SNR_db": -1.0, "vSNR_db": -1.0,
                    "mismatches": -1, "total_elements": -1,
                })
                continue
            stats = compare_data(a_fields[var], b_fields[var], abs_tol=atol, rel_tol=rtol)
            results.append({
                "phys": phys, "field": field, "sub_field": sf,
                "ss_a": ss_a, "ss_b": ss_b,
                "status": "OK" if stats["ok"] else "DIFF",
                "max_abs": stats["max_abs"], "max_rel": stats["max_rel"],
                "mae": stats["mae"], "rmse": stats["rmse"],
                "SNR_db": stats["SNR_db"], "vSNR_db": stats["vSNR_db"],
                "mismatches": stats["mismatches"], "total_elements": stats["count"],
            })

    return results


def _parse_text(text: str) -> Dict[str, np.ndarray]:
    """Parse serde text (already decompressed) into {field: ndarray}."""
    results: Dict[str, np.ndarray] = {}
    tags = []
    for m in _TAG_RE.finditer(text):
        eol = text.find("\n", m.end())
        data_start = eol + 1 if eol != -1 else len(text)
        tags.append((m.start(), m.group(1), data_start))

    if not tags:
        try:
            arr = np.fromstring(text, sep="\n", dtype=np.float64)
            if arr.size > 0:
                results["root"] = arr
        except ValueError:
            pass
        return results

    current_field: Optional[str] = None
    for i, (_, tag_name, data_start) in enumerate(tags):
        data_end = tags[i + 1][0] if i + 1 < len(tags) else len(text)
        if tag_name == "entries":
            if current_field is not None:
                block = text[data_start:data_end]
                arr = np.fromstring(block, sep="\n", dtype=np.float64)
                if arr.size > 0:
                    if current_field in results:
                        results[current_field] = np.concatenate([results[current_field], arr])
                    else:
                        results[current_field] = arr
        elif tag_name not in METADATA_TAGS:
            current_field = tag_name
    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Pairwise comparison of ICON velocity fields across substep counts. "
                    "Results stored in SQLite for checkpointing."
    )
    parser.add_argument("root", nargs="?", default="experiments/exclaim_ape_R02B04_dt8",
                        help="Directory containing .data[.zst] files")
    parser.add_argument("-o", "--output", default=None,
                        help="SQLite output path (default: <root>/convergence.db)")
    parser.add_argument("-p", "--phys", type=int, nargs="*", default=None,
                        help="Physics steps to compare (default: all)")
    parser.add_argument("--atol", type=float, default=1e-12)
    parser.add_argument("--rtol", type=float, default=1e-12)
    parser.add_argument("-j", "--workers", type=int, default=DEFAULT_WORKERS)
    args = parser.parse_args()

    root = Path(args.root)
    db_path = Path(args.output) if args.output else root / "convergence.db"

    records = discover_files(root)
    if not records:
        print(f"No matching .data files in {root}")
        return

    all_ss = sorted(set(r["ss"] for r in records))
    all_phys = sorted(set(r["phys"] for r in records))
    phys_steps = args.phys if args.phys else all_phys

    print(f"Root: {root}")
    print(f"DB:   {db_path}")
    print(f"Available ss: {all_ss}")
    print(f"Available phys: {all_phys}")
    print(f"Comparing phys: {phys_steps}")
    ss_pairs = list(itertools.combinations(all_ss, 2))
    print(f"SS pairs: {len(ss_pairs)} unique per (phys, field)")

    # Build lookup: (phys, field, ss) -> path
    lookup: Dict[Tuple[int, str, int], Path] = {}
    for r in records:
        lookup[(r["phys"], r["field"], r["ss"])] = r["path"]

    conn = init_db(db_path)

    # Build group tasks: one per (phys, field), containing all needed ss pairs
    group_tasks = []
    skipped = 0
    for phys in phys_steps:
        fields_at_phys = sorted(set(
            r["field"] for r in records if r["phys"] == phys
        ))
        for field in fields_at_phys:
            # Collect available ss paths for this (phys, field)
            ss_paths: Dict[int, Path] = {}
            for ss in all_ss:
                key = (phys, field, ss)
                if key in lookup:
                    ss_paths[ss] = lookup[key]

            # Filter to pairs that still need computing
            pairs_needed = []
            for ss_a, ss_b in ss_pairs:
                if ss_a not in ss_paths or ss_b not in ss_paths:
                    continue
                if pair_exists(conn, phys, field, ss_a, ss_b):
                    skipped += 1
                    continue
                pairs_needed.append((ss_a, ss_b))

            if pairs_needed:
                group_tasks.append((phys, field, ss_paths, pairs_needed))

    total_pairs = sum(len(t[3]) for t in group_tasks)
    if not group_tasks:
        print(f"All pairs already computed ({skipped} skipped). Nothing to do.")
    else:
        print(f"Groups: {len(group_tasks)} (phys, field) groups, {total_pairs} pairs to compute")
        if skipped:
            print(f"Skipped: {skipped} already in DB")
        print(f"Workers: {args.workers}\n")

        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futures = {
                ex.submit(partial(group_worker, atol=args.atol, rtol=args.rtol), t): t
                for t in group_tasks
            }
            pbar = tqdm(total=len(futures), desc="Comparing", unit="group")
            batch = []
            for fut in futures:
                _, field, _, _ = futures[fut]
                pbar.set_postfix_str(field, refresh=True)
                try:
                    rows = fut.result()
                    batch.extend(rows)
                except Exception as e:
                    tqdm.write(f"Error: {field}: {e}")
                pbar.update(1)
                # Flush periodically for checkpointing
                if len(batch) >= 200:
                    insert_rows(conn, batch)
                    batch = []
            if batch:
                insert_rows(conn, batch)
            pbar.close()

    # Print summary
    print_summary(conn, phys_steps)
    conn.close()
    print(f"\nResults in: {db_path}")


def print_summary(conn: sqlite3.Connection, phys_steps: List[int]):
    import polars as pl

    rows = conn.execute(
        "SELECT phys, field, sub_field, ss_a, ss_b, status, "
        "max_abs, max_rel, rmse, SNR_db, vSNR_db, mismatches, total_elements "
        "FROM comparisons ORDER BY phys, field, sub_field, ss_a, ss_b"
    ).fetchall()

    if not rows:
        print("No results yet.")
        return

    df = pl.DataFrame(rows, schema=[
        "phys", "field", "sub_field", "ss_a", "ss_b", "status",
        "max_abs", "max_rel", "rmse", "SNR_db", "vSNR_db",
        "mismatches", "total_elements",
    ], orient="row")

    def fmt_count(n: int) -> str:
        if n < 0:
            return "?"
        if n >= 1_000_000:
            return f"{n / 1_000_000:.2f}M"
        if n >= 1_000:
            return f"{n / 1_000:.1f}K"
        return str(n)

    display = (
        df.with_columns(
            pl.when(pl.col("sub_field") == "-")
            .then(pl.col("field"))
            .otherwise(pl.col("field") + " % " + pl.col("sub_field"))
            .alias("name"),
            (pl.lit("ss") + pl.col("ss_a").cast(pl.Utf8) + pl.lit(" v ") +
             pl.lit("ss") + pl.col("ss_b").cast(pl.Utf8)).alias("pair"),
            pl.col("total_elements").map_elements(fmt_count, return_dtype=pl.Utf8).alias("N"),
            pl.col("mismatches").map_elements(fmt_count, return_dtype=pl.Utf8).alias("mis-\nmatch"),
        )
        .select([
            "phys", "name", "pair", "status",
            "max_abs", "max_rel", "rmse",
            "SNR_db", "vSNR_db",
            "mis-\nmatch", "N",
        ])
        .rename({
            "max_abs": "max\nabs",
            "max_rel": "max\nrel",
            "SNR_db": "SNR\n(dB)",
            "vSNR_db": "vSNR\n(dB)",
        })
        .sort(["phys", "name", "pair"])
    )

    print("\n" + "=" * 100)
    print("PAIRWISE CONVERGENCE SUMMARY")
    print("=" * 100)
    with pl.Config(tbl_rows=500, tbl_width_chars=220, tbl_cols=-1, float_precision=2):
        if display["phys"].n_unique() == 1:
            print(display.drop("phys"))
        else:
            print(display)
    print("=" * 100)


if __name__ == "__main__":
    main()

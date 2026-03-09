import argparse
import os
import re
import sys
from typing import Dict, Iterable, List, Optional, Tuple, Any
import numpy as np
from pathlib import Path
import polars as pl
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor
from functools import partial

sys.stdout.reconfigure(line_buffering=True)

DEFAULT_WORKERS = int(
    os.environ.get(
        "SLURM_CPUS_PER_TASK", os.environ.get("SLURM_CPUS_ON_NODE", os.cpu_count() or 1)
    )
)

KNOWN_METADATA = {
    "assoc", "rank", "size", "lbound", "entries", "alloc", "missing"
}

def _stream_lines(path: Path) -> Iterable[str]:
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            ls = line.rstrip("\n\r")
            if not ls.strip(): continue
            yield ls

def extract_variable(path: Path, target_var: str) -> Tuple[np.ndarray, Tuple[int, ...]]:
    found_target = False
    rank, shape, data = 0, [], []
    current_tag = None
    for line in _stream_lines(path):
        if line.startswith("# "):
            tag = line.lstrip("# ").strip().split()[0]
            if tag and tag not in KNOWN_METADATA:
                if found_target: break
                if tag == target_var: found_target = True
            if found_target: current_tag = tag
            continue
        if found_target:
            try:
                val = float(line)
                if current_tag == "rank":
                    rank = int(val)
                elif current_tag == "size":
                    shape.append(int(val))
                elif current_tag in {"entries", target_var}:
                    data.append(val)
            except ValueError:
                continue
    if not data:
        raise ValueError(f"Variable '{target_var}' not found in {path}")
    return np.array(data), tuple(shape) if shape else (len(data),)

def compute_base_stats(g_data: Tuple[np.ndarray, Any], w_data: Tuple[np.ndarray, Any], atol: float, rtol: float) -> Tuple[Dict[str, Any], np.ndarray]:
    g, g_shape = g_data
    w, w_shape = w_data
    diff_flat = g - w
    abs_diff = np.abs(diff_flat)
    max_err = np.max(abs_diff) if abs_diff.size > 0 else 0.0
    mean_abs_err = np.mean(abs_diff) if abs_diff.size > 0 else 0.0
    mean_signed_err = np.mean(diff_flat) if diff_flat.size > 0 else 0.0
    err_std = np.std(diff_flat, ddof=1) if diff_flat.size > 1 else 0.0

    # Spatial sums for rigorous correlation
    spat_sum_sq, spat_sum_prod, spat_count = 0.0, 0.0, 0
    spatial_r_values = []
    if diff_flat.size > 1 and np.std(diff_flat) > 0:
        try:
            error_field = diff_flat.reshape(w_shape)
            for axis in range(error_field.ndim):
                if error_field.shape[axis] > 1:
                    s1, s2 = [slice(None)] * error_field.ndim, [slice(None)] * error_field.ndim
                    s1[axis], s2[axis] = slice(1, None), slice(0, -1)
                    v1, v2 = error_field[tuple(s1)].flatten(), error_field[tuple(s2)].flatten()
                    spat_sum_prod += np.sum(v1 * v2)
                    spat_sum_sq += np.sum(v1**2 + v2**2)
                    spat_count += v1.size * 2
                    if np.std(v1) > 0 and np.std(v2) > 0:
                        spatial_r_values.append(np.corrcoef(v1, v2)[0, 1])
        except:
            pass

    # Effective Sample Size (N_eff) correction for spatial correlation
    rho_spatial = np.mean(spatial_r_values) if spatial_r_values else 0.0
    rho_spatial = np.clip(rho_spatial, -0.99, 0.99)
    n_eff_spatial = diff_flat.size * (1 - rho_spatial) / (1 + rho_spatial)

    t_stat = (mean_signed_err / (err_std / np.sqrt(max(1, n_eff_spatial)))) if (err_std > 0 and n_eff_spatial > 1) else 0.0
    
    scale = np.fmax(np.fmax(np.abs(g), np.abs(w)), atol)
    is_close = abs_diff <= np.fmax(rtol * np.fmax(np.abs(g), np.abs(w)), atol)
    mismatches = np.count_nonzero(~is_close)
    norm_diff = np.linalg.norm(abs_diff)
    norm_w = np.linalg.norm(w)
    
    max_abs_w = np.max(np.abs(w)) if w.size > 0 else 0.0
    w_range = (np.max(w) - np.min(w)) if w.size > 0 else 0.0
    snr_db = 20 * np.log10(norm_w / norm_diff) if norm_diff > 0 and norm_w > 0 else (float("inf") if norm_diff == 0 else -float("inf"))

    p_corr = np.corrcoef(g, w)[0, 1] if np.std(g) > 0 and np.std(w) > 0 else (1.0 if np.allclose(g, w, atol=atol) else 0.0)
        
    stats = {
        "ok": int(mismatches) == 0,
        "max_abs": max_err,
        "mean_abs": mean_abs_err,
        "rmse": norm_diff / np.sqrt(g.size) if g.size > 0 else 0.0,
        "snr_db": snr_db,
        "nme_max": np.max(abs_diff) / max_abs_w if max_abs_w > 0 else 0.0,
        "nme_range": np.max(abs_diff) / w_range if w_range > 0 else 0.0,
        "bias": mean_signed_err,
        "bias_ratio": mean_signed_err / mean_abs_err if mean_abs_err > 0 else 0.0,
        "error_std": err_std,
        "bias_stddev_ratio": mean_signed_err / err_std if err_std > 0 else 0.0,
        "t_stat": t_stat,
        "spat_r": np.mean(spatial_r_values) if spatial_r_values else 0.0,
        "agg_sums": (spat_sum_prod, spat_sum_sq, spat_count),
        "pearson_corr": p_corr,
        "total": g.size,
        "has_nan": not (np.all(np.isfinite(g)) and np.all(np.isfinite(w)))
    }
    return stats, diff_flat

def worker(task: Tuple[str, int, int, Path, Path], var: str, atol: float, rtol: float) -> Tuple[str, int, int, Dict[str, Any], np.ndarray]:
    c_type, time_t, tstep, curr_path, ref_path = task
    try:
        d_curr, d_ref = extract_variable(curr_path, var), extract_variable(ref_path, var)
        if d_curr[0].size != d_ref[0].size:
            m = min(d_curr[0].size, d_ref[0].size)
            d_curr = (d_curr[0][:m], d_curr[1])
            d_ref  = (d_ref[0][:m], d_ref[1])
        stats, err_field = compute_base_stats(d_curr, d_ref, atol, rtol)
        return c_type, time_t, tstep, stats, err_field
    except Exception as e:
        print(f"Worker Error: {e}")
        return c_type, time_t, tstep, {}, np.array([])

METRIC_DESCRIPTIONS = {
    "snr_db": "Signal-to-Noise Ratio in decibels (20 * log10(||ref|| / ||err||)). Higher is better.",
    "rmse": "Root Mean Square Error. Lower is better.",
    "max_abs": "Maximum Absolute Error (max|g - w|). Lower is better.",
    "mean_abs": "Mean Absolute Error (mean|g - w|). Lower is better.",
    "nme_max": "Normalized Maximum Error w.r.t. Max Absolute Reference (max|g - w| / max|w|).",
    "nme_range": "Normalized Maximum Error w.r.t. Reference Range (max|g - w| / (max(w) - min(w))).",
    "bias": "Mean Signed Error (mean(g - w)). Positive means overestimation, negative means underestimation.",
    "bias_ratio": "Bias Ratio (mean(g - w) / mean|g - w|). Range [-1, 1]. 1 means all error is positive bias, 0 means centered around zero.",
    "error_std": "Standard Deviation of Error (std(g - w)). Measures the 'jitter' or random noise around the mean bias.",
    "t_stat": "One-sample T-statistic using Effective Sample Size (N_eff) to account for correlation. Significant if |t| > 2.",
    "bias_stddev_ratio": "Ratio of systematic bias to random noise (mean_bias / std_error).",
    "error_autocorr": "Autocorrelation of Error. In Time-Series: Mean spatial lag-1 correlation. In --aggregate: Rigorous 4D (Space+Time) correlation.",
    "pearson_corr": "Pearson Correlation Coefficient between signals. High values (~1.0) indicate identical patterns."
}

def show_table(data, title, tsteps, metric, row_label="T"):
    rows = []
    for lbl in sorted(data.keys()):
        row = {row_label: str(lbl)}
        row.update({str(ts): data[lbl][ts] for ts in tsteps})
        rows.append(row)
    df = pl.DataFrame(rows)
    df = df.select([c for c in df.columns if c == row_label or not (df[c] == "-").all()])
    if len(df.columns) > 1:
        comp_cols = [c for c in df.columns if c != row_label]
        mask = pl.any_horizontal([pl.col(c).is_in(["-", "ref"]).not_() for c in comp_cols])
        df = df.filter(mask | (pl.col(row_label).is_in(["Aggregate", "Expectation"])))
    print(f"\n{'='*80}\n{title} (Metric: {metric})\n{'='*80}")
    with pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_width_chars=1000):
        print(df)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root")
    parser.add_argument("-v", "--var", default="vn")
    parser.add_argument("-m", "--metric", default="snr_db")
    parser.add_argument("--aggregate", action="store_true")
    parser.add_argument("--atol", type=float, default=1e-12)
    parser.add_argument("--rtol", type=float, default=1e-12)
    parser.add_argument("-j", "--workers", type=int, default=DEFAULT_WORKERS)
    args = parser.parse_args()
    
    root = Path(args.root)
    if not root.is_dir() or args.metric not in METRIC_DESCRIPTIONS:
        print(f"Error: Invalid root or metric. Available: {', '.join(METRIC_DESCRIPTIONS.keys())}")
        sys.exit(1)
    
    print(f"\n>>> Metric: {args.metric}")
    print(f">>> Description: {METRIC_DESCRIPTIONS[args.metric]}\n")

    files = []
    t_pat = re.compile(r"tstep\.(\d+)")
    f_pat = re.compile(r"p_prog\.t0\.(\d+)\.data")
    for d in root.iterdir():
        if t_pat.match(d.name) and d.is_dir():
            ts = int(t_pat.match(d.name).group(1))
            for f in d.glob("p_prog.t0.*.data"):
                fm = f_pat.match(f.name)
                if fm:
                    files.append({
                        "tstep": ts,
                        "T": ((int(fm.group(1)) - 1) // 5) * ts,
                        "path": f
                    })

    if not files:
        print("No matching files found.")
        return
        
    df_f = pl.DataFrame(files)
    tsteps = sorted(df_f["tstep"].unique().to_list())
    times = sorted(df_f["T"].unique().to_list())

    tasks = []
    for t in times:
        row = df_f.filter(pl.col("T") == t).sort("tstep")
        if not row.is_empty():
            ref = row[0]
            for i, ts in enumerate(tsteps):
                curr = row.filter(pl.col("tstep") == ts)
                if not curr.is_empty():
                    tasks.append(("baseline", t, ts, curr["path"][0], ref["path"][0]))
                    if i > 0:
                        prev = row.filter(pl.col("tstep") == tsteps[i-1])
                        if not prev.is_empty():
                            tasks.append(("consecutive", t, ts, curr["path"][0], prev["path"][0]))

    res_b = {t: {ts: "-" for ts in tsteps} for t in times}
    res_c = {t: {ts: "-" for ts in tsteps} for t in times}
    agg_stats = {(ctype, ts): [0.0, 0.0, 0] for ctype in ["baseline", "consecutive"] for ts in tsteps}
    # (ctype, ts) -> [err_sum, err_sq_sum, total_N, abs_err_sum]
    bias_stats = {(ctype, ts): [0.0, 0.0, 0, 0.0] for ctype in ["baseline", "consecutive"] for ts in tsteps}
    last_err = {}

    tasks.sort(key=lambda x: (x[0], x[2], x[1])) # Type, Tstep, T
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for ctype, t, ts, stats, errf in tqdm(ex.map(partial(worker, var=args.var, atol=args.atol, rtol=args.rtol), tasks), total=len(tasks)):
            if not stats:
                continue
            key = (ctype, ts)
            sp_prod, sp_sq, sp_cnt = stats["agg_sums"]
            agg_stats[key][0] += sp_prod
            agg_stats[key][1] += sp_sq
            agg_stats[key][2] += sp_cnt
            
            if args.metric in {"error_autocorr", "t_stat"}:
                if key in last_err and last_err[key][0] < t:
                    if np.std(errf) > 0 and np.std(last_err[key][1]) > 0:
                        agg_stats[key][0] += np.sum(errf * last_err[key][1])
                        agg_stats[key][1] += np.sum(errf**2 + last_err[key][1]**2)
                        agg_stats[key][2] += errf.size * 2
                last_err[key] = (t, errf)
            
            n = stats["total"]
            bias_stats[key][0] += stats["bias"] * n
            # Sum of Squares = (std^2 * (n-1)) + n*mean^2 
            # This allows us to recompute the global variance later
            ss = (stats["error_std"]**2 * (n - 1)) + (n * stats["bias"]**2) if n > 1 else (n * stats["bias"]**2)
            bias_stats[key][1] += ss
            bias_stats[key][2] += n
            bias_stats[key][3] += stats["mean_abs"] * n
            
            val = stats["spat_r"] if "autocorr" in args.metric else stats.get(args.metric, 0)
            if ctype == "baseline":
                res_b[t][ts] = f"{val:.4f}"
            else:
                res_c[t][ts] = f"{val:.4f}"

    if args.aggregate:
        def get_val(ctype, ts):
            if "autocorr" in args.metric:
                prod, sq, cnt = agg_stats[(ctype, ts)]
                return f"{(2.0 * prod / sq):.4f}" if sq > 0 else "0.0000"
            e_sum, sq_sum, n, a_sum = bias_stats[(ctype, ts)]
            if n == 0:
                return "-"
            m_bias = e_sum / n
            if args.metric == "bias":
                return f"{m_bias:.4f}"
            if args.metric == "bias_ratio":
                m_abs = a_sum / n
                return f"{(m_bias / m_abs):.4f}" if m_abs > 0 else "0.0000"
            if args.metric == "bias_stddev_ratio":
                global_var = (sq_sum / n) - (m_bias**2)
                g_std = np.sqrt(max(0, global_var))
                return f"{(m_bias / g_std):.4f}" if g_std > 0 else "0.0000"
            if args.metric == "error_std":
                global_var = (sq_sum / n) - (m_bias**2)
                return f"{np.sqrt(max(0, global_var)):.4f}"
            if args.metric == "t_stat":
                global_var = (sq_sum / n) - (m_bias**2)
                # Rigorous N_eff calculation for the 4D volume using the 4D rho
                prod, sq, _ = agg_stats[(ctype, ts)]
                global_rho = (2.0 * prod / sq) if sq > 0 else 0.0
                global_rho = np.clip(global_rho, -0.99, 0.99)
                global_n_eff = n * (1 - global_rho) / (1 + global_rho)
                
                stderr = np.sqrt(max(0, global_var)) / np.sqrt(max(1, global_n_eff))
                return f"{(m_bias / stderr):.4f}" if stderr > 0 else "0.0000"
            return f"{m_bias:.4f}"
            
        lbl = "Aggregate" if "autocorr" in args.metric else "Expectation"
        show_table({lbl: {ts: get_val("baseline", ts) for ts in tsteps}}, f"{args.metric} w.r.t. most refined", tsteps, args.metric, "Type")
        show_table({lbl: {ts: get_val("consecutive", ts) for ts in tsteps}}, f"{args.metric} w.r.t. next finer level", tsteps, args.metric, "Type")
    else:
        for t in times:
            r = df_f.filter(pl.col("T") == t).sort("tstep")
            if not r.is_empty():
                res_b[t][r["tstep"][0]] = "ref"
        show_table(res_b, f"{args.metric} w.r.t. most refined", tsteps, args.metric)
        show_table(res_c, f"{args.metric} w.r.t. next finer level", tsteps, args.metric)

if __name__ == "__main__":
    main()

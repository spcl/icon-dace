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

def extract_variable(path: Path, target_var: str) -> np.ndarray:
    data = []
    current_var = re.sub(r"_\d+$", "", path.stem)
    found_target = False
    
    line_stream = _stream_lines(path)
    for line in line_stream:
        if line.startswith("# "):
            tag = line.lstrip("# ").strip().split()[0]
            if tag and tag not in KNOWN_METADATA:
                if found_target: break
                current_var = tag
                if current_var == target_var: found_target = True
            continue
            
        if found_target:
            try:
                data.append(float(line))
            except ValueError:
                continue
                
    if not data:
        raise ValueError(f"Variable '{target_var}' not found in {path}")
    return np.array(data)

METRIC_DESCRIPTIONS = {
    "snr_db": "Signal-to-Noise Ratio in decibels (20 * log10(||ref|| / ||err||)). Higher is better.",
    "rmse": "Root Mean Square Error. Lower is better.",
    "max_abs": "Maximum Absolute Error (max|g - w|). Lower is better.",
    "mean_abs": "Mean Absolute Error (mean|g - w|). Lower is better.",
    "max_rel": "Maximum Pointwise Relative Error (max(|g - w| / max(|g|, |w|, eps))). Can be unstable for small values.",
    "nme_max": "Normalized Maximum Error w.r.t. Max Absolute Reference (max|g - w| / max|w|).",
    "nme_range": "Normalized Maximum Error w.r.t. Reference Range (max|g - w| / (max(w) - min(w))).",
    "bias": "Mean Signed Error (mean(g - w)). Positive means overestimation, negative means underestimation.",
    "bias_ratio": "Bias Ratio (mean(g - w) / mean|g - w|). Range [-1, 1]. 1 means all error is positive bias, 0 means centered around zero.",
    "error_autocorr": "Lag-1 Autocorrelation of Error. High values (~1.0) indicate spatially clustered/smooth errors.",
    "pearson_corr": "Pearson Correlation Coefficient between signals. High values (~1.0) indicate identical patterns."
}

def compute_stats(g: np.ndarray, w: np.ndarray, atol: float, rtol: float) -> Dict[str, Any]:
    diff = g - w
    abs_diff = np.abs(diff)
    max_err = np.max(abs_diff) if abs_diff.size > 0 else 0.0
    mean_abs_err = np.mean(abs_diff) if abs_diff.size > 0 else 0.0
    mean_signed_err = np.mean(diff) if diff.size > 0 else 0.0
    
    scale = np.fmax(np.fmax(np.abs(g), np.abs(w)), atol)
    rel_diff = abs_diff / scale
    is_close = abs_diff <= np.fmax(rtol * np.fmax(np.abs(g), np.abs(w)), atol)
    mismatches = np.count_nonzero(~is_close)
    norm_diff = np.linalg.norm(abs_diff)
    norm_w = np.linalg.norm(w)
    
    max_abs_w = np.max(np.abs(w)) if w.size > 0 else 0.0
    w_range = (np.max(w) - np.min(w)) if w.size > 0 else 0.0
    
    if norm_diff == 0: 
        snr_db = float("inf")
        e_autocorr = 1.0
        p_corr = 1.0
    elif norm_w == 0: 
        snr_db = -float("inf")
        e_autocorr = 0.0
        p_corr = 0.0
    else: 
        snr_db = 20 * np.log10(norm_w / norm_diff)
        # Compute correlations only if variance is non-zero
        if np.std(diff) > 0 and diff.size > 1:
            e_autocorr = np.corrcoef(diff[:-1], diff[1:])[0, 1]
        else:
            e_autocorr = 0.0
            
        if np.std(g) > 0 and np.std(w) > 0:
            p_corr = np.corrcoef(g, w)[0, 1]
        else:
            p_corr = 1.0 if np.allclose(g, w, atol=atol) else 0.0
        
    return {
        "ok": int(mismatches) == 0,
        "max_abs": max_err,
        "mean_abs": mean_abs_err,
        "max_rel": np.max(rel_diff) if rel_diff.size > 0 else 0.0,
        "rmse": norm_diff / np.sqrt(g.size) if g.size > 0 else 0.0,
        "snr_db": snr_db,
        "nme_max": max_err / max_abs_w if max_abs_w > 0 else (float("inf") if max_err > 0 else 0.0),
        "nme_range": max_err / w_range if w_range > 0 else (float("inf") if max_err > 0 else 0.0),
        "bias": mean_signed_err,
        "bias_ratio": mean_signed_err / mean_abs_err if mean_abs_err > 0 else 0.0,
        "error_autocorr": e_autocorr,
        "pearson_corr": p_corr,
        "mismatches": int(mismatches),
        "total": g.size,
        "has_nan": not (np.all(np.isfinite(g)) and np.all(np.isfinite(w)))
    }

def comparison_worker(task: Tuple[Any, Any, int, int, Path, Path], var: str, metric: str, atol: float, rtol: float) -> Tuple[Any, Any, int, int, str]:
    """Worker function for parallel processing. Returns (comp_type, time_t, tstep, value)"""
    comp_type, time_t, tstep, _, curr_path, ref_path = task
    try:
        d_curr = extract_variable(curr_path, var)
        d_ref  = extract_variable(ref_path, var)
        if d_curr.size != d_ref.size:
            m_size = min(d_curr.size, d_ref.size)
            d_curr, d_ref = d_curr[:m_size], d_ref[:m_size]
        stats = compute_stats(d_curr, d_ref, atol, rtol)
        val = stats.get(metric, "err")
        return comp_type, time_t, tstep, f"{val:.4f}" if isinstance(val, float) else str(val)
    except Exception:
        return comp_type, time_t, tstep, "skip"

def print_summary_table(times, tsteps, results_map, title, metric, var):
    final_rows = []
    for t in times:
        row = {"T": t}
        row.update({str(ts): results_map[t][ts] for ts in tsteps})
        final_rows.append(row)

    df = pl.DataFrame(final_rows)
    
    # Filter out empty columns
    cols_to_keep = ["T"]
    for ts in tsteps:
        col_name = str(ts)
        if not (df[col_name] == "-").all():
            cols_to_keep.append(col_name)
    df = df.select(cols_to_keep)

    # Filter out empty rows
    if len(cols_to_keep) > 1:
        comp_cols = [c for c in cols_to_keep if c != "T" and not (df[c] == "ref").all()]
        if comp_cols:
            mask = pl.any_horizontal([pl.col(c) != "-" for c in comp_cols])
            df = df.filter(mask)

    print("\n" + "="*80)
    print(f"{title} (Metric: {metric}, Var: {var})")
    print("="*80)
    with pl.Config(tbl_rows=100, tbl_cols=-1, tbl_width_chars=1000):
        print(df)
    print("="*80)

def main():
    parser = argparse.ArgumentParser(description="Compare variables across refinement levels.")
    parser.add_argument("root", help="Root directory containing tstep folders.")
    parser.add_argument("-v", "--var", default="vn", help="Variable name (default: vn)")
    parser.add_argument("-m", "--metric", default="snr_db", help="Metric to display (snr_db, rmse, max_abs, etc.)")
    parser.add_argument("--atol", type=float, default=1e-12)
    parser.add_argument("--rtol", type=float, default=1e-12)
    parser.add_argument("-j", "--workers", type=int, default=DEFAULT_WORKERS)
    args = parser.parse_args()
    
    root = Path(args.root)
    if not root.is_dir():
        print(f"Error: {root} is not a directory.")
        sys.exit(1)

    # Print metric description
    if args.metric not in METRIC_DESCRIPTIONS:
        print(f"Error: Unknown metric '{args.metric}'.")
        print(f"Available metrics: {', '.join(METRIC_DESCRIPTIONS.keys())}")
        sys.exit(1)
        
    desc = METRIC_DESCRIPTIONS[args.metric]
    print(f"\n>>> Selected Metric: {args.metric}")
    print(f">>> Description: {desc}\n")

    # 1. Discover files
    files = []
    tstep_pat = re.compile(r"tstep\.(\d+)")
    file_pat = re.compile(r"p_prog\.t0\.(\d+)\.data")

    print(f"Scanning {root} for data files...")
    all_dirs = list(root.iterdir())
    for tstep_dir in tqdm(all_dirs, desc="Directories"):
        ts_m = tstep_pat.match(tstep_dir.name)
        if not ts_m or not tstep_dir.is_dir(): continue
        tstep_val = int(ts_m.group(1))
        for data_file in tstep_dir.glob("p_prog.t0.*.data"):
            f_m = file_pat.match(data_file.name)
            if not f_m: continue
            call_num = int(f_m.group(1))
            clock_t = ((call_num - 1) // 5) * tstep_val
            files.append({"tstep": tstep_val, "T": clock_t, "path": data_file})

    if not files:
        print("No matching files found.")
        return

    df_files = pl.DataFrame(files)
    tsteps = sorted(df_files["tstep"].unique().to_list())
    times = sorted(df_files["T"].unique().to_list())

    # 2. Prepare tasks for BOTH comparison types
    tasks = []
    res_consecutive = {t: {ts: "-" for ts in tsteps} for t in times}
    res_baseline = {t: {ts: "-" for ts in tsteps} for t in times}
    
    for t in times:
        row_files = df_files.filter(pl.col("T") == t).sort("tstep")
        if row_files.is_empty(): continue
        
        # Absolute baseline for this T
        abs_baseline_f = row_files[0]
        abs_baseline_ts = abs_baseline_f["tstep"][0]
        abs_baseline_path = abs_baseline_f["path"][0]
        
        for i, ts in enumerate(tsteps):
            curr_f = row_files.filter(pl.col("tstep") == ts)
            if curr_f.is_empty(): continue
            curr_path = curr_f["path"][0]
            
            # Baseline Comparison
            if ts == abs_baseline_ts:
                res_baseline[t][ts] = "ref"
            else:
                tasks.append(("baseline", t, ts, abs_baseline_ts, curr_path, abs_baseline_path))
                
            # Consecutive Comparison
            if i == 0:
                res_consecutive[t][ts] = "-"
            else:
                prev_ts = tsteps[i-1]
                prev_f = row_files.filter(pl.col("tstep") == prev_ts)
                if not prev_f.is_empty():
                    tasks.append(("consecutive", t, ts, prev_ts, curr_path, prev_f["path"][0]))

    # 3. Parallel processing
    print(f"Performing {len(tasks)} comparisons using {args.workers} workers...")
    worker_fn = partial(comparison_worker, var=args.var, metric=args.metric, atol=args.atol, rtol=args.rtol)
    
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        for c_type, t, ts, result in tqdm(executor.map(worker_fn, tasks), total=len(tasks), desc="Comparisons"):
            if c_type == "baseline":
                res_baseline[t][ts] = result
            else:
                res_consecutive[t][ts] = result

    # 4. Final display
    print_summary_table(times, tsteps, res_baseline, "Absolute Error (vs most refined)", args.metric, args.var)
    print_summary_table(times, tsteps, res_consecutive, "Relative Error (vs next finer level)", args.metric, args.var)

if __name__ == "__main__":
    main()

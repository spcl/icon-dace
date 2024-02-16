import numpy as np

MODE_INIT = "init"
MODE_PERTURB = "perturb"
MODE_CHECK = "check"
MODE_STATS = "stats"
MODE_TOLERANCE = "tolerance"
MODE_RUN_ENSEMBLE = "run"
MODE_VISUALIZE = "visualize"
MODE_PERFORMANCE = "performance"
MODE_CDO = "cdo"
CHECK_THRESHOLD = 1e-15

dataframe_type_dict = {
    "min": np.float64,
    "mean": np.float64,
    "max": np.float64,
    "ntime": np.int32,
    "time": np.float64,
    "height": np.float64,
    "name": str,
}

datetime_format = "%Y-%m-%dT%H:%M:%S"

compute_statistics = ["mean", "max", "min"]  # correspond to xarray.DataSet methods

cdo_bins = [1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1]

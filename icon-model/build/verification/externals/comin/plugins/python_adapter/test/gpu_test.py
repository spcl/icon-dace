import comin
import sys
import numpy as np

print("Hallo Welt!", file=sys.stderr)
glob = comin.descrdata_get_global()
print(f"{glob.device_name=}", file=sys.stderr)
print(f"{glob.device_vendor=}", file=sys.stderr)
print(f"{glob.device_driver=}", file=sys.stderr)

if "NVIDIA" in glob.device_vendor.upper():
    try:
        print("Using cupy!", file=sys.stderr)
        import cupy as xp
    except ImportError:
        print("Cannot import cupy, falling back to numpy", file=sys.stderr)
        import numpy as xp
else:
    print("No NVIDIA device found falling back to numpy", file=sys.stderr)
    import numpy as xp


@comin.register_callback(comin.EP_SECONDARY_CONSTRUCTOR)
def sec_ctor():
    global ta, ta_device, ta_host
    ta = comin.var_get([comin.EP_ATM_WRITE_OUTPUT_BEFORE], ("temp", 1), comin.COMIN_FLAG_READ)
    ta_device = comin.var_get([comin.EP_ATM_PHYSICS_AFTER], ("temp", 1), comin.COMIN_FLAG_WRITE | comin.COMIN_FLAG_DEVICE)
    ta_host = comin.var_get([comin.EP_ATM_NUDGING_BEFORE], ("temp", 1), comin.COMIN_FLAG_READ)


@comin.register_callback(comin.EP_ATM_WRITE_OUTPUT_BEFORE)
def foo():
    print(f"{ta.__cuda_array_interface__=}", file=sys.stderr)
    ta_arr = np.asarray(ta)
    if hasattr(ta_arr, "__cuda_array_interface__"):
        print(f"{ta_arr.__cuda_array_interface__=}", file=sys.stderr)
    print(f"{type(ta_arr)=}", file=sys.stderr)
    if hasattr(ta_arr, "device"):
        print(f"{ta_arr.device=}", file=sys.stderr)
    print(f"{ta_arr.base}", file=sys.stderr)
    print(f"{ta_arr.mean()=}", file=sys.stderr)
    tas = ta_arr[:, -1, :, 0, 0]
    print(f"{tas.mean()=}", file=sys.stderr)


@comin.register_callback(comin.EP_ATM_PHYSICS_AFTER)
def set_to_42():
    ta_xp = xp.asarray(ta_device)
    ta_xp[:] = 42.


@comin.register_callback(comin.EP_ATM_NUDGING_BEFORE)
def print_element():
    ta_np = np.asarray(ta_host)
    assert np.allclose(ta_np, 42.)
    print("check successful", file=sys.stderr)

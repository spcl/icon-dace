try:
    from ._yac import *
except Exception as e:
    print("Your PYTHONPATH probably points to the source directory of yac/python instead of the build directory.")
    raise

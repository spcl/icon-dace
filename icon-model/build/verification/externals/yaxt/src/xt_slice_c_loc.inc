#ifdef HAVE_FC_C_LOC_OF_SLICE
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_loc
#define XT_SLICE_C_LOC(slice, cptr) cptr = C_LOC(slice)
#else
  USE xt_core, ONLY: xt_slice_c_loc
#define XT_SLICE_C_LOC(slice, cptr) CALL xt_slice_c_loc(slice, cptr)
#endif

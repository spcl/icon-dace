This example shows how to offload computation to a external process
using yaxt. For that the `decomp_domain` and `glb_index` information
provided in the descrdata is used to construct the xmap for yaxt. The
external process is currently not parallelized, i.e. it uses the
indices 1..ncells.

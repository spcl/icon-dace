# ComIn Python API


## Global ComIn functions, variables and constants
- `COMIN_FLAG_DEVICE` (int)
- `COMIN_FLAG_NONE` (int)
- `COMIN_FLAG_READ` (int)
- `COMIN_FLAG_WRITE` (int)
- `COMIN_HGRID_UNSTRUCTURED_CELL` (int)
- `COMIN_HGRID_UNSTRUCTURED_EDGE` (int)
- `COMIN_HGRID_UNSTRUCTURED_VERTEX` (int)
- `COMIN_ZAXIS_2D` (int)
- `COMIN_ZAXIS_3D` (int)
- `COMIN_ZAXIS_3D_HALF` (int)
- `COMIN_ZAXIS_NONE` (int)
- `COMIN_ZAXIS_UNDEF` (int)
- `EP_ATM_ADVECTION_AFTER` (int)
- `EP_ATM_ADVECTION_BEFORE` (int)
- `EP_ATM_CHECKPOINT_AFTER` (int)
- `EP_ATM_CHECKPOINT_BEFORE` (int)
- `EP_ATM_CONVECTION_AFTER` (int)
- `EP_ATM_CONVECTION_BEFORE` (int)
- `EP_ATM_GWDRAG_AFTER` (int)
- `EP_ATM_GWDRAG_BEFORE` (int)
- `EP_ATM_INIT_FINALIZE` (int)
- `EP_ATM_INTEGRATE_AFTER` (int)
- `EP_ATM_INTEGRATE_BEFORE` (int)
- `EP_ATM_INTEGRATE_END` (int)
- `EP_ATM_INTEGRATE_START` (int)
- `EP_ATM_MICROPHYSICS_AFTER` (int)
- `EP_ATM_MICROPHYSICS_BEFORE` (int)
- `EP_ATM_NUDGING_AFTER` (int)
- `EP_ATM_NUDGING_BEFORE` (int)
- `EP_ATM_PHYSICS_AFTER` (int)
- `EP_ATM_PHYSICS_BEFORE` (int)
- `EP_ATM_RADHEAT_AFTER` (int)
- `EP_ATM_RADHEAT_BEFORE` (int)
- `EP_ATM_RADIATION_AFTER` (int)
- `EP_ATM_RADIATION_BEFORE` (int)
- `EP_ATM_SURFACE_AFTER` (int)
- `EP_ATM_SURFACE_BEFORE` (int)
- `EP_ATM_TIMELOOP_AFTER` (int)
- `EP_ATM_TIMELOOP_BEFORE` (int)
- `EP_ATM_TIMELOOP_END` (int)
- `EP_ATM_TIMELOOP_START` (int)
- `EP_ATM_TURBULENCE_AFTER` (int)
- `EP_ATM_TURBULENCE_BEFORE` (int)
- `EP_ATM_WRITE_OUTPUT_AFTER` (int)
- `EP_ATM_WRITE_OUTPUT_BEFORE` (int)
- `EP_ATM_YAC_DEFCOMP_AFTER` (int)
- `EP_ATM_YAC_DEFCOMP_BEFORE` (int)
- `EP_ATM_YAC_ENDDEF_AFTER` (int)
- `EP_ATM_YAC_ENDDEF_BEFORE` (int)
- `EP_ATM_YAC_SYNCDEF_AFTER` (int)
- `EP_ATM_YAC_SYNCDEF_BEFORE` (int)
- `EP_DESTRUCTOR` (int)
- `EP_FINISH` (int)
- `EP_SECONDARY_CONSTRUCTOR` (int)
- callable `callback_get_ep_name`: C function signature: void comin_callback_get_ep_name(int iep, char out_ep_name[MAX_LEN_EP_NAME+1])
- callable `current_get_datetime`: C function signature: `void comin_current_get_datetime(char const**,int*,int*);`
- callable `current_get_domain_id`: C function signature: `int comin_current_get_domain_id()`
- callable `current_get_plugin_info`: returns object describing the current plugin
- callable `descrdata_get_block`: C function signature: `int comin_descrdata_get_block(int j)`
- callable `descrdata_get_cell_indices`: C function signature: `void comin_descrdata_get_cell_indices(int jg, int i_blk, int i_startblk, int i_endblk, int* i_startidx, int* i_endidx, int irl_start, int irl_end)`
- callable `descrdata_get_cell_npromz`: C function signature: `int comin_descrdata_get_cell_npromz(int jg)`
- callable `descrdata_get_domain`: returns descriptive data for a given domain, arguments: jg
- callable `descrdata_get_edge_npromz`: C function signature: `int comin_descrdata_get_edge_npromz(int jg)`
- callable `descrdata_get_global`: returns global descriptive data object
- callable `descrdata_get_index`: C function signature: `int comin_descrdata_get_index(int j)`
- callable `descrdata_get_simulation_interval`: "returns simulation intervals: exp_start, exp_stop, run_start, run_stop
- callable `descrdata_get_timesteplength`: C function signature: `void double comin_descrdata_get_timesteplength(int jg)`
- callable `descrdata_get_vert_npromz`: C function signature: `int comin_descrdata_get_vert_npromz(int jg)`
- callable `descrdata_index_lookup_glb2loc_cell`: C function signature: `int comin_descrdata_index_lookup_glb2loc_cell(int jg, int global_idx)`
- callable `finish`: C function signature: `void comin_plugin_finish(const char* routine, const char* text)`
- callable `metadata_get`: retrieve metadata, arguments: (name string, domain id) , metadata key
- callable `metadata_set`: sets metadata for a requested field, arguments: name string, domain id, metadata key, metadata value
- callable `parallel_get_host_mpi_comm`: C function signature: `int comin_parallel_get_host_mpi_comm()`
- callable `parallel_get_host_mpi_rank`: C function signature: `int comin_parallel_get_host_mpi_rank()`
- callable `parallel_get_plugin_mpi_comm`: C function signature: `int comin_parallel_get_host_mpi_comm()`
- callable `setup_get_verbosity_level`: C function signature: `int comin_setup_get_verbosity_level()`
- callable `setup_get_version`: returns (major, minor, patch) version info
ismodule sys
- callable `var_descr_list`: List of exposed variables (descriptors)
- callable `var_get`: get variable object, arguments: [entry point], (name string, domain id), access flag)
- callable `var_request_add`: Request the host model to add a variable, arguments: (name string, domain id), lmodexclusive


## Members of data type plugin_info
- `plugin_info.args` (list)
- `plugin_info.comm` (str)
- `plugin_info.id` (int)
- `plugin_info.options` (str)


## Global descriptive data
- `global.device_driver` (str)
- `global.device_name` (str)
- `global.device_vendor` (str)
- `global.grf_bdywidth_c` (int)
- `global.grf_bdywidth_e` (int)
- `global.has_device` (bool)
- `global.host_git_branch` (str)
- `global.host_git_remote_url` (str)
- `global.host_git_tag` (str)
- `global.host_revision` (str)
- `global.lrestartrun` (bool)
- `global.max_dom` (int)
- `global.max_rlcell` (int)
- `global.max_rledge` (int)
- `global.max_rlvert` (int)
- `global.min_rlcell` (int)
- `global.min_rlcell_int` (int)
- `global.min_rledge` (int)
- `global.min_rledge_int` (int)
- `global.min_rlvert` (int)
- `global.min_rlvert_int` (int)
- `global.n_dom` (int)
- `global.nproma` (int)
- `global.vct_a(:)` (float64)
- `global.wp` (int)
- `global.yac_instance_id` (int)


## Members of data type simulation_interval
- `simulation_interval.exp_start` (str)
- `simulation_interval.exp_stop` (str)
- `simulation_interval.run_start` (str)
- `simulation_interval.run_stop` (str)


## Descriptive data for domains
- `domain.cells` (_descrdata)
- `domain.dom_end` (float)
- `domain.dom_start` (float)
- `domain.edges` (_descrdata)
- `domain.grid_filename` (str)
- `domain.grid_uuid` (memoryview)
- `domain.id` (int)
- `domain.n_childdom` (int)
- `domain.nlev` (int)
- `domain.nshift` (int)
- `domain.nshift_total` (int)
- `domain.number_of_grid_used(:)` (int32)
- `domain.verts` (_descrdata)


## Members of data type domain.cells
- `domain.cells.area(:,:)` (float64)
- `domain.cells.child_blk(:,:,:)` (int32)
- `domain.cells.child_id(:,:)` (int32)
- `domain.cells.child_idx(:,:,:)` (int32)
- `domain.cells.clat(:,:)` (float64)
- `domain.cells.clon(:,:)` (float64)
- `domain.cells.decomp_domain(:,:)` (int32)
- `domain.cells.edge_blk(:,:,:)` (int32)
- `domain.cells.edge_idx(:,:,:)` (int32)
- `domain.cells.end_block(:)` (int32)
- `domain.cells.end_index(:)` (int32)
- `domain.cells.glb_index(:)` (int32)
- `domain.cells.hhl(:,:,:)` (float64)
- `domain.cells.max_connectivity` (int)
- `domain.cells.nblks` (int)
- `domain.cells.ncells` (int)
- `domain.cells.ncells_global` (int)
- `domain.cells.neighbor_blk(:,:,:)` (int32)
- `domain.cells.neighbor_idx(:,:,:)` (int32)
- `domain.cells.num_edges(:,:)` (int32)
- `domain.cells.parent_glb_blk(:,:)` (int32)
- `domain.cells.parent_glb_idx(:,:)` (int32)
- `domain.cells.refin_ctrl(:,:)` (int32)
- `domain.cells.start_block(:)` (int32)
- `domain.cells.start_index(:)` (int32)
- `domain.cells.vertex_blk(:,:,:)` (int32)
- `domain.cells.vertex_idx(:,:,:)` (int32)


## Members of data type domain.edges
- `domain.edges.cell_blk(:,:,:)` (int32)
- `domain.edges.cell_idx(:,:,:)` (int32)
- `domain.edges.child_blk(:,:,:)` (int32)
- `domain.edges.child_id(:,:)` (int32)
- `domain.edges.child_idx(:,:,:)` (int32)
- `domain.edges.elat(:,:)` (float64)
- `domain.edges.elon(:,:)` (float64)
- `domain.edges.end_block(:)` (int32)
- `domain.edges.end_index(:)` (int32)
- `domain.edges.nblks` (int)
- `domain.edges.nedges` (int)
- `domain.edges.nedges_global` (int)
- `domain.edges.parent_glb_blk(:,:)` (int32)
- `domain.edges.parent_glb_idx(:,:)` (int32)
- `domain.edges.refin_ctrl(:,:)` (int32)
- `domain.edges.start_block(:)` (int32)
- `domain.edges.start_index(:)` (int32)
- `domain.edges.vertex_blk(:,:,:)` (int32)
- `domain.edges.vertex_idx(:,:,:)` (int32)


## Members of data type domain.verts
- `domain.verts.cell_blk(:,:,:)` (int32)
- `domain.verts.cell_idx(:,:,:)` (int32)
- `domain.verts.edge_blk(:,:,:)` (int32)
- `domain.verts.edge_idx(:,:,:)` (int32)
- `domain.verts.end_block(:)` (int32)
- `domain.verts.end_index(:)` (int32)
- `domain.verts.nblks` (int)
- `domain.verts.neighbor_blk(:,:,:)` (int32)
- `domain.verts.neighbor_idx(:,:,:)` (int32)
- `domain.verts.num_edges(:,:)` (int32)
- `domain.verts.nverts` (int)
- `domain.verts.nverts_global` (int)
- `domain.verts.refin_ctrl(:,:)` (int32)
- `domain.verts.start_block(:)` (int32)
- `domain.verts.start_index(:)` (int32)
- `domain.verts.vlat(:,:)` (float64)
- `domain.verts.vlon(:,:)` (float64)

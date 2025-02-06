#!/usr/bin/env python3

from yac import YAC, Reg2dGrid, Location, Field, TimeUnit, Reduction, InterpolationStack, NNNReductionType, UnstructuredGrid
import comin
import numpy as np
import matplotlib.pyplot as plt

glob = comin.descrdata_get_global()
assert glob.yac_instance_id != -1, "The host-model is not configured with yac"

yac = YAC.from_id(glob.yac_instance_id)
source_comp = yac.predef_comp("comin_example_source")

domain = comin.descrdata_get_domain(1)
connectivity = (np.asarray(domain.cells.vertex_blk)-1)*glob.nproma + (np.asarray(domain.cells.vertex_idx)-1)
icon_grid = UnstructuredGrid("comin_icon_grid",
                             np.ones(domain.cells.ncells)*3,
                             np.array(np.ravel(np.transpose(domain.verts.vlon))[:domain.verts.nverts]),
                             np.array(np.ravel(np.transpose(domain.verts.vlat))[:domain.verts.nverts]),
                             np.ravel(np.swapaxes(connectivity, 0, 1))[:3*domain.cells.ncells])
icon_cell_centers = icon_grid.def_points(Location.CELL,
                                         np.ravel(domain.cells.clon)[:domain.cells.ncells],
                                         np.ravel(domain.cells.clat)[:domain.cells.ncells])

rank = comin.parallel_get_host_mpi_rank()
if rank == 0:
    target_comp = yac.predef_comp("comin_example_target")

    vlon = np.linspace(-np.pi, np.pi, 360)
    vlat = np.linspace(np.pi/2, -np.pi/2, 180)
    grid = Reg2dGrid("comin_example_grid",
                     vlon, vlat,
                     cyclic=(True, False))
    cell_centers = grid.def_points(Location.CELL,
                                   vlon + 0.5*(vlon[1]-vlon[0]),
                                   vlat[:-1] + 0.5*(vlat[1]-vlat[0]))


@comin.register_callback(comin.EP_SECONDARY_CONSTRUCTOR)
def sec_ctr():
    global comin_pres_sfc
    comin_pres_sfc = comin.var_get([comin.EP_ATM_WRITE_OUTPUT_BEFORE], ("pres_sfc", 1),
                                   comin.COMIN_FLAG_READ)


@comin.register_callback(comin.EP_ATM_YAC_DEFCOMP_AFTER)
def define_fields():
    global yac_pres_sfc_target, yac_pres_sfc_source
    yac_pres_sfc_source = Field.create("pres_sfc", source_comp, icon_cell_centers, 1,
                                       str(int(comin.descrdata_get_timesteplength(1))),
                                       TimeUnit.SECOND)
    if rank == 0:
        yac_pres_sfc_target = Field.create("pres_sfc", target_comp, cell_centers, 1,
                                           str(int(comin.descrdata_get_timesteplength(1))),
                                           TimeUnit.SECOND)
        nnn = InterpolationStack()
        nnn.add_nnn(NNNReductionType.AVG, 1, 1.0)
        yac.def_couple("comin_example_source", "comin_icon_grid", "pres_sfc",
                       "comin_example_target", "comin_example_grid", "pres_sfc",
                       str(int(comin.descrdata_get_timesteplength(1))),
                       TimeUnit.SECOND, Reduction.TIME_NONE, nnn)


@comin.register_callback(comin.EP_ATM_TIMELOOP_END)
def process():
    yac_pres_sfc_source.put(np.ravel(comin_pres_sfc)[:domain.cells.ncells])
    if rank == 0:
        data, info = yac_pres_sfc_target.get()
        plt.imshow(np.reshape(data, [179, 360]))
        print(f"pres_sfc_{comin.current_get_datetime()}.png")
        plt.savefig(f"pres_sfc_{comin.current_get_datetime()}.png")

import dace

sdfg = dace.SDFG.from_file("velocity_tendencies_simplified_f.sdfgz")

l = ["cz_c", "cn_e", "deepatmo_t2mc", "deepatmo_t1ifc",
     "dzgpot_mc", "zgpot_mc", "zgpot_ifc", "deepatmo_t2mc", "deepatmo_t1mc",
     "ddt_vn_cen", "vor_v", "vor_q", "ddt_ua_cen_is_associated",
     "ddt_va_cen_is_associated", "ddt_va_cen", "ddt_ua_cen", "fbk_dom_volume",
     "vor_u"]

def find2(name, dt : dace.data.Structure):
    items = dt.members.items()
    to_delete = []
    for a_name, arr in items:
        if name in a_name:
            print(a_name, arr.shape, arr.dtype)
            to_delete.append(a_name)
        if isinstance(arr, dace.data.Structure):
            find2(name, arr)
    for a_name in to_delete:
        del dt.members[a_name]

def find(name):
    items = list(sdfg.arrays.items())
    to_delete = []
    for a_name, arr in items:
        if name in a_name:
            print(a_name, arr.shape, arr.dtype)
            to_delete.append(a_name)
        if isinstance(arr, dace.data.Structure):
            find2(name, arr)
    for a_name in to_delete:
        del sdfg.arrays[a_name]

    items = list(sdfg.symbols.items())
    for sym_name, sym in items:
        if name in sym_name and "f2dace" in sym_name:
            print(sym_name)

    for s in sdfg.states():
        for n in s.nodes():
            if isinstance(n, dace.nodes.NestedSDFG):
                find(name)

for name in l:
    find(name)

sdfg.validate()

sdfg.save("velocity_tendencies_simplified_f_pruned.sdfgz")
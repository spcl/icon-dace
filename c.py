import dace

sdfg = dace.SDFG.from_file("radiation_simplified_dbg22.sdfgz")
sdfg.compile()
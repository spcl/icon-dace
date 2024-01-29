#!/usr/bin/env python3

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

import faulthandler, os

# enable functions to dump Python tracebacks explicitly:
faulthandler.enable()

# limit number of OpenMP threads;
#this fixes a NumPy SegFault related to the OpenBLAS threadpool
# (see https://github.com/numpy/numpy/issues/15691):
os.environ["OMP_NUM_THREADS"] = "1"

print("""
==================================================================================
COMPARE TEST DATA TO REFERENCE FILE
==================================================================================

# Perform a number of comparisons between two GRIB data files.
# 03/20 : F. Prill (DWD)
""")
# --------------------------------------------------------------------------------


# --------------------------------------------------------------------------------
# > define a list of checks and their conditions for execution,
#   where the conditions may combine several metadata.
async def run_checks_glb(test_data, reference):
    await asyncio.gather(
        # relative difference of average:
        check_rel_msd(test_data, reference, 0.001,  grb_metadata(test_data.grb, "shortName") == "PMSL"),    \
        check_rel_avg(test_data, reference, 0.00004, grb_metadata(test_data.grb, "shortName") == "T_G"),     \
        check_rel_avg(test_data, reference, 0.0001, grb_metadata(test_data.grb, "shortName") == "T"),       \
        check_rel_avg(test_data, reference, 0.00003, grb_metadata(test_data.grb, "shortName") == "T_2M"),    \
        check_rel_avg(test_data, reference, 0.0002, grb_metadata(test_data.grb, "shortName") == "RELHUM_2M"),    \
        check_rel_avg(test_data, reference, 0.0015, grb_metadata(test_data.grb, "shortName") == "TOT_PREC"),     \
        check_rel_avg(test_data, reference, 0.0002, grb_metadata(test_data.grb, "shortName") == "TQV"),     \
        check_rel_avg(test_data, reference, 0.008,  grb_metadata(test_data.grb, "shortName") == "TQC_DIA"),  \
        check_rel_avg(test_data, reference, 0.005, grb_metadata(test_data.grb, "shortName") == "TQI_DIA"),  \
        check_rel_avg(test_data, reference, 0.005, grb_metadata(test_data.grb, "shortName") == "CLCL"),     \
        check_rel_avg(test_data, reference, 0.005, grb_metadata(test_data.grb, "shortName") == "CLCM"),     \
        check_rel_avg(test_data, reference, 0.005, grb_metadata(test_data.grb, "shortName") == "CLCH"),     \
        check_rel_avg(test_data, reference, 0.002, grb_metadata(test_data.grb, "shortName") == "CLCT"),     \
        check_rel_avg(test_data, reference, 0.0003, grb_metadata(test_data.grb, "shortName") == "ASOB_S"),   \
        check_rel_avg(test_data, reference, 0.0003, grb_metadata(test_data.grb, "shortName") == "ASOB_T"),   \
        check_rel_avg(test_data, reference, 0.0001, grb_metadata(test_data.grb, "shortName") == "ATHB_S"),   \
        check_rel_avg(test_data, reference, 0.0001, grb_metadata(test_data.grb, "shortName") == "ATHB_T"),   \
        check_rel_avg(test_data, reference, 0.0005, grb_metadata(test_data.grb, "shortName") == "ALHFL_S"),  \
        check_rel_avg(test_data, reference, 0.0005, grb_metadata(test_data.grb, "shortName") == "ASHFL_S"),  \
        check_rel_avg(test_data, reference, 0.001, grb_metadata(test_data.grb, "shortName") == "W_SO"),     \
        check_rel_avg(test_data, reference, 0.002, grb_metadata(test_data.grb, "shortName") == "RAIN_GSP"), \
        check_rel_avg(test_data, reference, 0.003, grb_metadata(test_data.grb, "shortName") == "RAIN_CON"), \
        check_rel_avg(test_data, reference, 0.003, grb_metadata(test_data.grb, "shortName") == "SNOW_GSP"), \
        check_rel_avg(test_data, reference, 0.005, grb_metadata(test_data.grb, "shortName") == "SNOW_CON"), \
        # GRIB2 meta-data checks:
        check_grb_metadata(test_data, reference, True),                                                  \
    )

async def run_checks_loc(test_data, reference):
    await asyncio.gather(
        # relative difference of average:
        check_rel_msd(test_data, reference, 0.001,  grb_metadata(test_data.grb, "shortName") == "PMSL"),    \
        check_rel_avg(test_data, reference, 0.00006, grb_metadata(test_data.grb, "shortName") == "T_G"),     \
        check_rel_avg(test_data, reference, 0.0001, grb_metadata(test_data.grb, "shortName") == "T"),       \
        check_rel_avg(test_data, reference, 0.00005, grb_metadata(test_data.grb, "shortName") == "T_2M"),    \
        check_rel_avg(test_data, reference, 0.00055, grb_metadata(test_data.grb, "shortName") == "RELHUM_2M"),    \
        check_rel_avg(test_data, reference, 0.002,  grb_metadata(test_data.grb, "shortName") == "TOT_PREC"),     \
        check_rel_avg(test_data, reference, 0.0005, grb_metadata(test_data.grb, "shortName") == "TQV"),     \
        check_rel_avg(test_data, reference, 0.02,  grb_metadata(test_data.grb, "shortName") == "TQC_DIA"),  \
        check_rel_avg(test_data, reference, 0.015,  grb_metadata(test_data.grb, "shortName") == "TQI_DIA"),  \
        check_rel_avg(test_data, reference, 0.015,  grb_metadata(test_data.grb, "shortName") == "CLCL"),     \
        check_rel_avg(test_data, reference, 0.015,  grb_metadata(test_data.grb, "shortName") == "CLCM"),     \
        check_rel_avg(test_data, reference, 0.015,  grb_metadata(test_data.grb, "shortName") == "CLCH"),     \
        check_rel_avg(test_data, reference, 0.008,  grb_metadata(test_data.grb, "shortName") == "CLCT"),     \
        check_rel_avg(test_data, reference, 0.001, grb_metadata(test_data.grb, "shortName") == "ASOB_S"),   \
        check_rel_avg(test_data, reference, 0.001, grb_metadata(test_data.grb, "shortName") == "ASOB_T"),   \
        check_rel_avg(test_data, reference, 0.0005, grb_metadata(test_data.grb, "shortName") == "ATHB_S"),   \
        check_rel_avg(test_data, reference, 0.0005, grb_metadata(test_data.grb, "shortName") == "ATHB_T"),   \
        check_rel_avg(test_data, reference, 0.001, grb_metadata(test_data.grb, "shortName") == "ALHFL_S"),  \
        check_rel_avg(test_data, reference, 0.001, grb_metadata(test_data.grb, "shortName") == "ASHFL_S"),  \
        check_rel_avg(test_data, reference, 0.002, grb_metadata(test_data.grb, "shortName") == "W_SO"),     \
        check_rel_avg(test_data, reference, 0.0055, grb_metadata(test_data.grb, "shortName") == "RAIN_GSP"), \
        check_rel_avg(test_data, reference, 0.005, grb_metadata(test_data.grb, "shortName") == "RAIN_CON"), \
        check_rel_avg(test_data, reference, 0.01, grb_metadata(test_data.grb, "shortName") == "SNOW_GSP"), \
        check_rel_avg(test_data, reference, 0.01, grb_metadata(test_data.grb, "shortName") == "SNOW_CON"), \
        # GRIB2 meta-data checks:
        check_grb_metadata(test_data, reference, True),                                                  \
    )

async def run_checks_0K(test_data, reference):
    await asyncio.gather(
        check_gt_0(test_data, grb_metadata(test_data.grb, "shortName") == "T"),
        check_gt_0(test_data, grb_metadata(test_data.grb, "shortName") == "T_2M"),
        check_gt_0(test_data, grb_metadata(test_data.grb, "shortName") == "T_G"),
    )

# --------------------------------------------------------------------------------
# > TEST: check relative average (AVG(TEST) - AVG(REF)) / AVG(REF)
async def check_rel_avg(test_data, reference, tol, condition):
    if (condition):
        await read_reference_record(test_data.grb, reference)
        avg_test = numpy.array(grb_values(test_data.grb)).mean()
        avg_ref  = numpy.array(grb_values(reference.grb)).mean()
        val      = numpy.abs( (avg_test - avg_ref) / avg_ref )
        if (val > tol):
            raise Exception(">>> '{}' : {} > {} tol".
                            format(grb_metadata(reference.grb,"shortName"), val, tol))
        else:
            print(">>> check_rel_avg '{}' : passed. {} <= {} tol".
                  format(grb_metadata(reference.grb,"shortName"), val, tol))

# --------------------------------------------------------------------------------
# > TEST: check relative mean squared difference (AVG(TEST - REF)**2) / AVG(REF)
async def check_rel_msd(test_data, reference, tol, condition):
    if (condition):
        await read_reference_record(test_data.grb, reference)
        avg_test = numpy.array((grb_values(test_data.grb)-grb_values(reference.grb))**2).mean()
        avg_ref  = numpy.array(grb_values(reference.grb)).mean()
        val      = numpy.abs(avg_test/avg_ref)
        if (val > tol):
            raise Exception(">>> '{}' : {} > {} tol".
                            format(grb_metadata(reference.grb,"shortName"), val, tol))
        else:
            print(">>> check_rel_msd '{}' : passed. {} <= {} tol".
                  format(grb_metadata(reference.grb,"shortName"), val, tol))

# --------------------------------------------------------------------------------
# > TEST: check GRIB2 key/values (with a given "black list"):
async def check_grb_metadata(test_data, reference, condition):
    if (condition):
        await read_reference_record(test_data.grb, reference)

        GRBKEY_IGNORE_LIST = ["localCreationDate*", "codedValues", "md5Section*", "md5Headers", \
                              "localDecodeDate*", "maximum", "average", "standardDeviation",    \
                              "skewness", "kurtosis", "sectionNumber", "section2Padding", "referenceValue", \
                              "binaryScaleFactor", "totalLength", "bitsPerValue", "section7Length", \
                              "uuidOfVGrid"]

        ignore_combined = "(" + ")|(".join(GRBKEY_IGNORE_LIST) + ")" # combined regular expression
        iterid = eccodes.codes_keys_iterator_new(test_data.grb)
        eccodes.codes_skip_computed(iterid)
        differ = False
        while eccodes.codes_keys_iterator_next(iterid):
            keyname = eccodes.codes_keys_iterator_get_name(iterid)
            if not re.match(ignore_combined, keyname):
                keyval     = eccodes.codes_get_string(test_data.grb, keyname)
                keyval_ref = eccodes.codes_get_string(reference.grb, keyname)
                if (keyval != keyval_ref): 
                    print("metadata '{}' differs!".format(keyname))
                    print("  '{}' / '{}'".format(keyval, keyval_ref))
                    differ = True
        eccodes.codes_keys_iterator_delete(iterid)
        if differ:
            raise Exception("GRB metadata differs for '{}'!".format(grb_metadata(test_data.grb,"shortName")))
        else:
            print(">>> check_grb_metadata '{}' : passed.".format(grb_metadata(reference.grb,"shortName")))

# --------------------------------------------------------------------------------
# > TEST: check that variable is greater than 0
async def check_gt_0(test_data, condition):
    if (condition):
        min_test = numpy.array(grb_values(test_data.grb)).min()
        if (min_test <= 0.):
            raise Exception(">>> '{}' : min <= 0".
                            format(grb_metadata(test_data.grb,"shortName")))
        else:
            print(">>> check_rel_avg '{}' : passed. min <= 0".
                  format(grb_metadata(test_data.grb,"shortName")))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# --------------------------------------------------------------------------------
# DO NOT EDIT BELOW THIS LINE ----------------------------------------------------
# --------------------------------------------------------------------------------
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

import numpy, asyncio, argparse, sys, eccodes, re
grb_values   = eccodes.codes_get_values # alias
grb_metadata = eccodes.codes_get        # alias

# --------------------------------------------------------------------------------
# > read all GRIB2 records in data file; call check functions.
async def read_grib_record(test_data, reference, check_fct):

    while True:
        test_data.grb = eccodes.codes_grib_new_from_file(test_data.pfile)
        test_data.record += 1
        if (test_data.grb is None):
            break
        else:
            test_data.processed_records.append(grb_metadata(test_data.grb,"shortName"))
            await check_fct(test_data, reference)
            eccodes.codes_release(test_data.grb)
    return False

# --------------------------------------------------------------------------------
# > @return True if test and reference record correspond to the same data.
def grb_record_eqv(grbA,grbB):
    if (grbA is None) or (grbB is None):  return False
    for item in ["shortName", "scaleFactorOfFirstFixedSurface", "tileIndex", \
                 "scaledValueOfFirstFixedSurface", "dateTime", "partitionNumber"]:
        try:
            if (grb_metadata(grbA,item) != grb_metadata(grbB,item)):  return False
        except Exception:
            pass
    return True

# --------------------------------------------------------------------------------
# > iterate over reference file to find a specific record.
#   Since we cannot assume that the records in the reference file have
#   an identical ordering compared to the test data file, we need to
#   start from the beginning of the reference file when the first
#   search was without success (we assume that the reference file is
#   shorter or equal to the test data file).
async def read_reference_record(grb_tst, reference):
    for rep in [0,1]:
        while True:
            if (grb_record_eqv(grb_tst, reference.grb)):
                #print("read reference record '{}', level {}, record #{}".
                #      format(grb_metadata(reference.grb,"shortName"), \
                #             grb_metadata(reference.grb,"level"),     \
                #             reference.record))
                reference.processed_records.append(grb_metadata(reference.grb,"shortName"))
                return reference.grb
            if not (reference.grb is None): eccodes.codes_release(reference.grb)
            reference.grb = eccodes.codes_grib_new_from_file(reference.pfile)
            reference.record += 1
            if (reference.grb is None): 
                break
        if (rep == 0):
            reference.pfile.seek(0)
            reference.record = 0
    raise Exception("No reference found for record '{}'!".format(grb_metadata(grb_tst,"shortName")))
    return False

# --------------------------------------------------------------------------------
# > main program (wrapped by a subroutine)
def main():
    # parse command-line options
    parser = argparse.ArgumentParser()
    parser.add_argument("--datafile",          required=True,       help="GRIB2 file containing test data.")
    parser.add_argument("--reffile",                                help="GRIB2 file containing reference data.")
    parser.add_argument("--check_0K",          action="store_true", help="Check that T, T_2M and T_G are > 0 K.")
    parser.add_argument("--global_data", "-g", action="store_true", help="indicates if dataset is global.")

    args = parser.parse_args()

    print("test data file : {} \nreference file : {}\n\n".format(args.datafile, args.reffile))

    class Object(object):
        pass

    test_data = Object()
    test_data.pfile  = open(args.datafile)
    test_data.record = 0
    test_data.grb    = None
    test_data.processed_records = []

    if args.check_0K:
        try:
            print("Check that T, T_2M and T_G are > 0 K.")
            asyncio.run(read_grib_record(test_data, None, run_checks_0K))
        except Exception as e:
            print("\nERROR!\n")
            print("checked data records so far: {}".format(list(dict.fromkeys(test_data.processed_records))))
            print("\n\n{}\n\n".format(e))
            sys.exit(1)

        print("Done.")
        print("visited data records: {}".format(list(dict.fromkeys(test_data.processed_records))))
        sys.exit(0)

    if args.reffile is None:
        parser.error("the following arguments are required if check_0K is not set: --reffile")

    reference = Object()
    reference.pfile  = open(args.reffile)
    reference.record = 0
    reference.grb    = None
    reference.processed_records = []

    try:
        if (args.global_data):
            print("perform global checks.")
            asyncio.run(read_grib_record(test_data, reference, run_checks_glb))
        else:
            print("perform local checks.")
            asyncio.run(read_grib_record(test_data, reference, run_checks_loc))
    except Exception as e:
        print("\nERROR!\n")
        print("checked data records so far: {}".format(list(dict.fromkeys(test_data.processed_records))))
        print("\n\n{}\n\n".format(e))
        sys.exit(1)

    print("Done.")
    print("visited data records: {}".format(list(dict.fromkeys(test_data.processed_records))))
    print("processed records: {}".format(list(dict.fromkeys(reference.processed_records))))

if __name__ == '__main__':
    main()


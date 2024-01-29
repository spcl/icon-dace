#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ICON-Land
#
# ---------------------------------------
# Copyright (C) 2013-2024, MPI-M, MPI-BGC
#
# Contact: icon-model.org
# Authors: AUTHORS.md
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------

#
#  dsl4jsb.py
#
#  Fortran preprocessor for JSBACH DSL language
#
#  Call "dsl4jsb.py -h" for usage.
#  To test, call "dsl4jsb.py" without options in the source
#  directory (containing "dsl4jsb.py" and "test.f90")
#
#  Nov 10, 2016: Reiner Schnur, MPIM
#

import sys
import getopt
import re
import os
import fnmatch
import io
import hashlib

# Search pattern to find lines with dsl4jsb macros and extract macro name
re_function = re.compile(r'dsl4jsb_(.*?) *( *$|,|\)| |\()', re.IGNORECASE)

# Search pattern to find hash value of original file in processed file
re_hash = re.compile(r'with md5 hash:\s*(\S*)\s*', re.IGNORECASE)

# Table of macros
#   Each entry is of the form
#       macro name: [re.compile(match pattern, re.I), replacement pattern]
#   (macro name without the "dsl4jsb_" and in lower case)
#
macros = {
    'use_processes': [re.compile(
        r'dsl4jsb_Use_processes (.*)', re.IGNORECASE),
        r'USE mo_jsb_process_class, ONLY: \1'],
    'use_memory':
        [re.compile(r'dsl4jsb_Use_memory *\((.*)_ *\)', re.IGNORECASE),
         r'USE mo_\1_memory_class, ONLY: t_\1_memory'],
    'def_memory': [re.compile(
        r'dsl4jsb_Def_memory *\((.*)_ *\)', re.IGNORECASE),
        r'TYPE(t_\1_memory), POINTER :: \1__mem'],
    'def_memory_tile': [re.compile(
        r'dsl4jsb_Def_memory_tile *\((.*?)_ *, *(.*?) *\)', re.IGNORECASE),
        r'TYPE(t_\1_memory), POINTER :: \1__mem__\2'],
    'get_memory': [re.compile(
        r'dsl4jsb_Get_memory *\((.*)_ *\)', re.IGNORECASE),
        r'SELECT TYPE (mm => tile%mem(\1_)%p); TYPE IS (t_\1_memory); \1__mem => mm ; END SELECT'],
    'get_memory_tile': [re.compile(
        r'dsl4jsb_Get_memory_tile *\((.*?)_ *, *(.*?) *\)', re.IGNORECASE),
        r'SELECT TYPE (mm => \2%mem(\1_)%p); TYPE IS (t_\1_memory); \1__mem__\2 => mm ; END SELECT'],
    'memory': [re.compile(
        r'dsl4jsb_Memory *\((.*?)_ *\)', re.IGNORECASE),
        r'\1__mem'],
    'use_config': [re.compile(
        r'dsl4jsb_Use_config *\((.*)_ *\)', re.IGNORECASE),
        r'USE mo_\1_config_class, ONLY: t_\1_config'],
    'def_config': [re.compile(
        r'dsl4jsb_Def_config *\((.*)_ *\)', re.IGNORECASE),
        r'TYPE(t_\1_config), POINTER :: \1__conf'],
    'get_config': [re.compile(
        r'dsl4jsb_Get_config *\((.*)_ *\)', re.IGNORECASE),
        r'SELECT TYPE (cc => model%processes(\1_)%p%config); TYPE IS (t_\1_config) ; \1__conf => cc ; END SELECT'],
    'config': [re.compile(
        r'dsl4jsb_Config *\((.*?)_ *\)', re.IGNORECASE),
        r'\1__conf'],
    'lctlib': [re.compile(
        r'dsl4jsb_Lctlib', re.IGNORECASE),
        r'model%lctlib(tile%lcts(1)%lib_id)'],
    'lctlib_param': [re.compile(
        r'dsl4jsb_Lctlib_param *\((.*?)\)', re.IGNORECASE),
        r'model%lctlib(tile%lcts(1)%lib_id)%\1'],
    'real2d_onchunk': [re.compile(
        r'dsl4jsb_Real2D_onChunk', re.IGNORECASE),
        r'REAL(wp), POINTER, DIMENSION(:)'],
    'real3d_onchunk': [re.compile(
        r'dsl4jsb_Real3D_onChunk', re.IGNORECASE),
        r'REAL(wp), POINTER, DIMENSION(:,:)'],
    'get_var2d_onchunk': [re.compile(
        r'dsl4jsb_Get_var2D_onChunk *\((.*?)_ *, *(.*?) *\)', re.IGNORECASE),
        r'\2 => \1__mem%\2%ptr(ics:ice,iblk)'],
    'get_var2d_onchunk_tile' : [re.compile(
        r'dsl4jsb_Get_var2D_onChunk_tile *\((.*?)_ *, *(.*?) *, *(.*?) *\)', re.IGNORECASE),
        r'\2 => \1__mem__\3%\2%ptr(ics:ice,iblk)'],
    'get_var2d_onchunk_tile_name' : [re.compile(
        r'dsl4jsb_Get_var2D_onChunk_tile_name *\((.*?)_ *, *(.*?) *, *(.*?) *\)', re.IGNORECASE),
        r'\2_\3 => \1__mem__\3%\2%ptr(ics:ice,iblk)'],
    'get_var3d_onchunk': [re.compile(
        r'dsl4jsb_Get_var3D_onChunk *\((.*?)_ *, *(.*?) *\)', re.IGNORECASE),
        r'\2 => \1__mem%\2%ptr(ics:ice,:,iblk)'],
    'var2d_onchunk': [re.compile(
        r'dsl4jsb_var2D_onChunk *\((.*?)_ *, *(.*?) *\)', re.IGNORECASE),
        r'\1__mem%\2%ptr(ics:ice,iblk)'],
    'var3d_onchunk': [re.compile(
        r'dsl4jsb_var3D_onChunk *\((.*?)_ *, *(.*?) *\)', re.IGNORECASE),
        r'\1__mem%\2%ptr(ics:ice,:,iblk)'],
    'real2d_ondomain': [re.compile(
        r'dsl4jsb_Real2D_onDomain', re.IGNORECASE),
        r'REAL(wp), POINTER, DIMENSION(:,:)'],
    'real3d_ondomain': [re.compile(
        r'dsl4jsb_Real3D_onDomain', re.IGNORECASE),
        r'REAL(wp), POINTER, DIMENSION(:,:,:)'],
    'get_var2d_ondomain': [re.compile(
        r'dsl4jsb_Get_var2D_onDomain *\((.*?)_ *, *(.*?) *\)', re.IGNORECASE),
        r'\2 => \1__mem%\2%ptr(:,:)'],
    'get_var2d_ondomain_tile': [re.compile(
        r'dsl4jsb_Get_var2D_onDomain_tile *\((.*?)_ *, *(.*?) *, *(.*?) *\)', re.IGNORECASE),
        r'\2 => \1__mem__\3%\2%ptr(:,:)'],
    'get_var2d_ondomain_tile_name': [re.compile(
        r'dsl4jsb_Get_var2D_onDomain_tile_name *\((.*?)_ *, *(.*?) *, *(.*?) *\)', re.IGNORECASE),
        r'\2_\3 => \1__mem__\3%\2%ptr(:,:)'],
    'get_var3d_ondomain': [re.compile(
        r'dsl4jsb_Get_var3D_onDomain *\((.*?)_ *, *(.*?) *\)', re.IGNORECASE),
        r'\2 => \1__mem%\2%ptr(:,:,:)'],
    'var2d_ondomain': [re.compile(
        r'dsl4jsb_var2D_onDomain *\((.*?)_ *, *(.*?) *\)', re.IGNORECASE),
        r'\1__mem%\2%ptr(:,:)'],
    'var3d_ondomain': [re.compile(
        r'dsl4jsb_var3D_onDomain *\((.*?)_ *, *(.*?) *\)', re.IGNORECASE),
        r'\1__mem%\2%ptr(:,:,:)'],
    'var_ptr': [re.compile(
        r'dsl4jsb_var_ptr *\((.*?)_ *, *(.*?) *\)', re.IGNORECASE),
        r'\1__mem%\2%ptr'],
    'aggregate_onchunk': [re.compile(
        r'dsl4jsb_Aggregate_onChunk *\((.*?)_ *, *(.*?) *, *(.*?)\)', re.IGNORECASE),
        r'CALL \3%Aggregate(tile, \1__mem%\2,ics,ice,iblk,routine)']
    }

help_text = """
    dsl4jsb.py [-h] [OPTION...]

        -h                usage
        -v, --verbose     print list of conversions (use twice for message for files not processed)
        -i, --ifile=FILE  input file (test.f90 if missing and input directory (-d) not given)
        -o, --ofile=FILE  output file (stdout if missing)
        -d, --idir=DIR    directory with source files
        -t, --odir=DIR    target directory for output files (must exist)
        -k, --keep-dirs   keep directory tree structure of files in input directory (if -d and -t given)
        -p, --pp-add=STR  add string STR before .suffix in processed files (pp = pre-processed)
        -n, --no-header   don't create header with hash of source file
    """


def get_input_files(root):
    for dir, subdirs, files in os.walk(root):
        for file in fnmatch.filter(files, '*.f90'):
            if dir == '.':
                yield(file)
            else:
                yield os.path.join(dir.replace('./', ''), file)


def tolower(match):
    return match.group(1).lower() + '_'


def process_file(inputfile, outputfile, header=True, verbose=False, verbose_more=False, cwd=None):
    # Compute hash value of input file
    hasher = hashlib.md5()
    with io.open(inputfile, 'rb') as infile:
        buf = infile.read()
        hasher.update(buf)
    # If output file exists, extract hash value of file
    # this one was generated from
    hash = ''
    if outputfile != '' and os.path.isfile(outputfile):
        with io.open(outputfile, 'r', encoding="utf-8") as outfile:
            for line in outfile:
                matches = re_hash.findall(line)
                if matches:
                    hash = matches[0]
                    # print(hash)
                    break
    # If both hashes are equal, don't generate new file
    if hasher.hexdigest() == hash:
        # Since input file hasn't changed, reset its modification time to the one
        # of the output file, so that make doesn't think the input file is newer
        os.utime(inputfile, (os.path.getatime(inputfile), os.path.getmtime(outputfile)))
        if verbose_more:
            print('File hasn"t changed: ' + os.path.realpath(os.path.abspath(inputfile)))
        return

    with io.open(inputfile, 'r', encoding="utf-8") as infile:
        with (io.open(outputfile, 'w', encoding="utf-8") if outputfile != '' else sys.stdout) as outfile:
            if header:
                outfile.write(u'!\n! This file was automatically generated from JSBACH source file:\n')
                outfile.write(u'!    ' + re.sub(r'.*/src', 'src', os.path.realpath(os.path.abspath(inputfile))) + u'\n')
                outfile.write(u'! with md5 hash: ' + hasher.hexdigest() + u'\n!\n')
                outfile.write(u'#line 1 "{}"\n'.format(os.path.realpath(os.path.abspath(inputfile))))
            if verbose:
                if cwd:
                    print('   Pre-processing: ' + os.path.relpath(os.path.realpath(os.path.abspath(inputfile)), cwd) + ' => ' + os.path.relpath(outputfile, cwd))
                else:
                    print('   Pre-processing: ' + os.path.realpath(inputfile) + ' => ' + outputfile)
            for line in infile:
                if outputfile == '' and inputfile == 'test.f90':
                    if (sys.version_info > (3, 0)):
                        outfile.write(u'! ' + line)
                    else:
                        outfile.write(u'! ' + line.encode('utf-8'))
                matches = re.finditer(re_function, line)
                if matches:
                    for match in matches:         # match.groupt(1) contains the macro name without "dsl4jsb_"
                        try:
                            key = match.group(1).lower()
                            pat, replace = macros[key]
                            line = pat.sub(replace, line)
                            line = re.sub(r'([A-Z2]*?)_', tolower, line)
                        except:
                            print('ERROR in ' + inputfile + ': dsl4jsb macro not found: ' + match.group(1).lower())
                            sys.exit()
                if outputfile == '':
                    if (sys.version_info > (3, 0)):
                        outfile.write(line)
                    else:
                        outfile.write(line.encode('utf-8'))
                else:
                    outfile.write(line)
    if outputfile != '':
        os.utime(outputfile, (os.path.getatime(outputfile), os.path.getmtime(inputfile)))

def main(argv):
    inputfile = ''
    outputfile = ''
    inputdir = ''
    outputdir = ''
    verbose = False
    verbose_more = False
    keep = False
    header = True
    pp_add = None
    try:
        opts, args = getopt.getopt(argv,"hvkni:o:d:t:p:",["ifile=","ofile=","idir=","odir=","pp-add="])
    except getopt.GetoptError:
        print(help_text)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(help_text)
            sys.exit()
        elif opt in ("-v", "--verbose"):
            if not verbose:
                verbose = True
            else:
                verbose_more = True
        elif opt in ("-k", "--keep-dirs"):
            keep = True
        elif opt in ("-n", "--no-header"):
            header = False
        elif opt in ("-p", "--pp-add"):
            pp_add = arg
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-d", "--idir"):
            inputdir = os.path.abspath(arg)     # Make inputdir absolute
        elif opt in ("-t", "--odir"):
            outputdir = os.path.abspath(arg)    # Make outputdir absolute
            if not os.path.exists(outputdir):
                print('Output directory ' + outputdir + ' does not exist.')
                sys.exit()
    if inputdir != '' and outputdir == '':
        print('Input directory specified, but no output directory')
        sys.exit()
    if inputfile != '' and inputdir != '':
        print('Specify either input file or input directory but not both!')
        sys.exit()
    if inputdir != '' and outputdir != '':
        inputfile = ''
        outputfile = ''
    if inputfile != '' and outputfile == '' and outputdir != '':
        if pp_add:
            root, ext = os.path.splitext(os.path.basename(inputfile))
            outputfile = os.path.join(outputdir, root + pp_add + ext)
        else:
            outputfile = os.path.join(outputdir, os.path.basename(inputfile))
    if inputfile == '' and inputdir == '':
        inputfile = 'test.f90'

    cwd = os.getcwd()    # Save current directory
    if inputdir != '':
        os.chdir(inputdir)   # Change to input directory in order to get relative directory names
        for ifile in get_input_files('.'):     # gets input files including (relative to inputdir) directory name
            if keep:
                if pp_add:
                    root, ext = os.path.splitext(ifile)
                    ofile = os.path.join(outputdir, root + pp_add + ext)
                else:
                    ofile = os.path.join(outputdir, ifile)
                dir = os.path.dirname(ofile)
                if not os.path.exists(dir):
                    os.makedirs(dir)
            else:
                if pp_add:
                    root, ext = os.path.splitext(os.path.basename(ifile))
                    ofile = os.path.join(outputdir, root + pp_add + ext)
                else:
                    ofile = os.path.join(outputdir, os.path.basename(ifile))
            process_file(os.path.join(inputdir, ifile), ofile, header, verbose, verbose_more, cwd)
        os.chdir(cwd)
    else:
        process_file(inputfile, outputfile, header, verbose, verbose_more, cwd)

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

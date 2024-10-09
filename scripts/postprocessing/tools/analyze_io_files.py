#!/usr/bin/python

# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
#
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

# ---------------------------------------------------------------------------
#
# Parses asynchronous ICON output and generates a LaTeX-format diagram of I/O
# activity.
#
# VERSION 2 - Additional output of file system information: For each process
#             the time stamps of the output files are plotted.
#             For this, the stdout output file is required (mapping between
#             I/O processes and file names) and another input file containing
#             the result of "ls -l".
# 
# IMPORTANT NOTE: We implicitly assume that all file time stamps are generated
#                 in the format (hh:mm) on the _same date_!
# NOTE #2:        We synchronize time stamps from "ls -l" and from the program
#                 output by the time stamp of the first file.
#
# Example usage
#  python ~/svn/test_programs/analyze_io/analyze_io_files.py case1_7333172.log.o case1_7333172.log.e ls.txt > test.tex; pdflatex test.tex; gv test.pdf
#
# or
#
# ---------------------------------------------------------------------------

import re, sys
from sets import Set

if (len(sys.argv) < 2):
    print "ANALYZE_IO_FILES.PY"
    print "Parses asynchronous ICON output and generates a LaTeX-format diagram of I/O activity."
    print "\n Example usage"
    print "  ./analyze_io_files.py <file1>.o > test.tex; pdflatex test.tex; gv test.pdf"
    print " or"
    print "  ./analyze_io_files.py <file1>.o <file1>.e ls.txt > test.tex; pdflatex test.tex; gv test.pdf\n"
    sys.exit(0)

# Get file names from command line arguments
#   Note: You may provide identical file names for
#         "file_intvls" and "file_filenm"
file_intvls = sys.argv[1] # file containing output intervals  (stderr)

if (len(sys.argv) > 2):
    file_filenm = sys.argv[2] # file containing output file names (stdout)
    file_ls     = sys.argv[3] # file containing "ls -l" output

# utility function
def compare_tstamp(x, y):
    idx = 0
    if (x[0]==y[0]): idx = 1
    if (x[idx]<y[idx]): return -1
    else:               return +1

# ---------------------------------------------------------------------------
# PART 1: Parsing output file: Output intervals
# ---------------------------------------------------------------------------

# prepare regular expressions:

# Starting line:
# for example: "2035:#################### I/O PE 2035 starting I/O at 174646.442"
#                                                 ^ group(1)            ^ group(2)
re_start=re.compile(r'I/O PE (\d*) starting I/O at ([\d.]*)')

# Ending line
# for example: "2046:#################### I/O PE 2046 done at 174650.158"
#                                                 ^ group(1)      ^ group(2)
re_done=re.compile(r'I/O PE (\d*) done at ([\d.]*)')

# create a set of I/O PEs involved
pe_list = Set([])

# read file, put PE number and time into lists
tstart = []
tdone  = []
for line in open(file_intvls):
    m = re_start.search(line)
    if  m:
        sec=   3600.*float(m.group(2)[0:2]) \
             +   60.*float(m.group(2)[2:4]) \
             +       float(m.group(2)[4:6])
        tstart.append( (int(m.group(1)), sec) )
        pe_list.add(m.group(1))
    m = re_done.search(line)
    if  m:
        sec=   3600.*float(m.group(2)[0:2]) \
             +   60.*float(m.group(2)[2:4]) \
             +       float(m.group(2)[4:6])
        tdone.append( (int(m.group(1)), sec) )

    
# get minimum,maximum time stamp
t0 = float(sorted(tstart, key=lambda tstamp: float(    tstamp[1]))[0][1])
t1 = float(sorted(tdone,  key=lambda tstamp: float(-1.*tstamp[1]))[0][1])

# sort lists by PE number, then by their time stamps:
tstart.sort(cmp=compare_tstamp)
tdone. sort(cmp=compare_tstamp)

# assign order to list of the I/O PEs
pe_dict = dict()
index = 0
for item in sorted(pe_list):
    pe_dict[int(item)] = index
    index=index+1

# define scale for picture
lastidx = len(tstart)-1
scale   = 1.8*index / (tstart[lastidx][1]-t0)
scaley  = 0.5


# ---------------------------------------------------------------------------
# PART 2: Parsing output file: get time stamps from "ls -l" output
# ---------------------------------------------------------------------------

if (len(sys.argv) > 2):
    # Line from "ls -l" result, for example
    # "-rw-r-----    1 dfi0     de       1051319792 Mar 01 17:35 case1_3_DOM03_PL_0003.nc"
    #                                                       ^ group(1),group(2)  ^ group(3)
    re_tstamps=re.compile(r'.* (\d\d):(\d\d) (.*.nc)')
    tstamps = re_tstamps.findall(open(file_ls).read())
    tstamps_dict  = dict()
    # convert time stamp into seconds:
    for item in tstamps:
        tstamps_dict[item[2]] =   3600.*float(item[0]) \
                                +   60.*float(item[1])
    #determine the minimum time stamp
    # file_t0=min([ (tstamps_dict[x],x) for x in tstamps_dict])[0]


    # Output message, for example:
    # "mo_name_list_output: Output to case1_1_DOM06_ML_0001.nc at simulation time  0.00000000     by PE   255"
    #                                           ^ group(1)                            ^ group(2)          ^ group(3)
    re_output=re.compile(r'mo_name_list_output: Output to (.*) at simulation time  (.*)     by PE   (\d*)')
    outmsg = re_output.findall(open(file_filenm).read())
        
    # built a list of tuples: (I/O process ; time stamp ; simulation time )
    tstamps_list = []
    for item in outmsg:
        if (not tstamps_dict.has_key(item[0])):
            print "Time stamp for \"", item[0], "\" not found!"
            sys.exit(0)
        tstamps_list.append( (int(item[2]),   \
                              float(tstamps_dict[item[0]]), \
                              float(item[1])) )
            

# ---------------------------------------------------------------------------
# PART 3: Loop over list and print intervals as LaTeX code
# ---------------------------------------------------------------------------

# generate a variable for consistency check: time of "start" message for a
# given I/O PE must not precede the "done" statement of last output stage.
time_dict = dict()
for item in sorted(pe_list):
    time_dict[int(item)] = 0.

#   LaTeX header
print("\\documentclass{article}")
print("\\usepackage[a4paper,landscape]{geometry}")
print("\\usepackage{tikz}")
print("\\usetikzlibrary{calc,arrows}")
print("\\begin{document}")
print "\\large\\begin{verbatim}IO status messages, \"",file_intvls, "\"\\end{verbatim}\\vspace*{-0.5em}"
#   LaTeX draw intervals
print("\\begin{tikzpicture}[scale=0.7]")
print("  \\begin{scope}[line width=0.5pt, >=latex]")
for index,item in enumerate(tstart):
    if (item[0] != tdone[index][0]):
        print "% Error: ", index
        sys.exit(0)
    else:
        print "    \\fill[black!20] (%5.2f" % (scale*(item[1]-t0)), ",", scaley*pe_dict[item[0]], \
              ") rectangle ($(%5.2f " % (scale*(tdone[index][1]-t0)), ",", \
              scaley*pe_dict[item[0]], ")+(0,-4pt)$);", \
              "  %  [", item[1], " , ", tdone[index][1], "]"                
        print "    \\draw[-] (%5.2f" % (scale*(item[1]-t0)), ",", scaley*pe_dict[item[0]], \
              ") -- ++(0,-4pt) -- ($(%5.2f " % (scale*(tdone[index][1]-t0)), ",", \
              scaley*pe_dict[item[0]], ")+(0,-4pt)$) -- ++(0,4pt);"
        print "    \\draw[-, red] (%5.2f" % (scale*(item[1]-t0)), ",", scaley*pe_dict[item[0]], \
              ") -- ++(0,-4pt);"
    # consistency check:
    if (item[1] < time_dict[int(item[0])]):
        print "% Error: ", index, item[1], time_dict[int(item[0])]
        sys.exit(0)
    time_dict[int(item[0])] = float(tdone[index][1])
print("  \\end{scope}")
#   LaTeX draw annotation
for item in sorted(pe_list):
    print "    \\node[font=\sffamily,blue,anchor=east,yshift=-2pt] at (-0.2,", scaley*pe_dict[int(item)], ") {\\tiny IO PE \#", item, "};"
#   LaTeX draw axis and vertical grid lines
print "    \\draw[blue,-latex,yshift=-1cm] (0,0) -- (%5.2f" % (scale*(t1-t0)), ", 0);"
print "    \\draw [blue,yshift=-1cm](0,3pt) -- (0,-3pt) node[rotate=60,font=\sffamily,anchor=east] {time(sec)=", 0, "};";
print "    \\draw [blue,yshift=-1cm](%5.2f" % (scale*(t1-t0)), ",3pt) -- ++(0,-6pt) node[rotate=60,font=\sffamily,anchor=east] {", (t1-t0), "};";
print "    \\foreach \p in {0,0.1,...,1} {"
print "      \\draw[line width=0.05,black,dashed] (%5.2f" % (scale*(t1-t0)), "*\\p,0) -- ++(0, %5.2f" % (scaley*len(pe_list)), ");"
print "    };"
for p in range(1,10):
    timeval = 0.1*(t1-t0)*p # + t0
    print "      \\draw [blue,yshift=-1cm](%5.2f" % (0.1*scale*(t1-t0)*p), \
          ",3pt) -- ++(0,-6pt) node[rotate=60,font=\sffamily,anchor=east] {%5.0f" % timeval,"};";
#   LaTeX print markers for file time stamps
if (len(sys.argv) > 2):
    for item in tstamps_list:
        print "    \\draw[yshift=-4pt,latex'-,thin,line width=0.5pt,blue] (%5.2f" % (scale*(item[1]-t0)), ",", scaley*pe_dict[item[0]], \
              ") -- ++(0,-3pt);"
#   LaTeX footer
print("\\end{tikzpicture}")
print("\\end{document}")




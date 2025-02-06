/**
 * @file parmetis_wrap.c
 * @brief ParMeTis wrappers for Fortran/C interfacing
 *
 * @copyright Copyright  (C)  2010  Thomas Jahns <jahns@dkrz.de>
 *
 * @version: 1.0
 * @author: Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords: ParMeTis wrapper
 * Maintainer: Thomas Jahns <jahns@dkrz.de>
 * URL: https://www.dkrz.de/redmine/projects/scales-ppm
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are  permitted provided that the following conditions are
 * met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the DKRZ GmbH nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Commentary:
 *
 *
 *
 * Code:
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <parmetis.h>
#include <cfortran.h>

#if METIS_VER_MAJOR >= 5
typedef idx_t idxtype;
#endif

/*
 * Usually ParMeTis brings its own Fortran to C wrapper, but it does
 * unfortunately not account for the MPI Communicator conversion that
 * is needed. Therefore this file provides an interposition of that
 * function.
 */
void ParMETIS_V3_PartKway_Wrapper(
  idxtype vtxdist[],
  idxtype xadj[],
  idxtype adjncy[],
  idxtype vwgt[],
  idxtype adjwgt[],
  int *wgtflag,
  int *numflag,
  int *ncon,
  int *nparts,
  float *tpwgts,
  float *ubvec,
  int options[],
  int *edgecut,
  idxtype part[],
  MPI_Fint *comm_f)
{
  MPI_Comm comm_c;

  comm_c = MPI_Comm_f2c((MPI_Fint)*comm_f);
  ParMETIS_V3_PartKway(vtxdist, xadj, adjncy, vwgt,
                       adjwgt, wgtflag, numflag, ncon,
                       nparts, tpwgts, ubvec, options,
                       edgecut, part, &comm_c);
}

FCALLSCSUB15(ParMETIS_V3_PartKway_Wrapper, PARMETIS_V3_PARTKWAY,
             parmetis_v3_partkway,
             INTV, INTV, INTV, INTV, INTV, PINT, PINT, PINT,
             PINT, FLOATVV, FLOATV, INTV, PINT, INTV, PVOID)

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-default: "bsd"
 * End:
 */

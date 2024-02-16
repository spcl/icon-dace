/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Rene Redler <rene.redler@mpimet.mpg.de>
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yac/
 *
 * This file is part of YAC.
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
 */
#include <math.h>
#include "generate_reg2d.h"

void generate_reg2d_decomp(int num_points[2], int total_num_procs, int * num_procs) {

  int flag = num_points[0] > num_points[1];
  float ratio = (float)(num_points[flag^1]) / (float)(num_points[flag]);

  int prev_num_procs[2];

  prev_num_procs[!flag] = total_num_procs;
  prev_num_procs[flag] = 1;

  float prev_proc_ratio = (float)total_num_procs;

  for (int i = total_num_procs-1; i > 0; --i) {

    if (total_num_procs%i == 0) {

      float curr_proc_ratio = (float)(i) / (float)(total_num_procs/i);

      double temp_a = fabs(ratio - curr_proc_ratio);
      double temp_b = fabs(ratio - prev_proc_ratio);
      if (temp_a < temp_b) {
        prev_num_procs[!flag] = i;
        prev_num_procs[flag] = total_num_procs/i;
        prev_proc_ratio = curr_proc_ratio;
      }
    }
  }

  num_procs[0] = prev_num_procs[0];
  num_procs[1] = prev_num_procs[1];
}


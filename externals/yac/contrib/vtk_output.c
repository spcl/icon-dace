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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "vtk_output.h"
#include "utils.h"

#define WRITE_ASCII

enum scalars_dt {
   INT = 0,
   UINT = 1,
   FLOAT = 2,
   DOUBLE = 3,
};

const char * scalars_dt_names[4] = {
   "int",
   "unsigned_int",
   "float",
   "double"
};

struct scalars_data {
   enum scalars_dt dt;
   void * data;
   char * name;
};

struct VTK_FILE_ {

   FILE * file;

   unsigned num_cells, num_points;

   struct scalars_data * scalars_cell_data, * scalars_point_data;
   unsigned num_scalars_cell_data, num_scalars_point_data;
};

VTK_FILE * vtk_open(const char * filename, const char * title) {

   VTK_FILE * vtk_file = xmalloc(1 * sizeof(*vtk_file));
   FILE * file;

   file = xfopen(filename, "w+");
   vtk_file->file = file;
   vtk_file->scalars_cell_data = NULL;
   vtk_file->num_scalars_cell_data = 0;
   vtk_file->scalars_point_data = NULL;
   vtk_file->num_scalars_point_data = 0;

   fprintf(file, "# vtk DataFile Version 2.0\n");
   fprintf(file, "%s\n", title);
#ifdef WRITE_ASCII
   fprintf(file, "ASCII\n");
#else
   fprintf(file, "BINARY\n");
#endif
   fprintf(file, "DATASET UNSTRUCTURED_GRID\n");

   return vtk_file;
}

void vtk_write_point_data(VTK_FILE * vtk_file, double * point_data, unsigned num_points) {

   FILE * file = vtk_file->file;

   vtk_file->num_points = num_points;

   unsigned i;

   fprintf(file, "POINTS %u double\n", num_points);
#ifdef WRITE_ASCII
   for (i = 0; i < num_points; ++i)
      fprintf(file, "%f %f %f\n", point_data[i*3+0], point_data[i*3+1], point_data[i*3+2]);
#else
   fwrite(point_data, sizeof(*point_data), 3 * num_points, file);
   fprintf(file, "\n");
#endif
}

void vtk_write_cell_data(VTK_FILE * vtk_file, unsigned * cell_corners,
                         unsigned * num_points_per_cell, unsigned num_cells) {

   FILE * file = vtk_file->file;

   vtk_file->num_cells = num_cells;

   unsigned i, j, k;
   unsigned total_num_cell_corners = 0;

   for (i = 0; i < num_cells; ++i)
      total_num_cell_corners += num_points_per_cell[i];

   fprintf(file, "CELLS %u %u\n", num_cells, total_num_cell_corners + num_cells);

   for (i = 0, k = 0; i < num_cells; ++i) {

#ifdef WRITE_ASCII
      fprintf(file, "%u", num_points_per_cell[i]);
      for (j = 0; j < num_points_per_cell[i]; ++j)
         fprintf(file, " %u", cell_corners[k++]);
      fprintf(file, "\n");
#else
      fwrite(num_points_per_cell+i, sizeof(*num_points_per_cell), 1, file);
      fwrite(cell_corners+k, sizeof(*cell_corners), num_points_per_cell[i], file);
      k += num_points_per_cell[i];
#endif
   }
#ifndef WRITE_ASCII
   fprintf(file, "\n");
#endif

   fprintf(file, "CELL_TYPES %u\n", num_cells);
   for (i = 0; i < num_cells; ++i) {

    unsigned type;

    switch(num_points_per_cell[i]) {
      case(1):
        type = 1;
        break;
      case(2):
        type = 3;
        break;
      default:
        type = 7;
    };

#ifdef WRITE_ASCII
      fprintf(file, "%u\n", type);
#else
      fwrite(&type, sizeof(type), 1, file);
#endif
   }
#ifndef WRITE_ASCII
   fprintf(file, "\n");
#endif

}

static void vtk_write_cell_scalars(
  VTK_FILE * vtk_file, void * scalars, unsigned num_cells, char * name,
  enum scalars_dt dt) {

  YAC_ASSERT(
    vtk_file->num_cells == num_cells,
    "ERROR: vtk_write_cell_scalars number of cells does not match")

   int scalar_index = vtk_file->num_scalars_cell_data++;

   vtk_file->scalars_cell_data = xrealloc(vtk_file->scalars_cell_data,
                                          vtk_file->num_scalars_cell_data *
                                          sizeof(*(vtk_file->scalars_cell_data)));

   size_t data_size = num_cells;

   switch(dt) {
      case(INT) :
         data_size *= sizeof(int);
         break;
      case(UINT) :
         data_size *= sizeof(unsigned);
         break;
      case(FLOAT) :
         data_size *= sizeof(float);
         break;
      case(DOUBLE) :
         data_size *= sizeof(double);
         break;
   }

   vtk_file->scalars_cell_data[scalar_index].dt = dt;
   vtk_file->scalars_cell_data[scalar_index].data = xmalloc(data_size);
   memcpy(vtk_file->scalars_cell_data[scalar_index].data, scalars, data_size);
   vtk_file->scalars_cell_data[scalar_index].name = xmalloc((strlen(name) + 1));
   memcpy(vtk_file->scalars_cell_data[scalar_index].name, name, strlen(name) + 1);
}

static void vtk_write_point_scalars(
  VTK_FILE * vtk_file, void * scalars, unsigned num_points, char * name,
  enum scalars_dt dt) {

  YAC_ASSERT(
    vtk_file->num_points == num_points,
    "ERROR: vtk_write_point_scalars number of points does not match")

   int scalar_index = vtk_file->num_scalars_point_data++;

   vtk_file->scalars_point_data = xrealloc(vtk_file->scalars_point_data,
                                           vtk_file->num_scalars_point_data *
                                           sizeof(*(vtk_file->scalars_point_data)));

   size_t data_size = num_points;

   switch(dt) {
      case(INT) :
         data_size *= sizeof(int);
         break;
      case(UINT) :
         data_size *= sizeof(unsigned);
         break;
      case(FLOAT) :
         data_size *= sizeof(float);
         break;
      case(DOUBLE) :
         data_size *= sizeof(double);
         break;
   }

   vtk_file->scalars_point_data[scalar_index].dt = dt;
   vtk_file->scalars_point_data[scalar_index].data = xmalloc(data_size);
   memcpy(vtk_file->scalars_point_data[scalar_index].data, scalars, data_size);
   vtk_file->scalars_point_data[scalar_index].name = xmalloc((strlen(name) + 1));
   memcpy(vtk_file->scalars_point_data[scalar_index].name, name, strlen(name) + 1);
}

void vtk_write_cell_scalars_uint(VTK_FILE * vtk_file, unsigned * scalars,
                                 unsigned num_cells, char * name) {

   vtk_write_cell_scalars(vtk_file, scalars, num_cells, name, UINT);
}

void vtk_write_cell_scalars_int(VTK_FILE * vtk_file, int * scalars,
                                unsigned num_cells, char * name) {

   vtk_write_cell_scalars(vtk_file, scalars, num_cells, name, INT);
}

void vtk_write_cell_scalars_float(VTK_FILE * vtk_file, float * scalars,
                                  unsigned num_cells, char * name) {

   vtk_write_cell_scalars(vtk_file, scalars, num_cells, name, FLOAT);
}

void vtk_write_cell_scalars_double(VTK_FILE * vtk_file, double * scalars,
                                   unsigned num_cells, char * name) {

   vtk_write_cell_scalars(vtk_file, scalars, num_cells, name, DOUBLE);
}

void vtk_write_point_scalars_uint(VTK_FILE * vtk_file, unsigned * scalars,
                                  unsigned num_points, char * name) {

   vtk_write_point_scalars(vtk_file, scalars, num_points, name, UINT);
}

void vtk_write_point_scalars_int(VTK_FILE * vtk_file, int * scalars,
                                 unsigned num_points, char * name) {

   vtk_write_point_scalars(vtk_file, scalars, num_points, name, INT);
}

void vtk_write_point_scalars_float(VTK_FILE * vtk_file, float * scalars,
                                   unsigned num_points, char * name) {

   vtk_write_point_scalars(vtk_file, scalars, num_points, name, FLOAT);
}

void vtk_write_point_scalars_double(VTK_FILE * vtk_file, double * scalars,
                                    unsigned num_points, char * name) {

   vtk_write_point_scalars(vtk_file, scalars, num_points, name, DOUBLE);
}

static void vtk_flush_scalar_data(FILE * file, struct scalars_data * data,
                                  unsigned num_scalars_data,
                                  unsigned num_data_points) {


   for (unsigned i = 0; i < num_scalars_data; ++i) {

      fprintf(file, "SCALARS %s %s 1\n", data[i].name,
              scalars_dt_names[data[i].dt]);
      fprintf(file, "LOOKUP_TABLE default\n");

      for (unsigned j = 0; j < num_data_points; ++j) {
#ifdef WRITE_ASCII
         switch (data[i].dt) {

            case (INT) :
               fprintf(file, "%d\n", ((int*)(data[i].data))[j]);
               break;
            case (UINT) :
               fprintf(file, "%u\n", ((unsigned*)(data[i].data))[j]);
               break;
            case (FLOAT) :
               fprintf(file, "%.8f\n", ((float*)(data[i].data))[j]);
               //fprintf(file, "%f\n", ((float*)(data[i].data))[j]);
               break;
            case (DOUBLE):
               fprintf(file, "%.18lf\n", ((double*)(data[i].data))[j]);
               //fprintf(file, "%lf\n", ((double*)(data[i].data))[j]);
               break;
         }
#else
         switch (data[i].dt) {

            case (INT) :
               fwrite(data[i].data, sizeof(int), 1, file);
               break;
            case (UINT) :
               fwrite(data[i].data, sizeof(unsigned), 1, file);
               break;
            case (FLOAT) :
               fwrite(data[i].data, sizeof(float), 1, file);
               break;
            case (DOUBLE):
               fwrite(data[i].data, sizeof(double), 1, file);
               break;
         }
         fprintf(file, "\n");
#endif
      }
      free(data[i].data);
      free(data[i].name);
   }
}

static void vtk_flush_scalar_cell_data(VTK_FILE * vtk_file) {

   FILE * file = vtk_file->file;

   fprintf(file, "CELL_DATA %u\n", vtk_file->num_cells);

   fprintf(file, "SCALARS cell_ids int 1\n");
   fprintf(file, "LOOKUP_TABLE default\n");
   for (unsigned i = 0; i < vtk_file->num_cells; ++i) {
#ifdef WRITE_ASCII
      fprintf(file, "%u\n", i);
#else
      fwrite(&i, sizeof(i), 1, file);
#endif
   }
#ifndef WRITE_ASCII
   fprintf(file, "\n");
#endif

   if (vtk_file->num_scalars_cell_data > 0)
      vtk_flush_scalar_data(file, vtk_file->scalars_cell_data,
                            vtk_file->num_scalars_cell_data, vtk_file->num_cells);
}

static void vtk_flush_scalar_point_data(VTK_FILE * vtk_file) {

   FILE * file = vtk_file->file;

   fprintf(file, "POINT_DATA %u\n", vtk_file->num_points);

   fprintf(file, "SCALARS point_ids int 1\n");
   fprintf(file, "LOOKUP_TABLE default\n");
   for (unsigned i = 0; i < vtk_file->num_points; ++i) {
#ifdef WRITE_ASCII
      fprintf(file, "%u\n", i);
#else
      fwrite(&i, sizeof(i), 1, file);
#endif
   }
#ifndef WRITE_ASCII
   fprintf(file, "\n");
#endif

   if (vtk_file->num_scalars_point_data > 0)
      vtk_flush_scalar_data(file, vtk_file->scalars_point_data,
                            vtk_file->num_scalars_point_data, vtk_file->num_points);
}

void vtk_close(VTK_FILE * vtk_file) {

   FILE * file = vtk_file->file;

   vtk_flush_scalar_cell_data(vtk_file);
   free(vtk_file->scalars_cell_data);
   vtk_flush_scalar_point_data(vtk_file);
   free(vtk_file->scalars_point_data);

   xfclose(file);

   free(vtk_file);
}



#line 1 "util_string_parse.rl"
// ICON
//
// ---------------------------------------------------------------
// Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
// Contact information: icon-model.org
//
// See AUTHORS.TXT for a list of authors
// See LICENSES/ for license information
// SPDX-License-Identifier: BSD-3-Clause
// ---------------------------------------------------------------

/* --------------------------------------------------------------------- *
 * String parser based on a Finite State Machine (FSM).                  *
 * This file must be processed with Ragel to produce the final C code.   *
 *                                                                       *
 * PARSER "parse_intlist"                                                *
 *    Parse a string with an integer number list with ranges.            *
 *                                                                       *
 *    Allowed is a comma- (or semicolon-) separated list of integers,    *
 *    and of integer ranges like "10...20".                              *
 *    One may also use the keyword "nlev" to denote the maximum integer  *
 *    (or, equivalently, "n" or "N").                                    *
 *    Furthermore, arithmetic expressions like "(nlev - 2)" are          *
 *    possible.                                                          *
 *                                                                       *
 *    Basic example:                                                     *
 *       parse_line = "1,3,5...10,20...nlev"                             *
 *    More complex example:                                              *
 *       parse_line = "1,2, 10 ...22;2;16-(3+11), nlev-2,16-(2+10);5"    *
 *                                                                       *
 * Ragel State Machine Compiler version 7.0.4 was used                   *
 * --------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_BUF_LEN  32768
#define MAX_STACK       10
#define SYM_PLUS      9998
#define SYM_MINUS     9999

/* All data needed to define the state of the Finite State Machine parser. */
struct t_parsedata {
  int  cs, act, buflen;           /* FSM state, buffer length     */
  char *ts, *te;                  /* token start, end             */
  char buffer[MAX_BUF_LEN + 1];   /* token string buffer          */
  int  value;                     /* integer value storage        */
  int  maxval;                    /* maximum value for integers   */
  int  nlev_val;                  /* value to replace for "nlev"  */
  int *out_values;                /* pointer to result list       */
  int  vidx, oidx;                /* value and operand stack size */
  int  vstack[MAX_STACK];         /* value stack                  */
  int  ostack[MAX_STACK];         /* operator stack               */
};


/* --------------------------------------------------------------------- *
 * DEFINITION OF RAGEL FINITE STATE MACHINE "parse_intlist"
 * --------------------------------------------------------------------- */


#line 159 "util_string_parse.rl"



#line 70 "util_string_parse.c"
static const char _parse_intlist_actions[] = {
	0, 1, 0, 1, 1, 1, 2, 1, 
	3, 1, 4, 1, 5, 1, 6, 1, 
	7, 1, 8, 1, 9, 1, 10, 1, 
	11, 1, 12, 1, 13, 2, 1, 0, 
	2, 2, 3, 2, 2, 5, 2, 2, 
	9, 2, 2, 10, 2, 3, 5, 2, 
	3, 9, 2, 3, 10, 2, 4, 3, 
	2, 4, 5, 2, 4, 9, 2, 4, 
	10, 2, 5, 3, 2, 8, 0, 2, 
	8, 1, 2, 8, 2, 2, 8, 3, 
	2, 8, 4, 3, 2, 3, 5, 3, 
	2, 3, 9, 3, 2, 3, 10, 3, 
	2, 5, 3, 3, 4, 3, 5, 3, 
	4, 3, 9, 3, 4, 3, 10, 3, 
	4, 5, 3, 3, 8, 1, 0, 3, 
	8, 2, 3, 3, 8, 4, 3
};

static const short _parse_intlist_key_offsets[] = {
	0, 0, 9, 17, 24, 32, 38, 44, 
	52, 53, 54, 62, 69, 77, 83, 88, 
	96, 103, 111, 117, 124, 132, 138, 144, 
	151, 152, 153, 159, 166, 167, 168, 176, 
	177, 178, 179, 180, 187, 195, 201, 207, 
	214, 215, 216, 222, 229, 230, 231, 232, 
	233, 234, 235, 242, 250, 256, 262, 269, 
	270, 271, 277, 284, 285, 286, 287, 288, 
	289, 290, 297, 305, 311, 318, 326, 332, 
	338, 345, 346, 347, 353, 360, 361, 362, 
	372, 375, 383, 389, 395, 403, 409, 414, 
	419, 426, 431, 438, 443, 448, 454, 459, 
	465, 472, 477, 482, 488, 494, 501, 509, 
	515, 521, 528, 534, 541
};

static const char _parse_intlist_trans_keys[] = {
	32, 44, 59, 78, 110, 9, 13, 49, 
	57, 32, 40, 78, 110, 9, 13, 49, 
	57, 32, 78, 110, 9, 13, 49, 57, 
	32, 41, 43, 45, 9, 13, 48, 57, 
	32, 41, 43, 45, 9, 13, 32, 43, 
	45, 46, 9, 13, 32, 40, 78, 110, 
	9, 13, 49, 57, 46, 46, 32, 40, 
	78, 110, 9, 13, 49, 57, 32, 78, 
	110, 9, 13, 49, 57, 32, 41, 43, 
	45, 9, 13, 48, 57, 32, 41, 43, 
	45, 9, 13, 32, 43, 45, 9, 13, 
	32, 40, 78, 110, 9, 13, 49, 57, 
	32, 78, 110, 9, 13, 49, 57, 32, 
	41, 43, 45, 9, 13, 48, 57, 32, 
	41, 43, 45, 9, 13, 32, 78, 110, 
	9, 13, 49, 57, 32, 41, 43, 45, 
	9, 13, 48, 57, 32, 41, 43, 45, 
	9, 13, 32, 41, 43, 45, 9, 13, 
	32, 41, 43, 45, 108, 9, 13, 101, 
	118, 32, 41, 43, 45, 9, 13, 32, 
	41, 43, 45, 108, 9, 13, 101, 118, 
	32, 40, 78, 110, 9, 13, 49, 57, 
	101, 118, 101, 118, 32, 78, 110, 9, 
	13, 49, 57, 32, 41, 43, 45, 9, 
	13, 48, 57, 32, 41, 43, 45, 9, 
	13, 32, 41, 43, 45, 9, 13, 32, 
	41, 43, 45, 108, 9, 13, 101, 118, 
	32, 41, 43, 45, 9, 13, 32, 41, 
	43, 45, 108, 9, 13, 101, 118, 101, 
	118, 101, 118, 32, 78, 110, 9, 13, 
	49, 57, 32, 41, 43, 45, 9, 13, 
	48, 57, 32, 41, 43, 45, 9, 13, 
	32, 41, 43, 45, 9, 13, 32, 41, 
	43, 45, 108, 9, 13, 101, 118, 32, 
	41, 43, 45, 9, 13, 32, 41, 43, 
	45, 108, 9, 13, 101, 118, 101, 118, 
	101, 118, 32, 78, 110, 9, 13, 49, 
	57, 32, 41, 43, 45, 9, 13, 48, 
	57, 32, 41, 43, 45, 9, 13, 32, 
	78, 110, 9, 13, 49, 57, 32, 41, 
	43, 45, 9, 13, 48, 57, 32, 41, 
	43, 45, 9, 13, 32, 41, 43, 45, 
	9, 13, 32, 41, 43, 45, 108, 9, 
	13, 101, 118, 32, 41, 43, 45, 9, 
	13, 32, 41, 43, 45, 108, 9, 13, 
	101, 118, 32, 40, 44, 59, 78, 110, 
	9, 13, 49, 57, 32, 9, 13, 32, 
	43, 45, 46, 9, 13, 48, 57, 32, 
	43, 45, 46, 9, 13, 32, 43, 45, 
	46, 9, 13, 32, 43, 45, 46, 9, 
	13, 48, 57, 32, 43, 45, 46, 9, 
	13, 32, 43, 45, 9, 13, 32, 43, 
	45, 9, 13, 32, 43, 45, 9, 13, 
	48, 57, 32, 43, 45, 9, 13, 32, 
	43, 45, 9, 13, 48, 57, 32, 43, 
	45, 9, 13, 32, 43, 45, 9, 13, 
	32, 43, 45, 108, 9, 13, 32, 43, 
	45, 9, 13, 32, 43, 45, 108, 9, 
	13, 32, 43, 45, 9, 13, 48, 57, 
	32, 43, 45, 9, 13, 32, 43, 45, 
	9, 13, 32, 43, 45, 108, 9, 13, 
	32, 43, 45, 46, 9, 13, 32, 43, 
	45, 46, 108, 9, 13, 32, 43, 45, 
	46, 9, 13, 48, 57, 32, 43, 45, 
	46, 9, 13, 32, 43, 45, 46, 9, 
	13, 32, 43, 45, 46, 108, 9, 13, 
	32, 43, 45, 46, 9, 13, 32, 43, 
	45, 46, 108, 9, 13, 32, 43, 45, 
	46, 9, 13, 0
};

static const char _parse_intlist_single_lengths[] = {
	0, 5, 4, 3, 4, 4, 4, 4, 
	1, 1, 4, 3, 4, 4, 3, 4, 
	3, 4, 4, 3, 4, 4, 4, 5, 
	1, 1, 4, 5, 1, 1, 4, 1, 
	1, 1, 1, 3, 4, 4, 4, 5, 
	1, 1, 4, 5, 1, 1, 1, 1, 
	1, 1, 3, 4, 4, 4, 5, 1, 
	1, 4, 5, 1, 1, 1, 1, 1, 
	1, 3, 4, 4, 3, 4, 4, 4, 
	5, 1, 1, 4, 5, 1, 1, 6, 
	1, 4, 4, 4, 4, 4, 3, 3, 
	3, 3, 3, 3, 3, 4, 3, 4, 
	3, 3, 3, 4, 4, 5, 4, 4, 
	4, 5, 4, 5, 4
};

static const char _parse_intlist_range_lengths[] = {
	0, 2, 2, 2, 2, 1, 1, 2, 
	0, 0, 2, 2, 2, 1, 1, 2, 
	2, 2, 1, 2, 2, 1, 1, 1, 
	0, 0, 1, 1, 0, 0, 2, 0, 
	0, 0, 0, 2, 2, 1, 1, 1, 
	0, 0, 1, 1, 0, 0, 0, 0, 
	0, 0, 2, 2, 1, 1, 1, 0, 
	0, 1, 1, 0, 0, 0, 0, 0, 
	0, 2, 2, 1, 2, 2, 1, 1, 
	1, 0, 0, 1, 1, 0, 0, 2, 
	1, 2, 1, 1, 2, 1, 1, 1, 
	2, 1, 2, 1, 1, 1, 1, 1, 
	2, 1, 1, 1, 1, 1, 2, 1, 
	1, 1, 1, 1, 1
};

static const short _parse_intlist_index_offsets[] = {
	0, 0, 8, 15, 21, 28, 34, 40, 
	47, 49, 51, 58, 64, 71, 77, 82, 
	89, 95, 102, 108, 114, 121, 127, 133, 
	140, 142, 144, 150, 157, 159, 161, 168, 
	170, 172, 174, 176, 182, 189, 195, 201, 
	208, 210, 212, 218, 225, 227, 229, 231, 
	233, 235, 237, 243, 250, 256, 262, 269, 
	271, 273, 279, 286, 288, 290, 292, 294, 
	296, 298, 304, 311, 317, 323, 330, 336, 
	342, 349, 351, 353, 359, 366, 368, 370, 
	379, 382, 389, 395, 401, 408, 414, 419, 
	424, 430, 435, 441, 446, 451, 457, 462, 
	468, 474, 479, 484, 490, 496, 503, 510, 
	516, 522, 529, 535, 542
};

static const unsigned char _parse_intlist_indicies[] = {
	0, 2, 2, 4, 5, 0, 3, 1, 
	7, 8, 10, 11, 7, 9, 6, 8, 
	13, 14, 8, 12, 6, 15, 16, 17, 
	17, 15, 18, 6, 19, 20, 21, 21, 
	19, 6, 22, 23, 23, 24, 22, 6, 
	25, 8, 27, 28, 25, 26, 6, 29, 
	6, 30, 6, 31, 32, 34, 35, 31, 
	33, 6, 36, 38, 39, 36, 37, 6, 
	40, 41, 42, 42, 40, 43, 6, 44, 
	45, 46, 46, 44, 6, 48, 49, 49, 
	48, 47, 50, 51, 53, 54, 50, 52, 
	47, 51, 56, 57, 51, 55, 47, 58, 
	59, 60, 60, 58, 61, 47, 62, 63, 
	64, 64, 62, 47, 65, 67, 68, 65, 
	66, 47, 69, 70, 71, 71, 69, 72, 
	47, 73, 74, 75, 75, 73, 47, 76, 
	77, 78, 78, 76, 47, 76, 77, 78, 
	78, 79, 76, 47, 80, 47, 67, 47, 
	81, 82, 83, 83, 81, 47, 81, 82, 
	83, 83, 84, 81, 47, 85, 47, 56, 
	47, 86, 51, 88, 89, 86, 87, 47, 
	90, 47, 88, 47, 91, 47, 53, 47, 
	92, 94, 95, 92, 93, 6, 96, 97, 
	98, 98, 96, 99, 6, 100, 101, 102, 
	102, 100, 6, 103, 104, 105, 105, 103, 
	6, 103, 104, 105, 105, 106, 103, 6, 
	107, 6, 94, 6, 108, 109, 110, 110, 
	108, 6, 108, 109, 110, 110, 111, 108, 
	6, 112, 6, 38, 6, 113, 47, 114, 
	47, 115, 6, 27, 6, 116, 118, 119, 
	116, 117, 6, 120, 121, 122, 122, 120, 
	123, 6, 124, 125, 126, 126, 124, 6, 
	127, 128, 129, 129, 127, 6, 127, 128, 
	129, 129, 130, 127, 6, 131, 6, 118, 
	6, 132, 133, 134, 134, 132, 6, 132, 
	133, 134, 134, 135, 132, 6, 136, 6, 
	13, 6, 137, 6, 10, 6, 138, 6, 
	4, 6, 139, 141, 142, 139, 140, 1, 
	143, 144, 145, 145, 143, 146, 1, 147, 
	148, 149, 149, 147, 1, 150, 152, 153, 
	150, 151, 1, 154, 155, 156, 156, 154, 
	157, 1, 158, 159, 160, 160, 158, 1, 
	161, 162, 163, 163, 161, 1, 161, 162, 
	163, 163, 164, 161, 1, 165, 1, 152, 
	1, 166, 167, 168, 168, 166, 1, 166, 
	167, 168, 168, 169, 166, 1, 170, 1, 
	141, 1, 171, 172, 2, 2, 173, 174, 
	171, 3, 1, 2, 2, 175, 177, 178, 
	178, 179, 177, 180, 176, 182, 183, 183, 
	24, 182, 181, 185, 186, 186, 187, 185, 
	184, 189, 190, 190, 191, 189, 192, 188, 
	193, 194, 194, 187, 193, 184, 48, 49, 
	49, 48, 195, 197, 198, 198, 197, 196, 
	200, 201, 201, 200, 202, 199, 203, 204, 
	204, 203, 196, 205, 206, 206, 205, 207, 
	199, 208, 209, 209, 208, 196, 211, 212, 
	212, 211, 210, 211, 212, 212, 213, 211, 
	210, 214, 215, 215, 214, 210, 214, 215, 
	215, 216, 214, 210, 218, 219, 219, 218, 
	220, 217, 221, 222, 222, 221, 195, 224, 
	225, 225, 224, 223, 224, 225, 225, 226, 
	224, 223, 228, 229, 229, 230, 228, 227, 
	228, 229, 229, 230, 231, 228, 227, 232, 
	233, 233, 191, 232, 234, 188, 235, 236, 
	236, 187, 235, 184, 237, 238, 238, 230, 
	237, 227, 237, 238, 238, 230, 239, 237, 
	227, 241, 242, 242, 243, 241, 240, 241, 
	242, 242, 243, 244, 241, 240, 22, 23, 
	23, 24, 22, 181, 0
};

static const char _parse_intlist_trans_targs[] = {
	1, 0, 80, 81, 106, 107, 79, 2, 
	3, 102, 104, 105, 4, 57, 58, 5, 
	83, 50, 4, 5, 83, 50, 6, 7, 
	8, 7, 84, 100, 101, 9, 10, 10, 
	11, 96, 98, 99, 11, 12, 42, 43, 
	13, 86, 35, 12, 13, 86, 35, 79, 
	14, 15, 15, 16, 88, 94, 95, 17, 
	26, 27, 18, 87, 19, 17, 18, 87, 
	19, 19, 20, 22, 23, 21, 87, 19, 
	20, 21, 87, 19, 22, 87, 19, 24, 
	25, 26, 87, 19, 28, 29, 30, 90, 
	92, 93, 32, 34, 35, 36, 38, 39, 
	37, 86, 35, 36, 37, 86, 35, 38, 
	86, 35, 40, 41, 42, 86, 35, 44, 
	45, 47, 98, 49, 50, 51, 53, 54, 
	52, 83, 50, 51, 52, 83, 50, 53, 
	83, 50, 55, 56, 57, 83, 50, 59, 
	60, 62, 64, 65, 66, 75, 76, 67, 
	108, 68, 66, 67, 108, 68, 68, 69, 
	71, 72, 70, 108, 68, 69, 70, 108, 
	68, 71, 108, 68, 73, 74, 75, 108, 
	68, 77, 78, 1, 65, 106, 107, 79, 
	79, 82, 2, 8, 81, 79, 82, 2, 
	79, 6, 7, 8, 79, 85, 2, 8, 
	84, 85, 2, 79, 79, 14, 15, 79, 
	89, 30, 88, 89, 30, 91, 30, 90, 
	91, 30, 79, 92, 30, 31, 94, 30, 
	33, 79, 97, 30, 96, 97, 30, 79, 
	98, 30, 46, 79, 100, 2, 8, 48, 
	103, 2, 102, 103, 2, 104, 2, 61, 
	79, 106, 2, 8, 63
};

static const char _parse_intlist_trans_actions[] = {
	0, 0, 0, 115, 17, 17, 25, 0, 
	0, 115, 17, 17, 29, 0, 0, 5, 
	74, 35, 1, 0, 17, 11, 0, 11, 
	0, 0, 115, 17, 17, 0, 0, 3, 
	3, 115, 71, 71, 0, 29, 0, 0, 
	5, 74, 35, 1, 0, 17, 11, 27, 
	0, 11, 0, 0, 115, 17, 17, 29, 
	0, 0, 5, 74, 35, 1, 0, 17, 
	11, 0, 29, 0, 0, 32, 119, 83, 
	1, 7, 77, 44, 53, 123, 99, 0, 
	0, 9, 80, 56, 0, 0, 0, 115, 
	17, 17, 0, 0, 0, 29, 0, 0, 
	32, 119, 83, 1, 7, 77, 44, 53, 
	123, 99, 0, 0, 9, 80, 56, 0, 
	0, 0, 17, 0, 0, 29, 0, 0, 
	32, 119, 83, 1, 7, 77, 44, 53, 
	123, 99, 0, 0, 9, 80, 56, 0, 
	0, 0, 0, 0, 29, 0, 0, 5, 
	74, 35, 1, 0, 17, 11, 0, 29, 
	0, 0, 32, 119, 83, 1, 7, 77, 
	44, 53, 123, 99, 0, 0, 9, 80, 
	56, 0, 0, 3, 3, 71, 71, 23, 
	38, 74, 35, 5, 68, 19, 17, 11, 
	47, 7, 44, 7, 87, 119, 95, 32, 
	68, 77, 65, 21, 50, 7, 44, 91, 
	119, 95, 68, 77, 65, 119, 83, 68, 
	77, 44, 107, 123, 99, 0, 123, 111, 
	0, 41, 74, 35, 68, 17, 11, 62, 
	80, 56, 0, 103, 123, 111, 53, 0, 
	119, 83, 68, 77, 44, 123, 99, 0, 
	59, 80, 56, 9, 0
};

static const char _parse_intlist_to_state_actions[] = {
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 13, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0
};

static const char _parse_intlist_from_state_actions[] = {
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 15, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0
};

static const short _parse_intlist_eof_trans[] = {
	0, 0, 7, 7, 7, 7, 7, 7, 
	7, 7, 7, 7, 7, 7, 48, 48, 
	48, 48, 48, 48, 48, 48, 48, 48, 
	48, 48, 48, 48, 48, 48, 48, 48, 
	48, 48, 48, 7, 7, 7, 7, 7, 
	7, 7, 7, 7, 7, 7, 48, 48, 
	7, 7, 7, 7, 7, 7, 7, 7, 
	7, 7, 7, 7, 7, 7, 7, 7, 
	7, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	176, 177, 182, 185, 189, 185, 196, 197, 
	200, 197, 200, 197, 211, 211, 211, 211, 
	218, 196, 224, 224, 228, 228, 189, 185, 
	228, 228, 241, 241, 182
};

static const int parse_intlist_start = 79;
static const int parse_intlist_first_final = 79;
static const int parse_intlist_error = 0;

static const int parse_intlist_en_main = 79;


#line 162 "util_string_parse.rl"


void do_parse_intlist(const char *in_parse_line, const int nvalues, const int nlev_val, 
                      int *out_values, int* ierr) 
{
  int i;
  char* parse_line =  strdup(in_parse_line);
  char *p  = parse_line;                                     /* input start            */
  char *pe = parse_line + strlen(parse_line), *eof = pe ;  /* pointer to input end.  */

  struct t_parsedata parsedata = {
    .maxval   = nvalues, .out_values = out_values,
    .nlev_val = nlev_val,
    .vidx     = 0,       .oidx       = 0};                 /* input scanner state    */
  struct t_parsedata *data = &parsedata;
  for (i=0; i<nvalues; ++i)  out_values[i] = 0;
  
  // initialize Finite State Machine
  
#line 448 "util_string_parse.c"
	{
	 data->cs = parse_intlist_start;
	 data->ts = 0;
	 data->te = 0;
	 data->act = 0;
	}

#line 181 "util_string_parse.rl"

  // execute Finite State Machine
  
#line 460 "util_string_parse.c"
	{
	int _klen;
	unsigned int _trans;
	const char *_acts;
	unsigned int _nacts;
	const char *_keys;

	if ( p == pe )
		goto _test_eof;
	if (  data->cs == 0 )
		goto _out;
_resume:
	_acts = _parse_intlist_actions + _parse_intlist_from_state_actions[ data->cs];
	_nacts = (unsigned int) *_acts++;
	while ( _nacts-- > 0 ) {
		switch ( *_acts++ ) {
	case 7:
#line 1 "NONE"
	{ data->ts = p;}
	break;
#line 481 "util_string_parse.c"
		}
	}

	_keys = _parse_intlist_trans_keys + _parse_intlist_key_offsets[ data->cs];
	_trans = _parse_intlist_index_offsets[ data->cs];

	_klen = _parse_intlist_single_lengths[ data->cs];
	if ( _klen > 0 ) {
		const char *_lower = _keys;
		const char *_mid;
		const char *_upper = _keys + _klen - 1;
		while (1) {
			if ( _upper < _lower )
				break;

			_mid = _lower + ((_upper-_lower) >> 1);
			if ( (*p) < *_mid )
				_upper = _mid - 1;
			else if ( (*p) > *_mid )
				_lower = _mid + 1;
			else {
				_trans += (unsigned int)(_mid - _keys);
				goto _match;
			}
		}
		_keys += _klen;
		_trans += _klen;
	}

	_klen = _parse_intlist_range_lengths[ data->cs];
	if ( _klen > 0 ) {
		const char *_lower = _keys;
		const char *_mid;
		const char *_upper = _keys + (_klen<<1) - 2;
		while (1) {
			if ( _upper < _lower )
				break;

			_mid = _lower + (((_upper-_lower) >> 1) & ~1);
			if ( (*p) < _mid[0] )
				_upper = _mid - 2;
			else if ( (*p) > _mid[1] )
				_lower = _mid + 2;
			else {
				_trans += (unsigned int)((_mid - _keys)>>1);
				goto _match;
			}
		}
		_trans += _klen;
	}

_match:
	_trans = _parse_intlist_indicies[_trans];
_eof_trans:
	 data->cs = _parse_intlist_trans_targs[_trans];

	if ( _parse_intlist_trans_actions[_trans] == 0 )
		goto _again;

	_acts = _parse_intlist_actions + _parse_intlist_trans_actions[_trans];
	_nacts = (unsigned int) *_acts++;
	while ( _nacts-- > 0 )
	{
		switch ( *_acts++ )
		{
	case 0:
#line 67 "util_string_parse.rl"
	{ if ((data->buflen < MAX_BUF_LEN) && ((*p)!='\n')) {
	  data->buffer[data->buflen++] = (*p); }
 }
	break;
	case 1:
#line 71 "util_string_parse.rl"
	{ data->buflen = 0; }
	break;
	case 2:
#line 73 "util_string_parse.rl"
	{ 
    /* buffer string with "\0" character appended. */
    if (data->buflen < MAX_BUF_LEN) data->buffer[data->buflen++] = 0; 
    data->vstack[data->vidx++] = atol(data->buffer); 
 }
	break;
	case 3:
#line 79 "util_string_parse.rl"
	{ 
   if (data->oidx != 0) {
     int
       numB = data->vstack[--data->vidx] ,
       op   = data->ostack[--data->oidx] ,
       numA = data->vstack[--data->vidx];
     switch (op) {
     case SYM_PLUS:
       numA += numB;
       break;
     case SYM_MINUS:
       numA -= numB;
       break;
     default:
       data->cs = parse_intlist_error;
     }
     data->vstack[data->vidx++] = numA;
   }
 }
	break;
	case 4:
#line 98 "util_string_parse.rl"
	{ data->vstack[data->vidx++] = data->nlev_val; }
	break;
	case 5:
#line 100 "util_string_parse.rl"
	{ 
   switch ((*p)) {
   case '+':
     data->ostack[data->oidx++] = SYM_PLUS;
     break;
   case '-':
     data->ostack[data->oidx++] = SYM_MINUS;
     break;
   default:
     data->cs = parse_intlist_error;
   }
 }
	break;
	case 8:
#line 1 "NONE"
	{ data->te = p+1;}
	break;
	case 9:
#line 113 "util_string_parse.rl"
	{ data->te = p;p--;{ 
   int value = data->vstack[--data->vidx];
   if (value > data->maxval) {data->cs = parse_intlist_error;}
   data->out_values[value] = 1;
 }}
	break;
	case 10:
#line 119 "util_string_parse.rl"
	{ data->te = p;p--;{ 
   int i,
     range_end   = data->vstack[--data->vidx],
     range_start = data->vstack[--data->vidx];
   if ((range_end   > data->maxval)  ||
       (range_start > data->maxval)  ||
       (range_start > range_end)) {data->cs = parse_intlist_error;}
   for (i=range_start; i<=range_end; ++i)  data->out_values[i] = 1;
 }}
	break;
	case 11:
#line 156 "util_string_parse.rl"
	{ data->te = p;p--;}
	break;
	case 12:
#line 113 "util_string_parse.rl"
	{{p = (( data->te))-1;}{ 
   int value = data->vstack[--data->vidx];
   if (value > data->maxval) {data->cs = parse_intlist_error;}
   data->out_values[value] = 1;
 }}
	break;
	case 13:
#line 119 "util_string_parse.rl"
	{{p = (( data->te))-1;}{ 
   int i,
     range_end   = data->vstack[--data->vidx],
     range_start = data->vstack[--data->vidx];
   if ((range_end   > data->maxval)  ||
       (range_start > data->maxval)  ||
       (range_start > range_end)) {data->cs = parse_intlist_error;}
   for (i=range_start; i<=range_end; ++i)  data->out_values[i] = 1;
 }}
	break;
#line 654 "util_string_parse.c"
		}
	}

_again:
	_acts = _parse_intlist_actions + _parse_intlist_to_state_actions[ data->cs];
	_nacts = (unsigned int) *_acts++;
	while ( _nacts-- > 0 ) {
		switch ( *_acts++ ) {
	case 6:
#line 1 "NONE"
	{ data->ts = 0;}
	break;
#line 667 "util_string_parse.c"
		}
	}

	if (  data->cs == 0 )
		goto _out;
	if ( ++p != pe )
		goto _resume;
	_test_eof: {}
	if ( p == eof )
	{
	if ( _parse_intlist_eof_trans[ data->cs] > 0 ) {
		_trans = _parse_intlist_eof_trans[ data->cs] - 1;
		goto _eof_trans;
	}
	}

	_out: {}
	}

#line 184 "util_string_parse.rl"

  (*ierr) = ( data->cs == parse_intlist_error );
  free(parse_line);
}


/* --------------------------------------------------------------------- *
 * Main routine
 * --------------------------------------------------------------------- */
/*
void main(void)
{
  char *parse_line = "1,3,5...10,20...nlev";

  const int maxval = 50;
  int   values[maxval+1];
  int   i,ierr   = 0;
  
  do_parse_intlist(parse_line, maxval, values, &ierr);
  if ( ierr != 0 ) fprintf(stderr, "PARSE ERROR\n" );

  printf("\"%s\" is translated into \n", parse_line);
  for (i=0; i<maxval; ++i)  printf("%i ", values[i]);
  printf("\n");
}
*/

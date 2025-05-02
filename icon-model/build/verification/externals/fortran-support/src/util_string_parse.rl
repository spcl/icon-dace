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

%%{
 machine parse_intlist;
 access data->;

 # append to the buffer
 action append { if ((data->buflen < MAX_BUF_LEN) && (fc!='\n')) {
	  data->buffer[data->buflen++] = fc; }
 }
 # clear out the buffer
 action clear       { data->buflen = 0; }
 # push an integer value onto stack
 action store_value { 
    /* buffer string with "\0" character appended. */
    if (data->buflen < MAX_BUF_LEN) data->buffer[data->buflen++] = 0; 
    data->vstack[data->vidx++] = atol(data->buffer); 
 }
 # 
 action eval_expr { 
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
 action store_nlev { data->vstack[data->vidx++] = data->nlev_val; }
 # push an operator onto stack
 action store_op { 
   switch (fc) {
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
 # mark a single value
 action add_val { 
   int value = data->vstack[--data->vidx];
   if (value > data->maxval) {data->cs = parse_intlist_error;}
   data->out_values[value] = 1;
 }
 # mark a whole range
 action add_range { 
   int i,
     range_end   = data->vstack[--data->vidx],
     range_start = data->vstack[--data->vidx];
   if ((range_end   > data->maxval)  ||
       (range_start > data->maxval)  ||
       (range_start > range_end)) {data->cs = parse_intlist_error;}
   for (i=range_start; i<=range_end; ++i)  data->out_values[i] = 1;
 }


 # ===========================================================================
 # parser definitions
 # ===========================================================================
 
 # token for a single integer value
 value           = space* [1-9] >clear $append [0-9]* $append %store_value space*;
 # token for "nlev"
 nlev            = space* ("nlev"|"n"|"N") space* %store_nlev ;
 # token for constant
 constant        = value|nlev ;
 # tokens for basic arithmetic expressions
 operator        = space* ('+'|'-') $store_op space* ;
 # token for integer expression, allowing for parentheses
 group           = constant (operator constant %eval_expr)* ;
 par_group       = group|('(' group ')' );
 expression      = par_group (operator par_group %eval_expr)* ;
 # token for a range of numbers
 rangesymbol     = space* "..." space* ;
 range           = expression >clear  rangesymbol  expression >clear ;  
 # token for separator
 separator       = space* [;,] space* ;

 # namelist tokens:
 main := |*
   expression => add_val;
   range      => add_range ;
   separator;
 *|;
 # ===========================================================================
}%%

%% write data;


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
  %% write init;

  // execute Finite State Machine
  %% write exec;

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

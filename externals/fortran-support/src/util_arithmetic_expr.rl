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
 * Evaluation of basic arithmetic infix expressions based on a           *
 * Finite State Machine (FSM) and Dijkstra's shunting yard algorithm.    *
 * This file must be processed with Ragel to produce the final C code.   *
 *                                                                       *
 * Ragel State Machine Compiler version 7.0.4 was used                   *
 * --------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util_arithmetic_expr.h"

const char*  const priorities = "(;><;+-;/*;^";
const char*  const left_assoc = "><+-/*^";
const char * const fct_name[NUM_FCT] = { "exp", "log", "sin", "cos", "min", "max", 
					 "if", "sqrt", "erf" };

/* All data needed to define the state of the Finite State Machine parser. */
struct t_parsedata {
  int ierr;                                /* Error handler                 */
  int  cs, act, buflen;                    /* FSM state, buffer length      */
  char *ts, *te;                           /* token start, end              */
  char buffer[MAX_BUF_LEN];                /* token string buffer           */
  struct t_list ostack, *queue;            /* operator stack, postfix queue */
};


/* Priority of an arithmetic operator */
int priority(char op) {
  const char* o = priorities;
  for (int val=0; (*o)!='\0'; o++)
    if ((*o) == ';')     val++;
    else if ((*o) == op) return val;
  return -1;
}

/* Associativity of an arithmetic operator */
int left_associative(char op) {
  const char* o = left_assoc;
  for (; (*o)!='\0'; o++) if ((*o) == op) return 1;
  return 0;
}


/* --------------------------------------------------------------------- *
 * DEFINITION OF RAGEL FINITE STATE MACHINE "parse_infix"
 * --------------------------------------------------------------------- */

%%{
 machine parse_infix;
 access data->;

 # append to the token buffer
 action append { if ((data->buflen < MAX_BUF_LEN) && (fc!='\n'))
   data->buffer[data->buflen++] = fc;
 }
 # clear out the token buffer
 action clear  { data->buflen = 0; }
 # push a value onto stack
 action store_value { 
   data->queue->list[data->queue->size].type = VALUE;
   if (data->buflen < MAX_BUF_LEN) data->buffer[data->buflen] = 0;	
   sscanf(data->buffer, "%lf", &data->queue->list[data->queue->size++].val);
 }
 # push a constant onto stack
 action store_constant { 
   if (data->buflen < MAX_BUF_LEN) data->buffer[data->buflen] = 0;	
   const char* constant = data->buffer;
   data->queue->list[data->queue->size].type = VALUE;
   if (strcmp(constant, "pi") == 0)
     data->queue->list[data->queue->size++].val = 4.*atan(1.);
 }
 # push a function onto stack
 action store_fct { 
   if (data->buflen < MAX_BUF_LEN) data->buffer[data->buflen] = 0;	
   const char* fct = data->buffer;
   struct t_item new_item =  { .type=FUNCTION };
   for (int i=0; i<NUM_FCT; i++) {
     char fct_i[5];
     data->ierr = get_fctname(i, fct_i);
     if (data->ierr != 1 && strcmp(fct, fct_i) == 0)
       { new_item.fct = i; data->ostack.list[data->ostack.size++] = new_item; break; }
   }
 }
# function argument separator
action store_sep {
  /* If the token is a function argument separator (i.e. a comma):
     Until the token at the top of the stack is a left parenthesis, pop
     operators off the stack onto the output queue. */
  while ((data->ostack.size>0) && (data->ostack.list[data->ostack.size-1].type==OPERATOR) && 
	 (data->ostack.list[data->ostack.size-1].op != '(')) 
    data->queue->list[data->queue->size++] = data->ostack.list[--data->ostack.size];
}
# store operator 
action store_op { 
  /*  While the stack is not empty and an operator is at the top and
      the operator at the top is higher (or, if left associative,
      equal) priority that the item then pop the operator on the top
      of the stack, add the popped operator to the queue. */
  char this_op = data->buffer[0];
  int lassoc   = left_associative(this_op);
  while (((data->ostack.size > 0) && (data->ostack.list[data->ostack.size-1].type==OPERATOR)) &&
	 ( ((lassoc==1) && (priority(this_op) <= priority(data->ostack.list[data->ostack.size-1].op)) ) ||
	   ((lassoc==0) && (priority(this_op) <  priority(data->ostack.list[data->ostack.size-1].op)) ) ) )
	data->queue->list[data->queue->size++] = data->ostack.list[--data->ostack.size];
  struct t_item new_item =  { .type=OPERATOR, .op=this_op };
  data->ostack.list[data->ostack.size++] = new_item;
 }
 # store parentheses
 action store_paren { 
  switch (data->buffer[0]) {
     case '(': {
         struct t_item new_item =  { .type=OPERATOR, .op='(' };
         data->ostack.list[data->ostack.size++] = new_item;
       }
       break;
  case ')':
       if (data->ostack.size > 0) {
         /* Until token at the top of the stack is a left parenthesis,
	    pop operators off the stack onto the output queue. */
         struct t_item top = data->ostack.list[--data->ostack.size];
         while ((data->ostack.size>0) && (top.type==OPERATOR) && (top.op != '(')) {
          data->queue->list[data->queue->size++] = top;
          top = data->ostack.list[--data->ostack.size];
         }
	 /* If the token at the top of the stack is a function token,
	    pop it onto the output queue. */
	 if (data->ostack.list[data->ostack.size-1].type==FUNCTION)
	   data->queue->list[data->queue->size++] = data->ostack.list[--data->ostack.size];
       } 
       break;        
  }
 }
 # push a function onto stack
 action store_variable { 
   data->queue->list[data->queue->size].type = VARIABLE;
   if (data->buflen < MAX_BUF_LEN) data->buffer[data->buflen] = 0;	
   strcpy(data->queue->list[data->queue->size++].field, data->buffer);
 }


 # ===========================================================================
 # parser definitions
 # ===========================================================================
 
 # token for a single value
 value           = space* ([0-9]|'-'[0-9]) >clear $append [0-9]* '.'? $append [0-9]* $append 'e'? $append [0-9]* $append space* ;
 # tokens for basic arithmetic expressions
 operator        = space* ('+'|'-'|'*'|'/'|'^'|'>'|'<') >clear $append space* ;
 parenthesis     = space* ('('|')') >clear $append space* ;
 func            = space* ('exp'|'log'|'sin'|'cos'|'min'|'max'|'if'|'sqrt'|'erf') >clear $append space*;
 constant        = space* ('pi'|"r") >clear $append space* ;
 # function argument separator
 sep             = space* ',' space* ;
 # reference to extern variable
 var_token       = space* '[' >clear [A-Za-z0-9_]+ $append ']' space*;

 # namelist tokens:
 main := |*
   value       => store_value;
   func        => store_fct;
   operator    => store_op;
   parenthesis => store_paren;
   sep         => store_sep;
   constant    => store_constant;
   var_token   => store_variable; 
 *|;
 # ===========================================================================
}%%

%% write data;


/* --------------------------------------------------------------------- */
int do_parse_infix(const char *in_parse_line, struct t_list *queue) {
  char *p0 = malloc(strlen (in_parse_line) + 1);     /* Space for length plus null */
  if (p0 == NULL) return -1;                         /* No memory                  */
  strcpy(p0, in_parse_line);                         /* Copy the characters        */

  char *p = p0, *pe = p0 + strlen(p0), *eof = pe ;   /* pointer to input state and end. */
  struct t_parsedata pdata;
  pdata.queue = queue;
  pdata.ostack.size = pdata.queue->size = 0;
  struct t_parsedata *data = &pdata;
  
  // initialize Finite State Machine
  %% write init;

  // execute Finite State Machine
  %% write exec;

  free(p0);
  if (data->cs == parse_infix_error) return 1;

  /* While the stack is not empty pop the stack, add to queue. */
  while (data->ostack.size > 0) 
    data->queue->list[data->queue->size++] = data->ostack.list[--data->ostack.size];
  return 0;
}


int get_fctname(const int id, char *string) {
  int ierr = 0;
  if ((id < 0) || (id >= NUM_FCT)) {
    ierr = 1;
  }
  else {
    strcpy(string, fct_name[id]);
  }
  return ierr;
}


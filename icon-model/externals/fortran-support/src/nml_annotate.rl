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
 * Fortran namelist scanner based on a Finite State Machine (FSM).       *
 * This file must be processed with Ragel to produce the final C code.   *
 *                                                                       *
 * Ragel State Machine Compiler version 7.0.4 was used                   *
 * --------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_BUF_LEN 32768
#define TRUNCATE_LEN 40

/* linked list data structure */
struct t_val {
    char *name;
    struct t_val *last, *next;
};
struct t_key {
    char *name;
    struct t_val *value, *dflt;
    struct t_key *last, *next;
};
struct t_nml {
    char *name;
    struct t_key *key_list;
    struct t_nml *last, *next;
};

/* All data needed to define the state of the Finite State Machine parser. */
struct t_nmldata {
    int cs, act, buflen;          /* FSM state, buffer length     */
    char *ts, *te;                /* token start, end             */
    char buffer[MAX_BUF_LEN + 1]; /* token string buffer          */
    struct t_nml *nml;            /* list of namelists            */

    struct t_nml *cur_nml; /* pointer to current list item */
    struct t_key *cur_key;
    struct t_val *cur_val;
};

/* return buffer string with "\0" character appended. */
char *terminate_string(struct t_nmldata *data) {
    if (data->buflen < MAX_BUF_LEN) data->buffer[data->buflen++] = 0;
    return data->buffer;
}

/* --------------------------------------------------------------------- *
 * Print-out of namelist data structure.
 * --------------------------------------------------------------------- */

/* traverse val list data structure and print out its contents; compare with
 * defaults */
void print_val(FILE *out, struct t_val *val, struct t_val *dflt,
               int *lchanged) {
    int lchanged_ = *lchanged;
    while (val) {
        fprintf(out, "%.*s%s%s%s%s", TRUNCATE_LEN, val->name,
                (strlen(val->name) >= TRUNCATE_LEN) ? " [...]" : "",
                (strlen(val->name) >= TRUNCATE_LEN)
                    ? (val->name + strlen(val->name) - 1)
                    : "",
                (strlen(val->name) >= TRUNCATE_LEN) ? " (truncated)" : "",
                (val->next != NULL) ? ", " : "\n");
        if (dflt != NULL) {
            lchanged_ |= (strcmp(val->name, dflt->name) != 0);
            dflt = dflt->next;
        }
        val = val->next;
    }
    *lchanged = lchanged_;
}
/* traverse key list and print out key-value pairs (and defaults). */
void print_key(FILE *out, struct t_key *key) {
    while (key) {
        fprintf(out, "    %-*s    ", TRUNCATE_LEN, key->name);
        int lchanged = 0;
        print_val(out, key->value, key->dflt, &lchanged);
        if (lchanged != 0) {
            lchanged = 0;
            fprintf(out, "        %*s>> DEFAULT: ", TRUNCATE_LEN, " ");
            print_val(out, key->dflt, key->dflt, &lchanged);
        }
        key = key->next;
    }
}
/* traverse namelist data structure and print out its contents. */
void print_nml(FILE *out, struct t_nml *nml) {
    while (nml) {
        fprintf(out, "\nNAMELIST %s\n", nml->name);
        print_key(out, nml->key_list);
        nml = nml->next;
    }
}

/* --------------------------------------------------------------------- *
 * Removing elements from namelist data structure.
 * --------------------------------------------------------------------- */

/* remove val list from data structure.                         */
int delete_val(struct t_val **val) {
    if ((*val) == NULL) return 0;
    struct t_val *tmp = (*val)->next;
    if ((*val)->last != NULL) (*val)->last->next = (*val)->next;
    if ((*val)->next != NULL) (*val)->next->last = (*val)->last;
    free((*val)->name);
    free(*val);
    (*val) = tmp;
    return 1;
}
/* remove key list from data structure.                         */
int delete_key(struct t_key **key) {
    if ((*key) == NULL) return 0;
    struct t_key *tmp = (*key)->next;
    if ((*key)->last != NULL) (*key)->last->next = (*key)->next;
    if ((*key)->next != NULL) (*key)->next->last = (*key)->last;
    while (delete_val(&(*key)->value))
        ; /* key: delete values   */
    while (delete_val(&(*key)->dflt))
        ; /* key: delete defaults */
    free((*key)->name);
    free(*key);
    (*key) = tmp;
    return 1;
}
/* remove namelist from data structure.                         */
int delete_nml(struct t_nml **nml) {
    if ((*nml) == NULL) return 0;
    struct t_nml *tmp = (*nml)->next;
    if ((*nml)->last != NULL) (*nml)->last->next = (*nml)->next;
    if ((*nml)->next != NULL) (*nml)->next->last = (*nml)->last;
    while (delete_key(&(*nml)->key_list))
        ;
    free((*nml)->name);
    free(*nml);
    (*nml) = tmp;
    return 1;
}

/* --------------------------------------------------------------------- *
 * Definition of default values.
 * --------------------------------------------------------------------- */

/* Define every odd namelist to be the defaults of the next namelist item.
   However, if the namelists occurs only once, we skip it. */
int define_odd_namelists_as_defaults(struct t_nml **nml) {
    struct t_nml *defaults = *nml;

    int i = 0;
    while (defaults != NULL) {
        if (defaults->next != NULL) {
            struct t_nml *settings = defaults->next;
            struct t_key *defaults_key = defaults->key_list,
                         *settings_key = settings->key_list;

            if (strcmp(defaults->name, settings->name) == 0) {
                while (defaults_key != NULL) {
                    if ((settings_key == NULL)
                        || (strcmp(settings_key->name, defaults_key->name)
                            != 0)) {
                        fprintf(stderr, "Keyword mismatch! %s / %s\n",
                                settings_key->name, defaults_key->name);
                        return 0;
                    } else if ((settings_key->dflt != NULL)
                               || (defaults_key->value == NULL)) {
                        fprintf(stderr, "Internal error!\n");
                        return 0;
                    } else {
                        settings_key->dflt = defaults_key->value;
                        defaults_key->value = NULL;
                    } /* else */
                    defaults_key = defaults_key->next;
                    settings_key = settings_key->next;
                } /* while */
                delete_nml(&defaults);
                if (i == 0) *nml = defaults;
            } /* if */
        }     /* if */
        defaults = defaults->next;
        i++;
    } /* while */
    return 1;
}

/* --------------------------------------------------------------------- *
 * DEFINITION OF RAGEL FINITE STATE MACHINE
 * --------------------------------------------------------------------- */

%%{
 machine nml;
 access data->;

 # append to the buffer
 action append { if ((data->buflen < MAX_BUF_LEN) && (fc!='\n')) 
                   data->buffer[data->buflen++] = fc; }
 # clear out the buffer
 action clear  { data->buflen = 0; }
 
 # append element to doubly linked list of namelists
 action add_nml { 
  struct t_nml *new_nml = (struct t_nml*) calloc(1, sizeof(struct t_nml));
  *new_nml = (struct t_nml) { .name = strdup(terminate_string(data)), 
                              .next=NULL, .last=data->cur_nml };
  if (data->cur_nml == NULL) data->nml           = data->cur_nml = new_nml;
  else                       data->cur_nml = data->cur_nml->next = new_nml;
  data->cur_key = NULL;
 }

 # append element to doubly linked list of keys
 action add_key {
  struct t_key *new_key = (struct t_key*) calloc(1, sizeof(struct t_key));
  *new_key = (struct t_key) { .name = strdup(terminate_string(data)), 
                              .next=NULL, .last=data->cur_key, .dflt = NULL };
  if (data->cur_key == NULL) data->cur_nml->key_list = data->cur_key = new_key;
  else                       data->cur_key = data->cur_key->next     = new_key;
  data->cur_val = NULL;
 }

 # append element to doubly linked list of values
 action add_val { 
  struct t_val *new_val = (struct t_val*) calloc(1, sizeof(struct t_val));
  *new_val = (struct t_val) { .name = strdup(terminate_string(data)), 
                              .next=NULL, .last=data->cur_val };
  if (data->cur_val == NULL) data->cur_key->value = data->cur_val = new_val;
  else                       data->cur_val = data->cur_val->next  = new_val;
 }

 # valid namelist key identifier (allowing also array subscripts xyz(...))
 identifier      = [a-zA-Z_] >clear $append [a-zA-Z_0-9()%]** $append ;
 # character for valid right hand side value (if not a string)
 value           = [a-zA-Z0-9_.+\-*] >clear $append [a-zA-Z0-9_.+\-*]* $append ;
 # string, delimited by " or ' (with preceding quantifier "<number>*")
 string          = ([0-9]* '*')? ('\'' [^']* '\'')|('"' [^"]* '"') ;

 # namelist name
 nml_name        = space* '&' identifier space+ ;
 # valid namelist key (followed by '=' sign)
 nml_key         = space* identifier space* '=' space* ;
 # right hand side namelist value (followed by comma and/or '/')
 nml_value       = space* (value|string >clear $append) space* (space+|','|(','? space* '/')) ;

 # namelist tokens:
 main := |*
          nml_name     => add_nml;   # namelist name
          nml_key      => add_key;   # a single key
          nml_value    => add_val;   # a single namelist value
         *|;
}%%

%% write data;


/* --------------------------------------------------------------------- *
 * Main routine: Read namelists from file, write annotated namelists.
 * --------------------------------------------------------------------- */

int util_annotate_nml(char *in_filename, char *out_filename) {
    FILE *nmlfile;                        /* input file object      */
    char buf[MAX_BUF_LEN + 1] = { "\0" }; /* input string buffer    */
    char *p = buf, *eof = 0;              /* input start, end       */
    struct t_nmldata nmldata;             /* input scanner state    */
    struct t_nmldata *data = &nmldata;
    data->nml = data->cur_nml = NULL; /* initialize linked list */

    if ((nmlfile = fopen(in_filename, "r")) == NULL) {
        fprintf(stderr, "Could not open file.\n");
        return -1;
    }

    // initialize Finite State Machine
    %% write init;

    /* process file line-by-line */
    int have = 0;
    while (fgets(p, MAX_BUF_LEN - have, nmlfile) != NULL) {
        char *pe = buf + strlen(buf); /* pointer to input end. */

        // execute Finite State Machine
        %% write exec;

        if (data->cs == nml_error) {
            fprintf(stderr, "nml_annotate: could not translate the collected "
                         "namelist data into a table!\n");
            break;
        }

        have = 0;
        if (data->ts != 0) {
            /* wrap around remaining part ("have" chars) of the buffer. */
            have = pe - data->ts;
            memmove(buf, data->ts, have);
            data->te = buf + (data->te - data->ts); /* pointer to token end.   */
            data->ts = buf;                         /* pointer to token start. */
        }
        p = buf + have;
    }
    fclose(nmlfile);

    /* use the contents of the namelists with odd indices as the defaults
    of their succeeding namelists; remove the default namelists
    afterwards. */
    define_odd_namelists_as_defaults(&data->nml);

    /* print out the contents of all namelists (with defaults). */
    FILE *out = fopen(out_filename, "w");
    print_nml(out, data->nml);
    fclose(out);

    /* clean up */
    while (delete_nml(&data->nml))
        ;

    return 0;
} /* int util_annotate_nml */

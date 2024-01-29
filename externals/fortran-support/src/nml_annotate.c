
#line 1 "nml_annotate.rl"
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


#line 260 "nml_annotate.rl"



#line 209 "nml_annotate.c"
static const char _nml_actions[] = {
	0, 1, 0, 1, 2, 1, 3, 1, 
	4, 1, 5, 1, 6, 1, 7, 1, 
	8, 1, 9, 2, 1, 0
};

static const unsigned char _nml_key_offsets[] = {
	0, 0, 17, 18, 23, 27, 32, 45, 
	46, 59, 71, 86, 104, 118, 122, 139, 
	144, 148, 151, 157
};

static const char _nml_trans_keys[] = {
	32, 34, 38, 39, 42, 43, 95, 9, 
	13, 45, 46, 48, 57, 65, 90, 97, 
	122, 34, 32, 44, 47, 9, 13, 32, 
	47, 9, 13, 95, 65, 90, 97, 122, 
	32, 37, 95, 9, 13, 40, 41, 48, 
	57, 65, 90, 97, 122, 39, 32, 39, 
	44, 47, 95, 9, 13, 42, 57, 65, 
	90, 97, 122, 32, 44, 47, 95, 9, 
	13, 42, 57, 65, 90, 97, 122, 32, 
	42, 44, 47, 95, 9, 13, 43, 46, 
	48, 57, 65, 90, 97, 122, 32, 37, 
	44, 47, 61, 95, 9, 13, 40, 41, 
	42, 46, 48, 57, 65, 90, 97, 122, 
	32, 37, 61, 95, 9, 13, 40, 41, 
	48, 57, 65, 90, 97, 122, 32, 61, 
	9, 13, 32, 34, 38, 39, 42, 43, 
	95, 9, 13, 45, 46, 48, 57, 65, 
	90, 97, 122, 32, 44, 47, 9, 13, 
	32, 47, 9, 13, 32, 9, 13, 32, 
	44, 47, 61, 9, 13, 32, 9, 13, 
	0
};

static const char _nml_single_lengths[] = {
	0, 7, 1, 3, 2, 1, 3, 1, 
	5, 4, 5, 6, 4, 2, 7, 3, 
	2, 1, 4, 1
};

static const char _nml_range_lengths[] = {
	0, 5, 0, 1, 1, 2, 5, 0, 
	4, 4, 5, 6, 5, 1, 5, 1, 
	1, 1, 1, 1
};

static const unsigned char _nml_index_offsets[] = {
	0, 0, 13, 15, 20, 24, 28, 37, 
	39, 49, 58, 69, 82, 92, 96, 109, 
	114, 118, 121, 127
};

static const char _nml_indicies[] = {
	0, 2, 3, 4, 5, 6, 8, 0, 
	6, 7, 8, 8, 1, 10, 9, 11, 
	12, 13, 11, 1, 15, 13, 15, 14, 
	16, 16, 16, 1, 17, 18, 18, 17, 
	18, 18, 18, 18, 1, 10, 19, 11, 
	19, 12, 13, 20, 11, 20, 20, 20, 
	1, 11, 12, 13, 20, 11, 20, 20, 
	20, 1, 11, 21, 12, 13, 20, 11, 
	20, 22, 20, 20, 1, 23, 24, 12, 
	13, 26, 25, 23, 24, 20, 25, 25, 
	25, 1, 27, 24, 26, 24, 27, 24, 
	24, 24, 24, 1, 27, 26, 27, 1, 
	0, 2, 3, 4, 5, 6, 8, 0, 
	6, 7, 8, 8, 1, 11, 12, 13, 
	11, 28, 15, 13, 15, 28, 17, 17, 
	29, 23, 12, 13, 26, 23, 28, 26, 
	26, 30, 0
};

static const char _nml_trans_targs[] = {
	1, 0, 2, 5, 7, 8, 9, 10, 
	11, 2, 3, 15, 16, 14, 14, 4, 
	6, 17, 6, 7, 9, 8, 10, 18, 
	12, 11, 19, 13, 14, 14, 14
};

static const char _nml_trans_actions[] = {
	0, 0, 19, 0, 19, 19, 19, 19, 
	19, 1, 1, 0, 7, 9, 17, 0, 
	19, 0, 1, 1, 1, 1, 1, 0, 
	1, 1, 0, 0, 15, 11, 13
};

static const char _nml_to_state_actions[] = {
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 3, 0, 
	0, 0, 0, 0
};

static const char _nml_from_state_actions[] = {
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 5, 0, 
	0, 0, 0, 0
};

static const unsigned char _nml_eof_trans[] = {
	0, 0, 0, 0, 15, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 29, 
	29, 30, 29, 31
};

static const int nml_start = 14;
static const int nml_first_final = 14;
static const int nml_error = 0;

static const int nml_en_main = 14;


#line 263 "nml_annotate.rl"


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
    
#line 345 "nml_annotate.c"
	{
	 data->cs = nml_start;
	 data->ts = 0;
	 data->te = 0;
	 data->act = 0;
	}

#line 284 "nml_annotate.rl"

    /* process file line-by-line */
    int have = 0;
    while (fgets(p, MAX_BUF_LEN - have, nmlfile) != NULL) {
        char *pe = buf + strlen(buf); /* pointer to input end. */

        // execute Finite State Machine
        
#line 362 "nml_annotate.c"
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
	_acts = _nml_actions + _nml_from_state_actions[ data->cs];
	_nacts = (unsigned int) *_acts++;
	while ( _nacts-- > 0 ) {
		switch ( *_acts++ ) {
	case 3:
#line 1 "NONE"
	{ data->ts = p;}
	break;
#line 383 "nml_annotate.c"
		}
	}

	_keys = _nml_trans_keys + _nml_key_offsets[ data->cs];
	_trans = _nml_index_offsets[ data->cs];

	_klen = _nml_single_lengths[ data->cs];
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

	_klen = _nml_range_lengths[ data->cs];
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
	_trans = _nml_indicies[_trans];
_eof_trans:
	 data->cs = _nml_trans_targs[_trans];

	if ( _nml_trans_actions[_trans] == 0 )
		goto _again;

	_acts = _nml_actions + _nml_trans_actions[_trans];
	_nacts = (unsigned int) *_acts++;
	while ( _nacts-- > 0 )
	{
		switch ( *_acts++ )
		{
	case 0:
#line 206 "nml_annotate.rl"
	{ if ((data->buflen < MAX_BUF_LEN) && ((*p)!='\n')) 
                   data->buffer[data->buflen++] = (*p); }
	break;
	case 1:
#line 209 "nml_annotate.rl"
	{ data->buflen = 0; }
	break;
	case 4:
#line 1 "NONE"
	{ data->te = p+1;}
	break;
	case 5:
#line 232 "nml_annotate.rl"
	{ data->te = p+1;{ 
  struct t_val *new_val = (struct t_val*) calloc(1, sizeof(struct t_val));
  *new_val = (struct t_val) { .name = strdup(terminate_string(data)), 
                              .next=NULL, .last=data->cur_val };
  if (data->cur_val == NULL) data->cur_key->value = data->cur_val = new_val;
  else                       data->cur_val = data->cur_val->next  = new_val;
 }}
	break;
	case 6:
#line 212 "nml_annotate.rl"
	{ data->te = p;p--;{ 
  struct t_nml *new_nml = (struct t_nml*) calloc(1, sizeof(struct t_nml));
  *new_nml = (struct t_nml) { .name = strdup(terminate_string(data)), 
                              .next=NULL, .last=data->cur_nml };
  if (data->cur_nml == NULL) data->nml           = data->cur_nml = new_nml;
  else                       data->cur_nml = data->cur_nml->next = new_nml;
  data->cur_key = NULL;
 }}
	break;
	case 7:
#line 222 "nml_annotate.rl"
	{ data->te = p;p--;{
  struct t_key *new_key = (struct t_key*) calloc(1, sizeof(struct t_key));
  *new_key = (struct t_key) { .name = strdup(terminate_string(data)), 
                              .next=NULL, .last=data->cur_key, .dflt = NULL };
  if (data->cur_key == NULL) data->cur_nml->key_list = data->cur_key = new_key;
  else                       data->cur_key = data->cur_key->next     = new_key;
  data->cur_val = NULL;
 }}
	break;
	case 8:
#line 232 "nml_annotate.rl"
	{ data->te = p;p--;{ 
  struct t_val *new_val = (struct t_val*) calloc(1, sizeof(struct t_val));
  *new_val = (struct t_val) { .name = strdup(terminate_string(data)), 
                              .next=NULL, .last=data->cur_val };
  if (data->cur_val == NULL) data->cur_key->value = data->cur_val = new_val;
  else                       data->cur_val = data->cur_val->next  = new_val;
 }}
	break;
	case 9:
#line 232 "nml_annotate.rl"
	{{p = (( data->te))-1;}{ 
  struct t_val *new_val = (struct t_val*) calloc(1, sizeof(struct t_val));
  *new_val = (struct t_val) { .name = strdup(terminate_string(data)), 
                              .next=NULL, .last=data->cur_val };
  if (data->cur_val == NULL) data->cur_key->value = data->cur_val = new_val;
  else                       data->cur_val = data->cur_val->next  = new_val;
 }}
	break;
#line 514 "nml_annotate.c"
		}
	}

_again:
	_acts = _nml_actions + _nml_to_state_actions[ data->cs];
	_nacts = (unsigned int) *_acts++;
	while ( _nacts-- > 0 ) {
		switch ( *_acts++ ) {
	case 2:
#line 1 "NONE"
	{ data->ts = 0;}
	break;
#line 527 "nml_annotate.c"
		}
	}

	if (  data->cs == 0 )
		goto _out;
	if ( ++p != pe )
		goto _resume;
	_test_eof: {}
	if ( p == eof )
	{
	if ( _nml_eof_trans[ data->cs] > 0 ) {
		_trans = _nml_eof_trans[ data->cs] - 1;
		goto _eof_trans;
	}
	}

	_out: {}
	}

#line 292 "nml_annotate.rl"

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

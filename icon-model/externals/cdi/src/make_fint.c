#ifdef HAVE_CONFIG_H
#include "config.h"
#if defined(HAVE_GETLINE) && !defined(_GNU_SOURCES)
// It seems that _GNU_SOURCES must be defined to enable getline(3) on some
// systems. However, it is unclear whether it is possible at all that
// HAVE_GETLINE gets defined on such systems by the configure script, which does
// not expand either AC_GNU_SOURCE or AC_USE_SYSTEM_EXTENSIONS and, therefore,
// checks for getline without _GNU_SOURCES defined.
#define _GNU_SOURCES 1
#endif
#else
#define VERSION "2.3.0"
#endif

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#endif

#include <ctype.h>
#include <errno.h>
#include <regex.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#ifndef HAVE_GETLINE
#ifdef HAVE_CONFIG_H
// Avoid warning: implicit declaration of function ‘getline’
ssize_t getline(char **linebuf, size_t *linebuf_size, FILE *fp);
#else
// This is only to cover the case when this file is compiled "manually", i.e.
// without the Makefile and without config.h, e.g. 'gcc make_fint.c', on a
// system where HAVE_GETLINE would not get defined or getline does not work
// without _GNU_SOURCES defined.
#include "getline.c"
#endif
#endif

// clang-format off

typedef struct
{
  size_t naline;
  char *fname;
  char *aline[99];
  char *text;
}
Docu;

// Example: ./make_fint -d ../doc/pio/ cdipio.h

static struct {
  const char *name;
  size_t len;
} fname_list[] = {
  { "c_quick_ref.txt", 0 },
  { "f_quick_ref.txt", 0 },
  { "tex/c_quick_ref.tex", 0 },
  { "tex/f_quick_ref.tex", 0 },
};
enum {
  NAME_C_QUICK_REF,
  NAME_F_QUICK_REF,
  NAME_C_QUICK_REF_TEX,
  NAME_F_QUICK_REF_TEX,
  fname_list_size = sizeof(fname_list)/sizeof(fname_list[0]),
};


static Docu cdoc[9999], fdoc[9999];
static size_t ncdoc = 0, nfdoc = 0;
static int debug = 0, verbose = 0;

static int doccmp(const void *s1, const void *s2)
{
  Docu *x = (Docu *) s1;
  Docu *y = (Docu *) s2;

  return strcmp(x->fname, y->fname);
}

static void doctotex(FILE *fp, const Docu *doc, const size_t ndoc)
{
  for (size_t i = 0; i < ndoc; ++i)
    {
      fprintf(fp, "\\section*{\\tt \\htmlref{%s}{%s}}\n\n", doc[i].fname, doc[i].fname);
      fprintf(fp, "\\begin{verbatim}\n");
      for (size_t k = 0; k < doc[i].naline; ++k)
	fprintf(fp, "    %s\n", doc[i].aline[k]);
      fprintf(fp, "\\end{verbatim}\n");
      fprintf(fp, "\n%s.\n\n\n", doc[i].text);
    }
}

static void doctotxt(FILE *fp, const Docu *doc, const size_t ndoc)
{
  for (size_t i = 0; i < ndoc; ++i)
    {
      fprintf(fp, "%s\n\n", doc[i].fname);
      for (size_t k = 0; k < doc[i].naline; ++k)
	fprintf(fp, "    %s\n", doc[i].aline[k]);
      fprintf(fp, "\n  %s.\n\n", doc[i].text);
    }
}

enum cftype {ISVOID, ISCONSTSTRING, ISINT, ISLOGICAL, ISREAL, ISDOUBLE, ISDATETYPE, ISSIZETYPE, ISMPI_COMM,
             ISXT_IDXLIST, ISXT_IDXLISTV, ISCHOICE, ISINTP, ISDATETYPEP, ISSIZETYPEP, ISFLOATV, ISFLOATVV,
             ISDOUBLEV, ISDOUBLEVV, ISINTV, ISINTVV, ISINTVVV, ISREALP,
             ISDOUBLEP, ISCBUF, ISUUID, ISUCHAR, ISSTRING, ISSTRINGP,
             VOIDFUNCVOID,
             NUM_KNOWN_ARG_TYPES};

static inline int
isArrayArgType(int argType);

enum conversionType { CONV_ARG, CONV_RET };


typedef int (*cfConversionEmitter)(FILE *outfp, const char *argName, size_t argNameLen, enum conversionType part);
typedef int (*cfPrologueEmitter)(FILE *outfp, size_t argNum);


static int cfMPICommConvert(FILE *outfp, const char *argName, size_t argNameLen, enum conversionType part);
static int cfDateTypeConvert(FILE *outfp, const char *argName, size_t argNameLen, enum conversionType part);
static int cfSizeTypeConvert(FILE *outfp, const char *argName, size_t argNameLen, enum conversionType part);
static int cfDateTypepConvert(FILE *outfp, const char *argName, size_t argNameLen, enum conversionType part);
static int cfSizeTypepConvert(FILE *outfp, const char *argName, size_t argNameLen, enum conversionType part);
static int cfXtIdxlistConvert(FILE *outfp, const char *argName, size_t argNameLen, enum conversionType part);
static int cfVoidFuncPrologue(FILE *outfp, size_t argNum);

struct symbol {
  const char *f77name, *cfint, *cfmt, *parseRE;
  // pair of parentheses which matches the argument name
  size_t nameMatch;
  bool needsExtraWrapper, needsPrologue;
  cfConversionEmitter convert;
  const char *convcfmt;
  cfPrologueEmitter prologue;
  regex_t preg;
};

// C symbol names
#define SYMRE "([A-Za-z_][A-Za-z_0-9]*)"
static inline int isSymStart(const int c)
{
  return (isalpha(c) || c == '_');
}

static inline int isSym(const int c)
{
  return (isalnum(c) || c == '_');
}

/* white-space */
#define WS "[[:blank:]\n]"
#define NWS "[^[:blank:]\n]"
#define ARRAY_BOUND "\\[([^]]*)\\]"WS"*"

#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 5) || defined (__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#endif
/* Note: size of this table must match the cftype enum */
static struct symbol funArgSym[]
  = { { "",                "",        "%svoid",
        "^"WS"*void"WS"*\\)", 0, 0, 0 },
      { "CHARACTER(80)",    "STRING",  "%schar *%.*s",
        "^"WS"*const"WS"+char"WS"+\\*"SYMRE WS"*\\(", 1, false, false },
      { "INTEGER",         "INT",     "%sint %.*s",
        "^"WS"*(const"WS"+)?int("WS"+"SYMRE")?"WS"*[,\\)]", 3, false, false },
      { "LOGICAL",         "LOGICAL", "%bool %.*s",
        "^"WS"*(const"WS"+)?bool("WS"+"SYMRE")?"WS"*[,\\)]", 3, false, false },
      { "REAL",            "FLOAT",   "%sfloat %.*s",
        "^"WS"*(const"WS"+)?float"WS"+"SYMRE"?"WS"*[,\\)]", 2, false, false },
      { "DOUBLEPRECISION", "DOUBLE",  "%sdouble %.*s",
        "^"WS"*(const"WS"+)?double"WS"+"SYMRE"?"WS"*[,\\)]", 2, false, false },
      { "INTEGER",         "INT",     "%sDateType %.*s",
        "^"WS"*(const"WS"+)?DateType("WS"+"SYMRE")?"WS"*[,\\)]", 3, true, false,
        cfDateTypeConvert, "%sint %.*s" },
      { "INTEGER",         "INT",     "%sSizeType %.*s",
        "^"WS"*(const"WS"+)?SizeType("WS"+"SYMRE")?"WS"*[,\\)]", 3, true, false,
        cfSizeTypeConvert, "%sint %.*s" },
      { "INTEGER",         "INT", "%sMPI_Comm %.*s",
        "^"WS"*MPI_Comm"WS"+"SYMRE"?"WS"*[,\\)]", 1, true, false,
        cfMPICommConvert, "%sint %.*s" },
      { "TYPE(XT_IDXLIST)", "PVOID", "%sXt_idxlist %.*s",
        "^"WS"*Xt_idxlist"WS"+"SYMRE"?"WS"*[,\\)]", 1, true, false,
        cfXtIdxlistConvert, "%svoid *%.*s" },
      { "TYPE(XT_IDXLIST)", "PVOID", "%sXt_idxlist %.*s[]",
        "^"WS"*(const"WS"+)?Xt_idxlist("WS"+"SYMRE")?"WS"*"ARRAY_BOUND
        "[,\\)]", 3, false, false },
      { "CHOICE", "PVOID", "%sconst void *%.*s",
        "^"WS"*const"WS"+void"WS"*\\*"WS"*"SYMRE"?"WS"*[,\\)]", 1, false, false },
      { "INTEGER",         "PINT",    "%sint *%.*s",
        "^"WS"*(const"WS"+)?int"WS"+\\*"SYMRE"?"WS"*[,\\)]", 2, false, false },
      { "INTEGER",         "PINT",    "%sDateType *%.*s",
        "^"WS"*(const"WS"+)?DateType"WS"+\\*"SYMRE"?"WS"*[,\\)]", 2, true, false,
        cfDateTypepConvert, "%sint *%.*s"},
      { "INTEGER",         "PINT",    "%sSizeType *%.*s",
        "^"WS"*(const"WS"+)?SizeType"WS"+\\*"SYMRE"?"WS"*[,\\)]", 2, true, false,
        cfSizeTypepConvert, "%sint *%.*s"},
      { "REAL",            "FLOATV",  "%sfloat %.*s[]",
        "^"WS"*(const"WS"+)?float("WS"+"SYMRE")?"WS"*"ARRAY_BOUND
        "[,\\)]", 3, false, false },
      { "REAL",            "FLOATVV", "%sfloat %.*s%s",
        "^"WS"*(const"WS"+)?float("WS"+"SYMRE")?"WS"*"ARRAY_BOUND
        ARRAY_BOUND"[,\\)]", 3, false, false },
      { "DOUBLEPRECISION", "DOUBLEV",  "%sdouble %.*s[]",
        "^"WS"*(const"WS"+)?double("WS"+"SYMRE")?"WS"*"ARRAY_BOUND
        "[,\\)]", 3, false, false },
      { "DOUBLEPRECISION", "DOUBLEVV", "%sdouble %.*s%s",
        "^"WS"*(const"WS"+)?double("WS"+"SYMRE")?"WS"*"ARRAY_BOUND
        ARRAY_BOUND"[,\\)]", 3, false, false },
      { "INTEGER",         "INTV",    "%sint  %.*s[]",
        "^"WS"*(const"WS"+)?int("WS"+"SYMRE")?"WS"*"ARRAY_BOUND
        "[,\\)]", 3, false, false },
      { "INTEGER",         "INTVV",    "%sint %.*s%s",
        "^"WS"*(const"WS"+)?int("WS"+"SYMRE")?"WS"*"ARRAY_BOUND
        ARRAY_BOUND "[,\\)]", 3, false, false },
      { "INTEGER",         "INTVVV",    "%sint %.*s%s",
        "^"WS"*(const"WS"+)?int("WS"+"SYMRE")?"WS"*"ARRAY_BOUND
        ARRAY_BOUND ARRAY_BOUND"[,\\)]", 3, false, false },
      { "REAL",            "PFLOAT",  "%sfloat *%.*s",
        "^"WS"*(const"WS"+)?float"WS"+\\*"SYMRE"?"WS"*[,\\)]", 2, false, false },
      { "DOUBLEPRECISION", "PDOUBLE", "%sdouble *%.*s",
        "^"WS"*(const"WS"+)?double"WS"+\\*"SYMRE"?"WS"*[,\\)]", 2, false, false },
      { "CHARACTER*(*)",   "PPSTRING",    "%schar *%.*s",
        "^"WS"*(const"WS"+)?char"WS"+\\*""([A-Za-z_][A-Za-z_0-9]*_cbuf)"
        WS"*[,\\)]", 2, false, false },
      { "INTEGER*1(16)",   "PVOID",    "%sunsigned char %.*s[16]",
        "^"WS"*(const"WS"+)?unsigned"WS"+char"WS"+"SYMRE"?\\[(16|CDI_UUID_SIZE)\\]"WS"*[,\\)]", 2, false, false },
      { "INTEGER*1(*)",   "PVOID",    "%sunsigned char *%.*s",
        "^"WS"*(const"WS"+)?unsigned"WS"+char"WS"+\\*"SYMRE"?"WS"*[,\\)]", 2, false, false },
      { "CHARACTER*(*)",   "STRING",  "%schar *%.*s",
        "^"WS"*const"WS"+char"WS"+\\*"WS"*"SYMRE"?"WS"*[,\\)]", 1, false, false },
      { "CHARACTER*(*)",   "PSTRING", "%schar *%.*s",
        "^"WS"*char"WS"+\\*"SYMRE"?"WS"*[,\\)]", 1, false, false },
      { "PROCEDURE", "ROUTINE", "%svoid (*%.*s)(void)",
        "^"WS"*void"WS"*\\("WS"*\\*"WS"*"SYMRE"?"WS"*\\)"
        WS"*\\("WS"*void"WS"*\\)"WS"*[,\\)]", 1, false, true,
        NULL, NULL, cfVoidFuncPrologue },
};

static struct symbol funRet[] = {
  { "",                "",        "%svoid %.*s",
    "void"WS"+"SYMRE WS"*\\(", 1, false, false },
  { "CHARACTER",       "STRING",  "%schar *%.*s",
    "char"WS"+\\*"WS"*"SYMRE WS"*\\(", 1, false, false },
  { "INTEGER",         "INT",     "%sint %.*s",
    "(const"WS"+)?int"WS"+"SYMRE WS"*\\(", 2, false, false },
  { "LOGICAL",         "LOGICAL", "%sbool %.*s",
    "(const"WS"+)?bool"WS"+"SYMRE WS"*\\(", 2, false, false },
  { "REAL",            "FLOAT",   "%sfloat %.*s",
    "(const"WS"+)?float"WS"+"SYMRE WS"*\\(", 2, false, false },
  { "DOUBLEPRECISION", "DOUBLE",  "%sdouble %.*s",
    "(const"WS"+)?double"WS"+"SYMRE WS"*\\(", 2, false, false },
  { "INTEGER",         "INT",     "%sMPI_Comm %.*s",
    "MPI_Comm"WS"+"SYMRE WS"*\\(", 1, true, false, cfMPICommConvert, "%sint %.*s" },
  { "INTEGER",         "INT",    "%sDateType %.*s",
    "(const"WS"+)?DateType"WS"+"SYMRE WS"*\\(", 2, true, false, cfDateTypeConvert, "%sint %.*s" },
  { "INTEGER",         "INT",    "%sSizeType %.*s",
    "(const"WS"+)?SizeType"WS"+"SYMRE WS"*\\(", 2, true, false, cfSizeTypeConvert, "%sint %.*s" },
};
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 5)
#pragma GCC diagnostic pop
#endif

enum { NUM_RET_TYPES = sizeof (funRet) / sizeof (funRet[0]) };
enum decl { UNKNOWN_DECL, FUNC_DECL, PARAM_DECL };

enum {
  MAX_FUNC_ARGS = 200,
  MAX_FUNC_NAME_LEN = 127,
};

static inline size_t
compress_whitespace(size_t len, char str[]);

static int
reCompile(regex_t *restrict RE, const char *restrict REstring,
          char * restrict *restrict lineBuf, size_t * restrict lineBufSize);

static size_t
symRegexCompile(size_t numSyms, struct symbol symList[],
                char **line, size_t *lineBufSize);

static void
build_header_name(size_t len, const char *fname, char *cppMacro);

static int detectComment(char **line_, ssize_t *lineLen, size_t *lineBufSize,
                         size_t maxMatch, regmatch_t reMatch[],
                         char *xname, size_t *xnameLen, char **xdes,
                         int *lineno, FILE *fpin, FILE *fpinc, FILE *fpint);

static regex_t commentStartRE, commentEndRE, commentRE, comment2RE, docCommentRE, docComment2RE;

static inline int
arrayArgRank(int argType);

static void
emitWrapper(char *restrict delegateNameBuf, FILE *fpint,
            const char *line,
            enum cftype functype, const char *funcname, size_t funcnameLen,
            size_t funcargc, const int funcargtype[],
            const regmatch_t funcargfull[], const regmatch_t funcargname[],
            size_t maxMatch, regmatch_t *restrict reMatch);

static void
sprintFortranArrayArgDims(size_t argDimsFSize, char argDimsF[argDimsFSize],
                          int argType, const char *argSpecC, size_t maxMatch, regmatch_t reMatch[]);

static void fortran_interface(char *fname, char *fnameinc, char *fnameint, const char *doc_root)
{
  char *line = NULL, *pline;
  size_t lineBufSize = 0;
  char sname[128], *parname;
  char xname[128];
  char *xdes = malloc(128);
  xname[0] = 0;
  size_t xnameLen = 0;
  enum cftype functype;
  int lineno = 0;

  char funcname[MAX_FUNC_NAME_LEN];
  regmatch_t funcargfull[MAX_FUNC_ARGS];
  regmatch_t funcargname[MAX_FUNC_ARGS];
  int  funcargtype[MAX_FUNC_ARGS];
  // char *strsort[99999];
  char timestr[30];
  struct tm *date_and_time;
  regmatch_t *reMatch = NULL;
  size_t maxMatch = 0;
  char **definedParnames = NULL;
  size_t numDefinedParnames = 0;

  time_t date_and_time_in_sec = time(NULL);
  timestr[0] = 0;

  if ( date_and_time_in_sec != -1 )
    {
      date_and_time = localtime(&date_and_time_in_sec);
      (void) strftime(timestr, sizeof(timestr), "%B %Y", date_and_time);
    }

  FILE *fpin = fopen(fname, "r");
  if ( fpin == NULL ) { perror(fname); return; }

  FILE *fpinc = fopen(fnameinc, "w");
  if ( fpinc == NULL ) { perror(fnameinc); return; }

  FILE *fpint = fopen(fnameint, "w");
  if ( fpint == NULL ) { perror(fnameint); return; }

  /* complete symbol table data */
  {
    maxMatch = symRegexCompile(NUM_KNOWN_ARG_TYPES, funArgSym, &line, &lineBufSize);
    size_t maxFunMatch = symRegexCompile(NUM_RET_TYPES, funRet, &line, &lineBufSize);
    if (maxFunMatch > maxMatch)
      maxMatch = maxFunMatch;
  }
  ++maxMatch;
  reMatch = (regmatch_t *)malloc((size_t)maxMatch * sizeof (reMatch[0]));
  /* compile comment start regular expression */
  {
    static const char commentStartREString[] = "^"WS"*/\\*"WS"*(.*"NWS")"WS"*";
    if (reCompile(&commentStartRE, commentStartREString, &line, &lineBufSize))
      exit(EXIT_FAILURE);
  }
  /* compile comment end regular expression */
  {
    static const char commentEndREString[] = "\\*/";
    if (reCompile(&commentEndRE, commentEndREString, &line, &lineBufSize))
      exit(EXIT_FAILURE);
  }
  /* compile complete comment regular expression */
  {
    static const char commentREString[]  = "^"WS"*/\\*"WS"*("NWS"+("WS"+"NWS"+)*)?"WS"*\\*/";
    static const char comment2REString[] = "^"WS"*//"WS"*(.*"NWS")";
    if (reCompile(&commentRE, commentREString, &line, &lineBufSize))
      exit(EXIT_FAILURE);
    if (reCompile(&comment2RE, comment2REString, &line, &lineBufSize))
      exit(EXIT_FAILURE);
  }
  /* compile documentation comment regular expression */
  {
    static const char docCommentREString[]  = "^"WS"*/\\*"WS"*"SYMRE":("WS"*).*"NWS"("WS"*)\\*/";
    static const char docComment2REString[] = "^"WS"*//"WS"*"SYMRE":("WS"*).*"NWS"("WS"*)";
    if (reCompile(&docCommentRE, docCommentREString, &line, &lineBufSize))
      exit(EXIT_FAILURE);
    if (reCompile(&docComment2RE, docComment2REString, &line, &lineBufSize))
      exit(EXIT_FAILURE);
  }
  regex_t cppCondRE;
  {
    static const char cppCondREString[]
      = "^"WS"*#"WS"*((ifn?def)"WS"+\\(?"SYMRE"\\)?|endif)"WS"*(/\\*[^*]*\\*/|//.*)?";
    if (reCompile(&cppCondRE, cppCondREString, &line, &lineBufSize))
      exit(EXIT_FAILURE);
  }
  regex_t cppElseRE;
  {
    static const char cppElseREString[]
      = "^"WS"*#"WS"*else"WS"*(/\\*[^*]*\\*/|//.*)?";
    if (reCompile(&cppElseRE, cppElseREString, &line, &lineBufSize))
      exit(EXIT_FAILURE);
  }
  regex_t emptyStringRE;
  {
    static const char emptyStringREString[] = "^"WS"*";
    if (reCompile(&emptyStringRE, emptyStringREString, &line, &lineBufSize))
      exit(EXIT_FAILURE);
  }

  // fortran include

  fprintf(fpinc, "! This file was automatically generated, don't edit!\n");
  fprintf(fpinc, "!\n");
  fprintf(fpinc, "! Fortran interface for CDI library version %s\n", VERSION);
  fprintf(fpinc, "!\n");
  fprintf(fpinc, "! Author:\n");
  fprintf(fpinc, "! -------\n");
  fprintf(fpinc, "! Uwe Schulzweida, MPI-MET, Hamburg,   %s\n", timestr);
  fprintf(fpinc, "!\n\n");

  // fortran interface

  fprintf(fpint, "// Automatically generated by make_fint.c, don't edit!\n");
  fprintf(fpint, "\n");
  fprintf(fpint, "// clang-format off\n");
  fprintf(fpint, "\n");
  fprintf(fpint, "#ifdef HAVE_CONFIG_H\n");
  fprintf(fpint, "#include \"config.h\"\n");
  fprintf(fpint, "#endif\n");
  fprintf(fpint, "\n");
  char *cppHeaderSentinelMacro;
  size_t cppHeaderSentinelMacroLen;
  {
    char *lastSlash = strrchr(fname, '/');
    char *fbasename = lastSlash ? (lastSlash+1) : fname;
    const size_t fbasenameLen = strlen(fbasename);
    cppHeaderSentinelMacroLen = fbasenameLen + 1;
    cppHeaderSentinelMacro = (char *)malloc(fbasenameLen + 2);
    build_header_name(fbasenameLen, fbasename, cppHeaderSentinelMacro);
    fprintf(fpint, "#ifndef %s\n"
            "#include \"%s\"\n"
            "#endif\n"
            "\n", cppHeaderSentinelMacro, fbasename);
  }
  fputs("#ifdef HAVE_CF_INTERFACE\n"
        "\n"
        "#include <limits.h>\n"
        "#include <assert.h>\n"
        "\n"
        "#ifndef __CFORTRAN_LOADED\n"
        "#  if defined __clang__\n"
        "#    pragma GCC diagnostic push\n"
        "#    pragma GCC diagnostic ignored \"-Wreserved-id-macro\"\n"
        "#  endif\n"
        "#  include \"cfortran.h\"\n"
        "#  if defined __clang__\n"
        "#    pragma GCC diagnostic pop\n"
        "#  endif\n"
        "#endif\n"
        "/* These functions are meant to be called from Fortran and don't\n"
        " * need an interface declaration in a C header. */\n"
        "#ifdef __clang__\n"
        "#  pragma GCC diagnostic push\n"
        "#  pragma GCC diagnostic ignored \"-Wmissing-prototypes\"\n"
        "#endif\n"
        "\n", fpint);


  // TODO: the string below should be generated using build_header_name
  //       for consistency
  fputs("#ifdef CDI_H_\n\n", fpint);
  /*
  fputs("static inline\n"
        "int DateType_c2f(DateType value_DateType)\n"
        "{\n"
        "  assert(value_DateType < INT_MAX);\n"
        "  return (int) value_DateType;\n"
        "}\n"
        "\n", fpint);
  */
  fputs("static inline\n"
        "int SizeType_c2f(SizeType value_SizeType)\n"
        "{\n"
        "  assert(value_SizeType < INT_MAX);\n"
        "  return (int) value_SizeType;\n"
        "}\n"
        "\n", fpint);
  fputs("#endif\n", fpint);

  bool startFortran = false;
  ssize_t lineLen;
  while ((lineLen = getline(&line, &lineBufSize, fpin)) >= 0)
    {
      static const char cplusplus_macro[] = "__cplusplus";
      lineno++;
      if (line[0] == '\n') continue;
      if (strncmp(line, "//FINT_ON", 9) == 0)
        {
          startFortran = true;
          continue;
        }
      if (!startFortran) continue;
      if (strncmp(line, "//FINT_OFF", 10) == 0) break;
      functype = ISVOID;
      pline = line;
      bool needsExtraWrapper = false, needsPrologue = false;
      size_t funcnameLen;
      enum decl declType = UNKNOWN_DECL;
      do {
        for (int retType = 0; retType < NUM_RET_TYPES; ++retType)
          if (!regexec(&funRet[retType].preg, pline, maxMatch, reMatch, 0))
            {
              functype = (enum cftype) retType;
              declType = FUNC_DECL;
              needsPrologue |= funRet[retType].needsPrologue;
              needsExtraWrapper |= funRet[retType].needsExtraWrapper;
              break;
            }
        if (declType == UNKNOWN_DECL)
          break;
        regmatch_t *nameMatch = reMatch + funRet[functype].nameMatch;
        if (debug)
          printf("Found: %.*s\n",
                 (int) (nameMatch->rm_eo - nameMatch->rm_so),
                 pline + nameMatch->rm_so);
        const ssize_t funNameLast = reMatch[0].rm_eo - 1;
        const ssize_t nameLen = nameMatch->rm_eo - nameMatch->rm_so;
        funcnameLen = (size_t)nameLen;
        if ( pline[funNameLast] != '(' )
          {
            fprintf(stderr, "%s\n>(< not found!", line);
            return;
          }
        memcpy(funcname, pline + nameMatch->rm_so, funcnameLen);
        funcname[funcnameLen] = 0;
        pline += reMatch[0].rm_eo;
      } while (0);
      int cppSwitchLen, cppSymLen;

      if (declType == FUNC_DECL)
        {
          size_t funcargc = 0;
	  funcargname[funcargc].rm_so = (regoff_t)(pline - line);
          {
            ssize_t i = 0;
            size_t innerParens = 0;
            do {
              const ssize_t restLen = lineLen - (ssize_t)(pline - line);
              for (; i < restLen; ++i)
                {
                  switch (pline[i])
                    {
                    case ',':
                      if (!innerParens)
                        {
                          funcargc++;
                          funcargname[funcargc].rm_so = (regoff_t)(pline - line + i + 1);
                        }
                      break;
                    case '(':
                      ++innerParens;
                      break;
                    case ')':
                      if (!innerParens)
                        {
                          funcargc++;
                          funcargname[funcargc].rm_so = (regoff_t)(pline - line + i + 1);
                          goto endOfArgSearch;
                        }
                      else
                        --innerParens;
                      break;
                    }
                }
              endOfArgSearch:
              if (i < restLen)
                break;
              char *lineExtension = NULL;
              size_t extSize = 0, plineOff = (size_t)(pline - line);
              ssize_t extLen;
              if ((extLen = getline(&lineExtension, &extSize, fpin)) <= 0)
                break;
              if ((size_t)(lineLen + extLen) >= lineBufSize)
                if (!(line = (char*) realloc(line, (size_t)(lineLen + extLen + 1))))
                  exit(EXIT_FAILURE);
              memcpy(line + lineLen, lineExtension, (size_t)extLen + 1);
              lineLen += extLen;
              pline = line + plineOff;
              ++lineno;
            } while (1);
          }

            /* test if argument list is actually empty */
          if (funcargc == 1
              && !regexec(&emptyStringRE, line + funcargname[0].rm_so, 1, reMatch, 0)
              && (funcargname[0].rm_so + reMatch[0].rm_eo == funcargname[funcargc].rm_so - 1))
            {
              funcargc = 0;
            }
          {
            size_t i;
            for (i = 0; i < funcargc; ++i)
              {
                pline = line + funcargname[i].rm_so;
                int argtype;
                const regoff_t argStart = (regoff_t)(pline - line);
                for (argtype = ISVOID; argtype < NUM_KNOWN_ARG_TYPES; ++argtype)
                  if (!regexec(&funArgSym[argtype].preg, pline, maxMatch, reMatch, 0))
                    {
                      funcargtype[i] = argtype;
                      funcargfull[i].rm_so = reMatch[0].rm_so + argStart;
                      funcargfull[i].rm_eo = reMatch[0].rm_eo + argStart;
                      regmatch_t *nameMatch = reMatch + funArgSym[argtype].nameMatch;
                      funcargname[i].rm_so = nameMatch->rm_so + argStart;
                      funcargname[i].rm_eo = nameMatch->rm_eo + argStart;
                      needsExtraWrapper |= funArgSym[argtype].needsExtraWrapper;
                      needsPrologue |= funArgSym[argtype].needsPrologue;
                      break;
                    }
                if (argtype == NUM_KNOWN_ARG_TYPES)
                  {
                    fprintf(stderr, "%s not implemented\n", line + funcargname[i].rm_so);
                    break;
                  }
              }
            if ( i != funcargc )
              {
                fprintf(stderr, "problem parsing line: %s\n", line);
                continue;
              }
          }

	  strcpy(sname, funcname);

	  // fortran include

	  if ( functype == ISVOID )
	    fprintf(fpinc, "!     %-16s", "");
	  else
	    fprintf(fpinc, "      %-16s", funArgSym[functype].f77name);

          fprintf(fpinc, "%s", sname);
	  fprintf(fpinc, "\n");
	  if ( (funcargc == 1 && funcargtype[0] == ISVOID) ) funcargc = 0;
	  for (size_t i = 0; i < funcargc; i++ )
	    {
              fprintf(fpinc, (i == 0) ? "!%36s(" : ",\n!%36s ", "");
              const int argType = funcargtype[i];
              const int isArray = isArrayArgType(argType);
              enum { argDimsFSize = 128 };
              char argDimsF[argDimsFSize];
              if (!isArray)
                argDimsF[0] = 0;
              else
                sprintFortranArrayArgDims(argDimsFSize, argDimsF, argType,
                                          line + funcargfull[i].rm_so, maxMatch, reMatch);
	      fprintf(fpinc, "%-16s%.*s%s", funArgSym[argType].f77name,
                      (int)(funcargname[i].rm_eo - funcargname[i].rm_so),
                      line + funcargname[i].rm_so, argDimsF);
	    }
	  if ( funcargc )
	    fprintf(fpinc, ")\n");
	  fprintf(fpinc, "      %-16s%s\n\n", "EXTERNAL", sname);

	  /* fortran interface */
          const char *delegateName;
          char delegateNameBuf[MAX_FUNC_NAME_LEN + 7];
          /* emit prologue if needed */
          if (needsPrologue)
            {
              if (funRet[functype].needsPrologue)
                funRet[functype].prologue(fpint, (size_t)-1);
              for (size_t i = 0; i < funcargc; ++i)
                {
                  if (funArgSym[funcargtype[i]].needsPrologue)
                    funArgSym[funcargtype[i]].prologue(fpint, i);
                }
            }
          /* emit wrapper for type conversions if needed */
          if (needsExtraWrapper)
            {
              emitWrapper(delegateNameBuf, fpint, line,
                          functype, funcname, funcnameLen,
                          funcargc, funcargtype, funcargfull, funcargname,
                          maxMatch, reMatch);
              delegateName = delegateNameBuf;
            }
          else
            delegateName = funcname;

	  fputs(functype == ISVOID ? "FCALLSCSUB" : "FCALLSCFUN", fpint);
	  fprintf(fpint, "%zd ", funcargc);
	  fputs("(", fpint);
	  if ( functype != ISVOID )
	    fprintf(fpint, "%s, ", funRet[functype].cfint);
	  fprintf(fpint, "%s, ", delegateName);
	  for (size_t i = 0; i < funcnameLen; ++i)
            sname[i] = (char)toupper((int) sname[i]);
	  fprintf(fpint, "%s, ", sname);
	  for (size_t i = 0; i < funcnameLen; ++i)
            sname[i] = (char)tolower((int) sname[i]);
	  fprintf(fpint, "%s", sname);
	  for (size_t i = 0; i < funcargc; ++i)
	    {
	      fprintf(fpint, ", %s", funArgSym[funcargtype[i]].cfint);
	    }
	  fputs(")\n", fpint);


	  if ( funcnameLen == xnameLen && memcmp(funcname, xname, funcnameLen) == 0 )
	    {
	      char xline[256];
              size_t xlineLen = 0;
	      int nch;

	      /* C Quick Guide */

	      cdoc[ncdoc].naline = 0;
	      cdoc[ncdoc].text   = NULL;
	      cdoc[ncdoc].fname  = strdup(funcname);

	      nch = sprintf(xline, funRet[functype].cfmt, "", (int)funcnameLen, funcname);
              xline[nch  ] = ' ';
              xline[nch+1] = '(';
              xline[nch+2] = '\0';
              nch += 2;
              xlineLen = (size_t)nch;

	      if ( (funcargc == 1 && funcargtype[0] == ISVOID) ) funcargc = 0;

	      for (size_t i = 0; i < funcargc; ++i)
		{
		  if (i)
                    {
                      xline[xlineLen  ] = ',';
                      xline[xlineLen+1] = ' ';
                      xline[xlineLen+2] = '\0';
                      xlineLen += 2;
                    }

                  /* extract full argument text from match */
                  char farg[128];
                  /* - 1 to omit closing paren ) or comma , */
                  int nchn = snprintf(farg, sizeof (farg), "%.*s",
                                      (int)(funcargfull[i].rm_eo - funcargfull[i].rm_so - 1),
                                      line + funcargfull[i].rm_so);
                  if (nchn < 0) abort();
                  /* compress white-space */
                  nchn = (int)compress_whitespace((size_t)nchn, farg);
		  if ( (xlineLen + (size_t)nchn) > (size_t)80 )
		    {
                      if (i) xline[--xlineLen] = 0;
		      cdoc[ncdoc].aline[cdoc[ncdoc].naline++] = strdup(xline);
		      sprintf(xline, "%*s", nch, "");
                      xlineLen = (size_t)nch;
		    }
		  strcat(xline, farg);
                  xlineLen += (size_t)nchn;
		}
              xline[xlineLen  ] = ')';
              xline[xlineLen+1] = ';';
              xline[xlineLen+2] = '\0';
              xlineLen += 2;
	      cdoc[ncdoc].aline[cdoc[ncdoc].naline++] = strdup(xline);
	      cdoc[ncdoc].text  = strdup(xdes);

	      ncdoc++;

	      /* Fortran Quick Guide */

	      fdoc[nfdoc].naline = 0;
	      fdoc[nfdoc].text   = NULL;
	      fdoc[nfdoc].fname  = strdup(funcname);

              nch = sprintf(xline, "%s%s %s", funArgSym[functype].f77name,
                            functype == ISVOID ? "SUBROUTINE" : " FUNCTION",
                            xname);

	      if ( (funcargc == 1 && funcargtype[0] == ISVOID) ) funcargc = 0;
              if (funcargc)
                {
                  xline[nch  ] = ' ';
                  xline[nch+1] = '(';
                  xline[nch+2] = '\0';
                  nch += 2;
                }

              xlineLen = (size_t)nch;

	      for (size_t i = 0; i < funcargc; ++i)
		{
		  if (i)
                    {
                      xline[xlineLen  ] = ',';
                      xline[xlineLen+1] = ' ';
                      xline[xlineLen+2] = '\0';
                      xlineLen += 2U;
                    }

                  enum { argDimsFSize = 128 };
                  char farg[128], argDimsF[argDimsFSize];
                  /* FIXME: optional empty argument name unhandled */
                  const int argType = funcargtype[i];
                  const int isArray = isArrayArgType(argType);
                  if (!isArray)
                    argDimsF[0] = 0;
                  else
                    sprintFortranArrayArgDims(argDimsFSize, argDimsF, argType,
                                              line + funcargfull[i].rm_so, maxMatch, reMatch);
		  const int nchn
                    = snprintf(farg, sizeof (farg), "%s %.*s%s", funArgSym[argType].f77name,
                               (int)(funcargname[i].rm_eo - funcargname[i].rm_so),
                               line + funcargname[i].rm_so, argDimsF);
                  if (nchn < 0) abort();
		  if ( (xlineLen + (size_t)nchn) > 80 )
		    {
                      if (i) xline[--xlineLen] = 0;
		      fdoc[nfdoc].aline[fdoc[nfdoc].naline++] = strdup(xline);
		      sprintf(xline, "%*s", nch, "");
                      xlineLen = (size_t)nch;
		    }
		  strcat(xline, farg);
                  xlineLen += (size_t)nchn;
		}
	      if ( funcargc )
                {
                  xline[xlineLen  ] = ')';
                  xline[xlineLen+1] = '\0';
                  xlineLen += 1;
                }
	      fdoc[nfdoc].aline[fdoc[nfdoc].naline++] = strdup(xline);
	      fdoc[nfdoc].text  = strdup(xdes);

	      nfdoc++;
	    }
	}
      else if ( memcmp(line, "#define", 7) == 0 )
	{
	  pline = line + 7;
	  while ( isspace((int) *pline) ) pline++;
	  parname = pline;
	  size_t len = strlen(pline);
          size_t parnameLen = 0;
	  while (parnameLen < len && !isspace((int)pline[parnameLen]))
            ++parnameLen;
	  if ( parnameLen == len ) continue;
	  pline += parnameLen+1;
	  while ( isspace((int) *pline) ) pline++;
	  if ( isdigit((int)*(pline + (*pline == '-'))) )
            {
              parname[parnameLen] = 0;
              const int parvalue = atoi(pline);
              /* fortran include */
              fprintf(fpinc, "      INTEGER    %s\n"
                      "      PARAMETER (%-22s = %2d)\n", parname, parname, parvalue);
            }
          else if (isSymStart((int)*pline) || (*pline == '-' && isSymStart(pline[1])))
            {
              parname[parnameLen] = 0;
              size_t parValueLen = 1 + (*pline == '-');
              while (isSym((int)pline[parValueLen]))
                ++parValueLen;
              char *parValue = pline;
              for (size_t i = 0; i < numDefinedParnames; ++i)
                if (!strcmp(definedParnames[i], parValue + (*pline == '-')))
                  goto foundParname;
              fprintf(stderr, "Found definition for %s to unknown value: %*s\n", parname, (int)parValueLen, parValue);
              foundParname:
              fprintf(fpinc, "      INTEGER    %s\n"
                      "      PARAMETER (%-22s = %1.*s)\n",
                      parname, parname, (int)parValueLen, parValue);
            }
          else
            {
              if ( strncmp(parname, cppHeaderSentinelMacro,
                           cppHeaderSentinelMacroLen) == 0 ) continue;
              fprintf(fpinc, "%s", line);
              continue;
            }
          definedParnames = realloc(definedParnames, sizeof (definedParnames[0]) * (numDefinedParnames+1));
          if (!definedParnames)
            abort();
          char *remembered = definedParnames[numDefinedParnames] = malloc(parnameLen+1);
          if (!remembered)
            abort();
          memcpy(remembered, parname, parnameLen+1);
          ++numDefinedParnames;
	}
      else if (!regexec(&cppCondRE, line, maxMatch, reMatch, 0)
               && ((cppSwitchLen = reMatch[2].rm_eo - reMatch[2].rm_so) == 5)
               && ((size_t)(cppSymLen = reMatch[3].rm_eo - reMatch[3].rm_so) == sizeof (cplusplus_macro) - 1)
               && !memcmp(line + reMatch[3].rm_so, cplusplus_macro, sizeof (cplusplus_macro) - 1))
	{
          fprintf(stderr, "Found conditional C++ block, skipping to #else\n");
          while ((lineLen = getline(&line, &lineBufSize, fpin)) >= 0)
            {
              ++lineno;
              static const char endif_str[] = "endif";
              if (!regexec(&cppElseRE, line, maxMatch, reMatch, 0)
                  || (!regexec(&cppCondRE, line, maxMatch, reMatch, 0)
                      && ((size_t)(cppSymLen = reMatch[1].rm_eo - reMatch[1].rm_so) == sizeof (endif_str) - 1)
                      && !memcmp(line + reMatch[1].rm_so, endif_str, sizeof (endif_str) - 1)))
                break;
            }
        }
      else if (detectComment(&line, &lineLen, &lineBufSize, maxMatch, reMatch,
                             xname, &xnameLen, &xdes, &lineno, fpin, fpinc, fpint))
        ;
      else
	{
	  if ( lineLen > 1 )
	    fprintf(stderr, "skip: line %3d  size %3zd  %s%s", lineno, lineLen, line,
                   line[lineLen-1]=='\n'?"":"missing new-line\n");
	}
    }

  fputs("\n"
        "#if defined __clang__\n"
        "#  pragma GCC diagnostic pop\n"
        "#endif\n"
        "\n"
        "// clang-format on\n"
        "\n"
        "#endif\n", fpint);

  fclose(fpin);
  fclose(fpinc);
  fclose(fpint);

  qsort(cdoc, ncdoc, sizeof(Docu), doccmp);
  qsort(fdoc, nfdoc, sizeof(Docu), doccmp);

  if ( NULL == doc_root ) return;

  char *filename;
  const size_t doc_root_len = strlen(doc_root);
  {
    size_t max_len = 0;
    for (size_t i = 0; i < (size_t)fname_list_size; ++i)
      {
        size_t len = strlen(fname_list[i].name);
        fname_list[i].len = len;
        if (len > max_len)
          max_len = len;
      }
    /* path to document root, separating /, max of path within root, terminating '\0'  */
    filename = (char*) malloc(doc_root_len + 1 + max_len + 1);
  }

  memcpy(filename, doc_root, doc_root_len);
  filename[doc_root_len] = '/';
  memcpy(filename + doc_root_len + 1,
         fname_list[NAME_C_QUICK_REF].name,
         fname_list[NAME_C_QUICK_REF].len + 1);
  FILE *fp = fopen(filename, "w");
  if ( fp )
    {
      fputs("C Quick Reference\n"
            "-----------------\n\n", fp);

      doctotxt(fp, cdoc, ncdoc);
      fclose(fp);
    }
  else
    {
      fprintf(stderr, "warning: cannot open documentation output file %s, %s",
              filename, strerror(errno));
    }

  memcpy(filename + doc_root_len + 1,
         fname_list[NAME_F_QUICK_REF].name,
         fname_list[NAME_F_QUICK_REF].len + 1);
  fp = fopen(filename, "w");
  if ( fp )
    {
      fputs("Fortran Quick Reference\n"
            "-----------------------\n\n", fp);

      doctotxt(fp, fdoc, nfdoc);
      fclose(fp);
    }
  else
    {
      fprintf(stderr, "warning: cannot open documentation output file %s, %s", filename, strerror(errno));
    }

  memcpy(filename + doc_root_len + 1,
         fname_list[NAME_C_QUICK_REF_TEX].name,
         fname_list[NAME_C_QUICK_REF_TEX].len + 1);
  fp = fopen(filename, "w");
  if ( fp )
    {
      fputs("\\chapter*{Quick Reference}\n"
            "\\addcontentsline{toc}{chapter}{Quick Reference}\n"
            "\n"
            "This appendix provide a brief listing of the C language bindings of the CDI library routines:\n"
            "\n", fp);

      doctotex(fp, cdoc, ncdoc);
      fclose(fp);
    }
  else
    {
      fprintf(stderr, "warning: cannot open documentation output file %s, %s", filename, strerror(errno));
    }

  memcpy(filename + doc_root_len + 1,
         fname_list[NAME_F_QUICK_REF_TEX].name,
         fname_list[NAME_F_QUICK_REF_TEX].len + 1);
  fp = fopen(filename, "w");
  if ( fp )
    {
      fputs("\\chapter*{Quick Reference}\n"
            "\\addcontentsline{toc}{chapter}{Quick Reference}\n"
            "\n"
            "This appendix provide a brief listing of the Fortran language bindings of the CDI library routines:\n"
            "\n", fp);

      doctotex(fp, fdoc, nfdoc);
      fclose(fp);
    }
  else
    {
      fprintf(stderr, "warning: cannot open documentation output file %s, %s", filename, strerror(errno));
    }
  free(filename);
}

static void
build_header_name(const size_t len, const char *fname, char *cppMacro)
{
  for (size_t i = 0; i < len; ++i)
    switch (fname[i])
      {
      case '.':
        cppMacro[i] = '_';
        break;
      default:
        cppMacro[i] = (char)toupper((int)fname[i]);
      }
  cppMacro[len] = '_';
  cppMacro[len + 1] = '\0';
}

int main(int argc, char *argv[])
{
  const char *doc_root = NULL, *out_root = NULL;

  int optargCount = 0;
  {
    int opt;
    while ((opt = getopt(argc, argv, "d:o:")) != -1)
      switch (opt) {
      case 'd':
        doc_root = optarg;
        optargCount += 2;
        break;
      case 'o':
        out_root = optarg;
        optargCount += 2;
        break;
      default: /* '?' */
        fprintf(stderr, "Usage: %s [-d DOCROOT] [-o OUTROOT] includefile\n", argv[0]);
        return EXIT_FAILURE;
      }
  }


  if ( argc != 2 + optargCount)
    {
      printf("Usage:  %s [-d DOCROOT] [-o OUTROOT] includefile\n", argv[0]);
      return 1;
    }

  char *fname = argv[1 + optargCount];

  const char *cs = strrchr(fname, '/');
  const char *fbasename = (cs == NULL) ? fname : cs + 1;

  const char *cp = strrchr(fbasename, '.');
  const size_t fbasename_len = (cp == NULL) ? strlen(fbasename) : (size_t) (cp - fbasename);

  const char *fprefix;
  size_t fprefix_len;

  if (out_root == NULL)
    {
      fprefix = fname;
      fprefix_len = (size_t) (fbasename - fname);
    }
  else
    {
      fprefix = out_root;
      fprefix_len = strlen(out_root);
    }

  char fnameinc[128], fnameint[128];

  if (fprefix_len > 0)
    {
      memcpy(fnameinc, fprefix, fprefix_len);
      memcpy(fnameint, fprefix, fprefix_len);

      if (fprefix[fprefix_len - 1] != '/')
        {
          fnameinc[fprefix_len] = fnameint[fprefix_len] = '/';
          fprefix_len += 1;
        }
    }

  size_t fname_idx = fprefix_len;
  memcpy(fnameinc + fname_idx, fbasename, fbasename_len);
  memcpy(fnameint + fname_idx, fbasename, fbasename_len);

  fname_idx += fbasename_len;
  strcpy(fnameinc + fname_idx, ".inc");
  strcpy(fnameint + fname_idx, "Fortran.c");

  fortran_interface(fname, fnameinc, fnameint, doc_root);

  return 0;
}

static inline size_t
compress_whitespace(size_t len, char str[])
{
  size_t wpos = 0;
  size_t i = 0;
  /* skip leading white-space */
  while (i < len && (isblank(str[i]) || str[i] == '\n'))
    ++i;
  /* after the leading white-space the following is
   * an alternation of white- and non-white-space
   * characters, where sequences of the latter will
   * be compressed to a single space */
  while (i < len)
    {
      /* handle white-space */
      while (i < len && !(isblank(str[i]) || str[i] == '\n'))
        str[wpos++] = str[i++];
      /* skip non-white-space */
      size_t wscount = 0;
      while (i < len && (isblank(str[i]) || str[i] == '\n'))
        ++i, ++wscount;
      if (wscount)
        str[wpos++] = ' ';
    }
  str[wpos] = '\0';
  return wpos;
}

enum {
  REGEX_MAX_ERRSTRLEN = 1024,
};


static size_t
symRegexCompile(size_t numSyms, struct symbol symList[], char **line, size_t *lineBufSize)
{
  size_t maxMatch = 0;
  for (size_t sym = 0; sym < numSyms; ++sym)
    {
      if (reCompile(&symList[sym].preg, symList[sym].parseRE, line, lineBufSize))
        exit(EXIT_FAILURE);
      const size_t numMatches = symList[sym].nameMatch + (size_t)arrayArgRank((int)sym);
      if (numMatches > maxMatch)
        maxMatch = numMatches;
    }
  return maxMatch;
}

static int
reCompile(regex_t *restrict RE, const char *restrict REstring,
          char * restrict *restrict lineBuf, size_t * restrict lineBufSize)
{
  int errcode;
  if ((errcode = regcomp(RE, REstring, REG_EXTENDED)))
    {
      char *restrict line;
      size_t resize;
      if (*lineBufSize < REGEX_MAX_ERRSTRLEN
          && (line = (char*) realloc(*lineBuf, resize = REGEX_MAX_ERRSTRLEN)))
        {
          *lineBuf = line;
          *lineBufSize = resize;
          regerror(errcode, RE, line, *lineBufSize);
          fprintf(stderr, "Error compiling regular expression: %s: %s\n", REstring, *lineBuf);
        }
    }
  return errcode;
}

// emit conversion code for DateType argument
static int cfDateTypeConvert(FILE *outfp, const char *argName,
                             size_t argNameLen, enum conversionType part)
{
  int retval = 0;
  switch (part)
    {
    case CONV_ARG:
      retval = fprintf(outfp, "(DateType)%.*s", (int)argNameLen, argName);
      break;
    case CONV_RET:
      retval = fprintf(outfp, "DateType_c2f(%.*s)", (int)argNameLen, argName);
      break;
    }
  return retval;
}

// emit conversion code for SizeType argument
static int cfSizeTypeConvert(FILE *outfp, const char *argName,
                             size_t argNameLen, enum conversionType part)
{
  int retval = 0;
  switch (part)
    {
    case CONV_ARG:
      retval = fprintf(outfp, "(SizeType)%.*s", (int)argNameLen, argName);
      break;
    case CONV_RET:
      retval = fprintf(outfp, "SizeType_c2f(%.*s)", (int)argNameLen, argName);
      break;
    }
  return retval;
}

// emit conversion code for DateType* argument
static int cfDateTypepConvert(FILE *outfp, const char *argName,
                              size_t argNameLen, enum conversionType part)
{
  int retval = 0;
  switch (part)
    {
    case CONV_ARG:
      retval = fprintf(outfp, "&%.*s_DateType", (int)argNameLen, argName);
      break;
    case CONV_RET:
      abort();

      break;
    }
  return retval;
}

// emit conversion code for SizeType* argument
static int cfSizeTypepConvert(FILE *outfp, const char *argName,
                              size_t argNameLen, enum conversionType part)
{
  int retval = 0;
  switch (part)
    {
    case CONV_ARG:
      retval = fprintf(outfp, "&%.*s_SizeType", (int)argNameLen, argName);
      break;
    case CONV_RET:
      abort();

      break;
    }
  return retval;
}

/* emit conversion code for MPI_Comm argument */
static int cfMPICommConvert(FILE *outfp, const char *argName,
                            size_t argNameLen, enum conversionType part)
{
  int retval = 0;
  switch (part)
    {
    case CONV_ARG:
      retval = fprintf(outfp, "MPI_Comm_f2c(%.*s)", (int)argNameLen, argName);
      break;
    case CONV_RET:
      retval = fprintf(outfp, "MPI_Comm_c2f(%.*s)", (int)argNameLen, argName);
      break;
    }
  return retval;
}

/* emit conversion code for Xt_idxlist argument */
static int cfXtIdxlistConvert(FILE *outfp, const char *argName,
                              size_t argNameLen, enum conversionType part)
{
  int retval = 0;
  switch (part)
    {
    case CONV_ARG:
      retval = fprintf(outfp, "(*(Xt_idxlist *)%.*s)", (int)argNameLen, argName);
      break;
    case CONV_RET:
      abort();
      break;
    }
  return retval;
}

static int cfVoidFuncPrologue(FILE *outfp, size_t argNum)
{
  const int retval
    = fprintf(outfp, "\n#undef ROUTINE_%zu\n#define ROUTINE_%zu %s\n",
              argNum+1, argNum+1, "(void (*)(void))");
  return retval;
}

enum {
  FOUND_NOTHING,
  FOUND_COMMENT,
  FOUND_DOCCOMMENT,
};

static int detectComment(char **line_, ssize_t *lineLen, size_t *lineBufSize,
                         size_t maxMatch, regmatch_t reMatch[],
                         char *xname, size_t *xnameLen, char **xdes_,
                         int *lineno, FILE *fpin, FILE *fpinc, FILE *fpint)
{
  char *restrict line = *line_;
  int matchType;
  do {
    if (!regexec(&docCommentRE, line, maxMatch, reMatch, 0) || !regexec(&docComment2RE, line, maxMatch, reMatch, 0))
      {
        /* found documentation comment */
        size_t nameMatchLen = (size_t)(reMatch[1].rm_eo - reMatch[1].rm_so),
          docMatchLen = (size_t)(reMatch[3].rm_so - reMatch[2].rm_eo);
        memcpy(xname, line + reMatch[1].rm_so, nameMatchLen);
        xname[nameMatchLen] = 0;
        *xnameLen = nameMatchLen;
        char *xdes = *xdes_ = realloc(*xdes_, docMatchLen + 1);
        memcpy(xdes, line + reMatch[2].rm_eo, docMatchLen);
        xdes[docMatchLen] = 0;
        {
          char *eol = xdes;
          while ((eol = strchr(eol, '\n')))
            {
              ++eol;
              /* delete whitespace following newline */
              const size_t squeezeLen = strspn(eol, " \t*");
              char *startoftext = eol + squeezeLen;
              memmove(eol, startoftext, docMatchLen - (size_t)(startoftext - xdes));
              docMatchLen -= squeezeLen;
              xdes[docMatchLen] = 0;
            }
        }
        if (verbose || debug)
          printf("Found documentation for \"%s\": \"%s\"\n", xname, xdes);
        matchType = FOUND_DOCCOMMENT;
        break;
      }
    else if (!regexec(&commentRE, line, maxMatch, reMatch, 0) || !regexec(&comment2RE, line, maxMatch, reMatch, 0))
      {
        const size_t commentLen = (size_t)(reMatch[1].rm_eo - reMatch[1].rm_so);
        char *comment = line + reMatch[1].rm_so;
        {
          const char savedCommentEnd = comment[commentLen];
          comment[commentLen] = '\0';
          /* fortran include */
          fputs("!\n", fpinc);
          char *cline = comment;
          do {
            cline += strspn(cline, " \t*");
            char *eol = strchr(cline, '\n');
            if (!eol) eol = comment + commentLen;
            size_t lineLen = (size_t)(eol - cline);
            fprintf(fpinc, "!  %.*s\n", (int)lineLen, cline);
            cline = (eol != comment + commentLen) ? eol + 1: NULL;
          } while (cline);
          fputs("!\n", fpinc);
          comment[commentLen] = savedCommentEnd;
        }
        /* fortran interface */
        fprintf(fpint, "\n/*  %.*s  */\n\n", (int)commentLen, comment);
        matchType = FOUND_COMMENT;
        break;
      }
    /* found comment start, read further lines and retry */
    else if (!regexec(&commentStartRE, line, maxMatch, reMatch, 0))
      {
        int foundCommentEnd = 0;
        char *lineExtension = NULL;
        size_t extSize = 0;
        do {
          ssize_t extLen;
          if ((extLen = getline(&lineExtension, &extSize, fpin)) <= 0)
            break;
          if ((size_t)(*lineLen + extLen) >= *lineBufSize)
            if (!(line = realloc(line, (size_t)(*lineLen + extLen + 1))))
              exit(EXIT_FAILURE);
          memcpy(line + *lineLen, lineExtension, (size_t)extLen + 1);
          *lineLen += extLen;
          ++(*lineno);
          foundCommentEnd
            = !regexec(&commentEndRE, lineExtension, maxMatch, reMatch, 0);
        } while (!foundCommentEnd);
      }
    else
      {
        /* found no comment at all */
        matchType = 0;
        break;
      }
  } while (1);
  *line_ = line;
  return matchType;
}

static inline int
isArrayArgType(const int argType)
{
  return argType == ISFLOATV
    || argType == ISFLOATVV
    || argType == ISDOUBLEV
    || argType == ISDOUBLEVV
    || argType == ISXT_IDXLISTV
    || argType == ISINTV
    || argType == ISINTVV
    || argType == ISINTVVV;
}

static inline int
arrayArgRank(const int argType)
{
  int rank = 0;
  switch (argType) {
  case ISFLOATV:
  case ISDOUBLEV:
  case ISINTV:
  case ISXT_IDXLISTV:
    rank = 1;
    break;
  case ISINTVV:
  case ISFLOATVV:
  case ISDOUBLEVV:
    rank = 2;
    break;
  case ISINTVVV:
    rank = 3;
    break;
  }
  return rank;
}

static void
sprintCArrayArgDims(size_t argDimsCSize, char argDimsC[argDimsCSize],
                    int argType, const char *argSpecC,
                    const regmatch_t reMatch[]);

static void
sprintFortranArrayArgDims(size_t argDimsFSize, char argDimsF[argDimsFSize],
                          int argType, const char *argSpecC,
                          size_t maxMatch, regmatch_t reMatch[])
{
  if (regexec(&funArgSym[argType].preg, argSpecC, maxMatch,
              reMatch, 0))
    {
      fprintf(stderr, "unexpected non-matching of array argument regexp!\n");
      exit(1);
    }
  argDimsF[0] = '(';
  size_t argDimsFPos = 1;
  const int arank = arrayArgRank(argType);
  for (int rank = arank; rank > 0; --rank)
    {
      const size_t rankBoundMatch
        = funArgSym[argType].nameMatch + (size_t)rank;
      const size_t boundStringLen
        = (size_t)(reMatch[rankBoundMatch].rm_eo
                   - reMatch[rankBoundMatch].rm_so);
      if (!boundStringLen)
        argDimsF[argDimsFPos++] = '*';
      else
        {
          memcpy(argDimsF + argDimsFPos,
                 argSpecC + reMatch[rankBoundMatch].rm_so,
                 boundStringLen);
          argDimsFPos += boundStringLen;
        }
      argDimsF[argDimsFPos++] = rank > 1 ? ',' : ')';
    }
  argDimsF[argDimsFPos++] = 0;
}

/*
 * Creates a wrapper function to call funcname indirectly.
 *
 * This is needed in case a function needs a wrapper going from
 * Fortran to C or back because a variable type cannot be used
 * directly like e.g. MPI_Comm
 */
static void
emitWrapper(char *restrict delegateNameBuf, FILE *fpint,
            const char *line,
            enum cftype functype, const char *funcname, size_t funcnameLen,
            size_t funcargc, const int funcargtype[],
            const regmatch_t funcargfull[], const regmatch_t funcargname[],
            size_t maxMatch, regmatch_t *restrict reMatch)
{
  static const char fwrap_suffix[] = "_fwrap";
  size_t delegateNameLen = funcnameLen;
  strcpy(delegateNameBuf, funcname);
  strcpy(delegateNameBuf+delegateNameLen, fwrap_suffix);
  delegateNameLen += sizeof (fwrap_suffix) - 1;
  const char *delegateName = delegateNameBuf;
  fputs("static ", fpint);
  fprintf(fpint, (funRet[functype].convert ? funRet[functype].convcfmt : funRet[functype].cfmt),
          "", (int)delegateNameLen, delegateName);
  fputs("(", fpint);
  for (size_t i = 0; i < funcargc; ++i)
    {
      if (i > 0) fputs(", ", fpint);
      enum { arrayDimsCSize = 128 };
      char arrayDimsC[arrayDimsCSize];
      static const char constStr[] = "const ", nonConstStr[] = "";
      const char *constStrP;
      const int argType = funcargtype[i];
      const int isArray = isArrayArgType(argType);
      if (!isArray)
        {
          arrayDimsC[0] = 0;
          constStrP = nonConstStr;
        }
      else
        {
          if (regexec(&funArgSym[argType].preg, line + funcargfull[i].rm_so, maxMatch, reMatch, 0))
            {
              fprintf(stderr, "unexpected non-matching of array argument regexp!\n");
              exit(1);
            }
          if (reMatch[1].rm_eo - reMatch[1].rm_so > 5
              && !strncmp(line + funcargfull[i].rm_so + reMatch[1].rm_so, constStr, sizeof (constStr) - 2))
            constStrP = constStr;
          else
            constStrP = nonConstStr;
          sprintCArrayArgDims(arrayDimsCSize, arrayDimsC, argType, line + funcargfull[i].rm_so, reMatch);
        }
      fprintf(fpint, (funArgSym[argType].convert ? funArgSym[argType].convcfmt : funArgSym[argType].cfmt),
              constStrP, (int)(funcargname[i].rm_eo - funcargname[i].rm_so),
              line + funcargname[i].rm_so, arrayDimsC);
    }
  fputs(")\n{\n", fpint);
  for (size_t i = 0; i < funcargc; ++i)
    {
      if (funArgSym[funcargtype[i]].convert && funcargtype[i] == ISDATETYPEP)
        {
          // create temporary DateType variable
          fprintf(fpint, "  DateType %.*s_DateType;\n",
                  (int)(funcargname[i].rm_eo - funcargname[i].rm_so), line + funcargname[i].rm_so);
        }
      if (funArgSym[funcargtype[i]].convert && funcargtype[i] == ISSIZETYPEP)
        {
          // create temporary SizeType variable
          fprintf(fpint, "  SizeType %.*s_SizeType;\n",
                  (int)(funcargname[i].rm_eo - funcargname[i].rm_so), line + funcargname[i].rm_so);
        }
    }
  if (functype != ISVOID)
    {
      fputs("  ", fpint);
      fprintf(fpint, funRet[functype].cfmt, "", 1, "v");
      fprintf(fpint, ";\n");
      fprintf(fpint, "  v = %s(", funcname);
    }
  else
    fprintf(fpint, "  %s(", funcname);
  for (size_t i = 0; i < funcargc; ++i)
    {
      if (i > 0) fputs(", ", fpint);
      if (funArgSym[funcargtype[i]].convert)
        {
          funArgSym[funcargtype[i]].convert(fpint, line + funcargname[i].rm_so,
                                            (size_t)(funcargname[i].rm_eo - funcargname[i].rm_so), CONV_ARG);
        }
      else
        fprintf(fpint, "%.*s", (int)(funcargname[i].rm_eo - funcargname[i].rm_so), line + funcargname[i].rm_so);
    }
  fputs(");\n", fpint);
  for (size_t i = 0; i < funcargc; ++i)
    {
      if (funArgSym[funcargtype[i]].convert && funcargtype[i] == ISDATETYPEP)
        {
          fprintf(fpint, "  assert(%.*s_DateType < INT_MAX);\n",
                  (int)(funcargname[i].rm_eo - funcargname[i].rm_so), line + funcargname[i].rm_so);
          // copy temporary DateType variable
          fprintf(fpint, "  *%.*s = %.*s_DateType;\n",
                  (int)(funcargname[i].rm_eo - funcargname[i].rm_so), line + funcargname[i].rm_so,
                  (int)(funcargname[i].rm_eo - funcargname[i].rm_so), line + funcargname[i].rm_so);
        }
      if (funArgSym[funcargtype[i]].convert && funcargtype[i] == ISSIZETYPEP)
        {
          fprintf(fpint, "  assert(%.*s_SizeType < INT_MAX);\n",
                  (int)(funcargname[i].rm_eo - funcargname[i].rm_so), line + funcargname[i].rm_so);
          // copy temporary SizeType variable
          fprintf(fpint, "  *%.*s = %.*s_SizeType;\n",
                  (int)(funcargname[i].rm_eo - funcargname[i].rm_so), line + funcargname[i].rm_so,
                  (int)(funcargname[i].rm_eo - funcargname[i].rm_so), line + funcargname[i].rm_so);
        }
    }
  if (functype != ISVOID)
    {
      fputs("  return ", fpint);
      if (funRet[functype].convert)
        funRet[functype].convert(fpint, "v", 1, CONV_RET);
      else
        fputc('v', fpint);
      fputs(";\n", fpint);
    }
  fputs("}\n", fpint);
}

static void
sprintCArrayArgDims(size_t argDimsCSize, char argDimsC[argDimsCSize],
                    int argType, const char *argSpecC,
                    const regmatch_t reMatch[])
{
  size_t argDimsCPos = 0;
  const int arank = arrayArgRank(argType);
  for (int rank = 1; rank <= arank; ++rank)
    {
      const size_t rankBoundMatch
        = funArgSym[argType].nameMatch + (size_t)rank;
      const size_t boundStringLen
        = (size_t)(reMatch[rankBoundMatch].rm_eo
                   - reMatch[rankBoundMatch].rm_so);
      argDimsC[argDimsCPos++] = '[';
      if (boundStringLen)
        {
          memcpy(argDimsC + argDimsCPos,
                 argSpecC + reMatch[rankBoundMatch].rm_so,
                 boundStringLen);
          argDimsCPos += boundStringLen;
        }
      argDimsC[argDimsCPos++] = ']';
    }
  argDimsC[argDimsCPos++] = 0;
}

// clang-format on

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */

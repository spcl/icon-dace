#!/usr/bin/python

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

# -*- coding: utf-8 -*-
#-----------------------------------------------------------------------------
#
# check_style.py is a script for reformating Fortran code
# according to basic style rules for ICON
#
# Revision History:
#    L.Linardakis, MPI-M 2010-01-14
#
# Usage: check_style inputName, outputName, outScriptName
#   inputName = source code filename
#   outputName = the cheked and possibly reformed file name
#   outScriptName = output script for possible further changes using sed
#
# The script for the moment performs very basic checks:
#   1. All fortran keywords turned into upper case
#   2. All other names are turned to lower case
#   3. Checks for name with length > maxNameLength = 32.
#      Gives a warning and outputs to outScriptName
#   4. Changes tabs to spaces = noOfTabSpaces (=2)
#   5. reformats to the correct indentation
# 
# NOTE:
#   1. Comments and Strings are UNCHANGED
#      This can potentially be problematic and should be checked
#   2. Line continuations (&) is not taken into account for processing a fortran line
#      They will be maintained, and will be reproduced in the new code
#
# Bug Report:
#
#
#-----------------------------------------------------------------------------

import getopt, sys

noOfTabSpaces=2   # replace tabs with that many spaces
maxNameLength=31  # max variabble name length
maxLineLength=99  # max line length

compile_comment_directives=[
#these should immediately follow the !
#add compile directives as necessary
'$OMP ',
'$POMP ',
'CDIR ',
'DIR$ ',
'OCL '
]


fortran_keywords=[
"INTEGER", # DATA TYPES
"REAL",
"COMPLEX",
"LOGICAL",
"CHARACTER",
"ACCESS",   #FORTRAN 77
"ASSIGN",
"BACKSPACE",
"BLANK",
"BLOCK",
"CALL",
"CLOSE",
"COMMON",
"CONTINUE",
"DATA",
"DIMENSION",
"DIRECT",
"DO",
"ELSE",
"ELSEIF",
"ENDIF",
"ENDDO",
"END",
"ENDFUNCTION",
"ENDWHERE",
"ENTRY",
"EOF",
"EQUIVALENCE",
"ERR",
"EXIST",
"EXTERNAL",
"FILE",
"FMT",
"FORM",
"FORMAT",
"FORMATTED",
"FUNCTION",
"GO",
"TO",
"IF",
"IMPLICIT",
"INCLUDE",
"INQUIRE",
"INTRINSIC",
"IOSTAT",
"LOGICAL",
"NAMED",
"NAMELIST",
"NEXTREC",
"NUMBER",
"OPEN",
"OPENED",
"PARAMETER",
"PAUSE",
"PRINT",
"PROGRAM",
"READ",
"REC",
"RECL",
"RETURN",
"REWIND",
"SEQUENTIAL",
"STATUS",
"STOP",
"SUBROUTINE",
"THEN",
"TYPE",
"UNFORMATTED",
"UNIT",
"WRITE",
"SAVE",
# Logical 
"TRUE",
"FALSE",
"AND",
"OR",
"NOT",
"GT",
"LT",
"GE",
"LE",
"EQ",
# FORTRAN 90
"ALLOCATE",  
"ALLOCATABLE",
"ALLOCATED",
"CASE",
"CONTAINS",
"CYCLE",
"DEALLOCATE",
"ELSEWHERE",
"EXIT",
"INTERFACE",
"INTENT",
"MODULE",
"ONLY",
"OPERATOR",
"OPTIONAL",
"POINTER",
"PRIVATE",
"PROCEDURE",
"PUBLIC",
"RECURSIVE",
"SELECT",
"SEQUENCE",
"TARGET",
"TRIM",
"USE",
"WHILE",
"WHERE",
"NONE",
"ELEMENTAL", # FORTRAN 95
"FORALL",
"PURE",
"ABSTRACT",  # FORTRAN 2003
"ASSOCIATE",
"CLASS",
"DECIMAL",
"DECORATE",
"DELEGATE",
"ENCODING",
"ENDFILE",
"ENUM",
"ENUMERATOR",
"EXTENDS",
"EXTENSIBLE",
"FLUSH",
"GENERIC",
"IOMSG",
"IMPORT",
"MOVE_ALLOC",
"NEXTREC",
"NON_OVERRIDABLE",
"PASS",
"PENDING",
"REFERENCE",
"ROUND",
"STATIC",
"TYPEALIAS",
"NULLIFY",
"ASYNCHRONOUS", # ATTRIBUTES
"BIND",
"PROTECTED",
"VOLATILE",
# intrinsic functions
"ABS", "ACHAR", "ACOS", "ADJUSTL", "ADJUSTR", "AIMAG", "AINT",
"ALL", "ALLOCATED", "ANINT", "ANY", "ASIN", "ASSOCIATED",
"ATAN", "ATAN2", "BIT_SIZE", "BTEST", "CEILING", "CHAR", "CMPLX",
"CONJG", "COS", "COSH", "COUNT", "CSHIFT", "DATE_AND_TIME", "DBLE",
"DIGITS", "DIM", "DOT_PRODUCT", "DPROD", "EOSHIFT", "EPSILON",
"EXP", "EXPONENT", "FLOOR", "FRACTION", "HUGE", "IACHAR", "IAND",
"IBCLR", "IBITS", "IBSET", "ICHAR", "IEOR", "INDEX", "INT", "IOR",
"ISHFT", "ISHFTC", "KIND", "LBOUND", "LEN", "LEN_TRIM", "LGE", "LGT",
"LLE", "LLT", "LOG", "LOG10", "LOGICAL", "MATMUL", "MAX",
"MAXEXPONENT", "MAXLOC", "MAXVAL", "MERGE", "MIN", "MINEXPONENT",
"MINLOC", "MINVAL", "MOD", "MODULO", "MVBITS", "NEAREST", "NINT",
"NOT", "PACK", "PRECISION", "PRESENT", "PRODUCT", "RADIX",
"RANDOM_NUMBER", "RANDOM_SEED", "RANGE",
"REPEAT", "RESHAPE", "RRSPACING", "SCALE", "SCAN",
"SELECTED_INT_KIND", "SELECTED_REAL_KIND", "SET_EXPONENT",
"SHAPE", "SIGN", "SIN", "SINH", "SIZE", "SPACING", "SPREAD", "SQRT",
"SUM", "SYSTEM_CLOCK", "TAN", "TANH", "TINY", "TRANSFER",
"TRANSPOSE", "TRIM", "UBOUND", "UNPACK", "VERIFY",
# F95 intrinsic functions.
"NULL", "CPU_TIME",
# F2003.
"MOVE_ALLOC", "COMMAND_ARGUMENT_COUNT", "GET_COMMAND",
"GET_COMMAND_ARGUMENT", "GET_ENVIRONMENT_VARIABLE",
"SELECTED_CHAR_KIND", "WAIT", "FLUSH", "NEW_LINE",
"EXTENDS", "EXTENDS_TYPE_OF", "SAME_TYPE_AS", "BIND",
# F2003 ieee_arithmetic intrinsic module.
"IEEE_SUPPORT_UNDERFLOW_CONTROL", "IEEE_GET_UNDERFLOW_MODE",
"IEEE_SET_UNDERFLOW_MODE",
# F2003 iso_c_binding intrinsic module.
"C_LOC", "C_FUNLOC", "C_ASSOCIATED", "C_F_POINTER",
"C_F_PROCPOINTER"
]


splitSeparators=[
' ',
',',
'=',
'+',
'-']

separators=[
' ',
',',
'=',
'+',
'-',
'*',
'/',
'(',
')',
'%',
'.',
'>',
'<',
':',
'!',
'&',
'"',
"'"]

newLevelKeyword=[
'PROGRAM',
'MODULE',
'SUBROUTINE',
'FUNCTION',
'INTERFACE',
'THEN',
'DO',
'WHERE']

endLevelKeyword=[
'END',
'ENDIF',
'ENDDO',
'ENDFUNCTION',
'ENDWHERE']

noLevelKeyword=[
"PROCEDURE"]

# wordtypes
is_empty           = 0
is_string          = 1
is_comment         = 2
is_fortran_keyword = 3
is_variable        = 4
is_separator       = 5
is_spaces          = 6
is_number          = 7

stringDelimeters=["'",'"']
stringDelim="'\""
spaces="                                                                                          "

wordsListCounter   = 0
wordsList = [' '] * 1000
wordsListReplace = [' '] * 1000


def warning(line,message):
  print '********************'
  print lineNumber,": ",line
  print "Warning: "+message
  print '********************'
  
def error(line,message):
  print '********************'
  print lineNumber,": ",line
  print "ERROR: "+message
  print '********************'
  sys.exit(1)

def exceedsNameLength(word):
  warning(word, ': word exceeds maxNameLength')
  # check if we have already written it
  if (sedScript):
    scriptfile = open(sedScriptFile, 'r')
    filetext = scriptfile.read()
    scriptfile.close()
    if (filetext.find(word) < 0):
      scriptfile = open(sedScriptFile, 'a')
      scriptfile.write("echo 's/"+word+"/newName/g' >> tt\n")
      scriptfile.close()

  return
  
def findEndOfString(inLine, delim, startCheck):
  global stringContinousDelim
  
  k =  inLine.find(delim, startCheck)
  if (k < 0):
    endOfLine=len(inLine)-1
    if (not(inLine[endOfLine] == '&')):
      error(inLine,"findEndOfString unmatched string delimiter")
    else:
      stringContinousDelim = delim
      return (endOfLine, '&')
  return(k, delim)
    

def fortranKeyword(inWord):
  global lineHasNewLevel
  global lineHasEndLevel
  global lineHasSelectKeyword
  global lineHasCaseKeyword
  global lineHasElseKeyword
  global lineHasTypeKeyword
  global lineHasPublicKeyword
  global lineHasContainsKeyword
  global lineHasFlushKeyword
  global lineHasWhereKeyword
  global lineHasNoLevel
  
  word = inWord.upper()
  if (word in fortran_keywords):
      # we have a fortran keyword
      # check the tabbing level
      if (word in newLevelKeyword):
          lineHasNewLevel = True
      if (word in endLevelKeyword):
          lineHasEndLevel = True
      if (word in noLevelKeyword):
          lineHasNoLevel = True

      # special tabbing for CASE
      if ("SELECT" == word):
        lineHasSelectKeyword = True
      elif ("CASE" == word):
        lineHasCaseKeyword = True
      # special tabbing for else
      elif ("ELSE" == word):
        lineHasElseKeyword = True
      elif ("ELSEIF" == word):
        lineHasElseKeyword = True
      elif ("CONTAINS" == word):
        lineHasContainsKeyword = True
      # special tabbing for Type
      elif ("TYPE" == word):
        lineHasTypeKeyword = True
      elif ("PUBLIC" == word):
        lineHasPublicKeyword = True
      elif ("FLUSH" == word):
        lineHasFlushKeyword = True
      if ("WHERE" == word):        
        lineHasWhereKeyword = True
        
      if (caseCheck):
        return (True, word)
      
      return (True, inWord)
        
  return (False, inWord)

  
def processVariableName(word):
  n=len(word)
  if (n < 1):
    return word
  if (caseCheck):
    return word.lower()
  
  return word

  # this was used only for CamelCase conversion
  if (word.isupper()):
    return word.lower()
  if (word.islower()):
    return word

  newName = word[0].lower()
  for i in range(1,n):
    if (word[i-1].islower() and word[i].isupper()):
      newName += "_"
    newName += word[i].lower()
  
  return newName
  


# main processing method
def processNextWord(inLine):
    
  global lineHasParenthesis
  
  word=''
  remainingLine = ''
  sep=''
  leadingSpaces = 0
  
  workLine = inLine.lstrip()
  workLineLength = len(workLine)
  leadingSpaces = len(inLine) - workLineLength

  if (workLineLength < 1):
    # this should not happen
    print 'Warning: word is empty:'+word+'->'+workLine
    wordType = is_empty
    leadingSpaces = 0
    return (word, remainingLine, sep, wordType, leadingSpaces)

  firstChar = workLine[0]
  # check if it is a comment, not used
  if (firstChar == '!'):
    # this should not happen
    print 'Warning: word is comment:'+word+'->'+workLine
    word = workLine
    remainingLine = ''
    wordType = is_comment
    return (word, remainingLine, sep, wordType, leadingSpaces)
  
  # check if we are in a string
  if ( stringDelim.find(firstChar) >= 0) :
    (n, sep) = findEndOfString(workLine, firstChar, 1)
    #  n = workLine.find(firstChar, 1)
    word = workLine[0:n]
    remainingLine = workLine[n+1:]
    wordType = is_string
    return (word, remainingLine, sep, wordType, leadingSpaces)


  # get next word
  n = len(workLine) + 1
  for sep in separators:
    k = workLine.find(sep)
    if ( k >= 0):
      n = min(n,k)

  if (n >= len(workLine)):
    # just one remaining word
    sep=''
    word = workLine
    remainingLine = ''
  else:
    sep = workLine[n]
    word = workLine[0:n]
    remainingLine = workLine[n+1:]
  
  if ("(" == sep):
    lineHasParenthesis = True
  
  if (len(word) < 1):
    # its just the separator
    wordType = is_separator
    return (word, remainingLine, sep, wordType, leadingSpaces)
 
  if (word.isspace()):
    # this should not happen
    print 'Warning: word is spaces:'+word+'->'+workLine
    wordType = is_spaces
    return (word, remainingLine, sep, wordType, leadingSpaces)


  (isFortranKeyword,word) = fortranKeyword(word)

  wordLength = len(word)
  if (isFortranKeyword):
    wordType = is_fortran_keyword
    if (' ' == sep and indentCheck):
      # make sure we have exaclty one space
      remainingLine = remainingLine.lstrip()
  elif (word[0].isdigit()):
    wordType = is_number  
  else:
    wordType = is_variable
    word = processVariableName(word)      
    if (len(word) > maxNameLength):
      print word[wordLength-3:]
      exceedsNameLength(word)

  return (word, remainingLine, sep, wordType, leadingSpaces)
# end of processNextWord()

def findClosedParenthesis(line):
  k0 = findNotInString(line, '(')
  if (k0 < 0):
    return k0

  k0 = k0 + 1
  level = 1
  while (level > 0):
    k1 = line.find('(',k0)
    k2 = line.find(')',k0)
    #k1 = findNotInString(line[k0:], '(')
    #k2 = findNotInString(line[k0:], ')')
    if (k2 < 0):
      # not closed parenthesis
      warning(line, "not closed parenthesis")
      return -1
    
    if (k1 < 0):
      level = level - 1
      k0 = k2 + 1
      
    elif(k2 < k1):
      level = level - 1
      k0 = k2 + 1
      
    else:
      level = level + 1
      k0 = k1 + 1
    # print k0,k1,k2
  return k0          


def findNotInString(inLine,findString):

  # print 'Checking:'+inLine
  inLineLength = len(inLine)
  n = inLineLength + 1
  for delim in stringDelimeters:
    # print 'checking for ',delim 
    k = inLine.find(delim)
    if (k >= 0):
      n = min(n,k)
      
  if (n > inLineLength):
    endCheck = inLineLength
  else:
    endCheck = n
          
  # print 'endCheck=',endCheck
  # print 'find String in:',inLine[0:endCheck-1]
  k =  inLine.find(findString, 0, endCheck)
  if (k >=0 or endCheck >= inLineLength):
    return k
  
  delim = inLine[endCheck]
  # print 'delim is:'+delim
  (k, sep) = findEndOfString(inLine, delim, endCheck+1)
  if (k >= len(inLine)-1):
    return -1

  # print 'next line is is:'+inLine[k+1:]
  n = findNotInString(inLine[k+1:],findString)
  # print 'n is=',n
  if (n < 0):
    return -1
  return n+k+1


def getLineComments(inLine):
  global lineHasCompileDirecitve
  
  # get rid of the commends
  n = findNotInString(inLine,"!")
  if ( n < 0):
    return (inLine, '', 0)
  
  remainingLine = inLine[0:n]
  remainingLine = remainingLine.rstrip()
  commentLine   = inLine[n:]
  commentSpaces = len(inLine) - (len(remainingLine) + len(commentLine))
  # check if comment is an omp directive
  if (n==0):
    commLength=len(commentLine) - 1
    for commWord in compile_comment_directives:
      keylen=len(commWord)
      if (commLength >= keylen):
        chk_keyword=commentLine[1:keylen+1].upper()
        if (chk_keyword == commWord):
          lineHasCompileDirecitve = True
          
  #if (n == 0 and len(commentLine) > OMPKeywordLength):
    #check_omp=commentLine[1:OMPKeywordLength+1].upper()
    #if (check_omp == OMPKeyword):
      #lineHasCompileDirecitve = True
  #if (n == 0 and len(commentLine) > POMPKeywordLength):
    #check_omp=commentLine[1:POMPKeywordLength+1].upper()
    #if (check_omp == POMPKeyword):
      #lineHasCompileDirecitve = True
    
  
  # print commentLine
  return (remainingLine, commentLine, commentSpaces)

# the main driving method
def initNewLine():
  global lineHasParenthesis
  global lineHasNewLevel
  global lineHasEndLevel
  global lineHasSelectKeyword
  global lineHasCaseKeyword
  global lineHasElseKeyword
  global lineHasContainsKeyword
  global lineHasTypeKeyword
  global lineHasPublicKeyword
  global lineHasFlushKeyword
  global lineHasCompileDirecitve
  global lineHasWhereKeyword
  global lineHasNoLevel
   
  lineHasNewLevel = False
  lineHasEndLevel = False
  lineHasSelectKeyword = False
  lineHasCaseKeyword = False
  lineHasElseKeyword = False
  lineHasContainsKeyword = False
  lineHasTypeKeyword = False
  lineHasPublicKeyword = False
  lineHasParenthesis = False
  lineHasFlushKeyword = False
  lineHasCompileDirecitve = False
  lineHasWhereKeyword = False
  lineHasNoLevel = False
  

#def splitCommentLine(commentLine)

  #line1 = ''
  #line2 = ''
  #n = len(workLine) + 1
  #for sep in separators:
    #k = workLine.find(sep)
    #if ( k >= 0):
      #n = min(n,k)
   
   
   #if (len(outLine) == 0):
     #return splitCommentLine(commentLine)



def splitLine(outLine, commentLine):
  line1 = commentLine
  line2 = outLine
      
  return (line1, line2)



# the main driving method
def processSourceFortran(inName,outName):

  global wordsListCounter
  global lineHasTypeKeyword
  
  print "------------------------------------"
  print "proccessing ", inName, "..."
  infile  = open(inName, 'r')
  if (not noWrite):
    outfile = open(outName, 'w')

  global lineNumber
  global stringContinousDelim
  
  prefixLength = len(prefix)
  lineNumber = 0
  caseRepeats    = 0
  caseLevels = 0
  indentLevel = 0
  sep=''
  addDummy=''
  inContinuation = False
  fileChanged = False
  
  for inLine in infile:
     lineNumber = lineNumber + 1
     stringContinousDelim = ''
     if (not(inContinuation)):
       initNewLine()
     hasPreprocessKeyword = False
     
     # replace tabs with spaces
     line = inLine.expandtabs(noOfTabSpaces)
     # strip leading,trailing spaces
     line = line.rstrip()
     
     if (indentCheck):
       line = addDummy+line.lstrip()
     addDummyLength=len(addDummy)
     
     outLine = ''

     (remainingLine,commentLine,commentSpaces) = getLineComments(line)

     if (len(remainingLine) > 0):
       # see if it's a preprocessing directive
       if ( "#" == remainingLine[0]):
         hasPreprocessKeyword = True
         outLine   = remainingLine
         remainingLine = ""
       # check if we are in continaution line
       elif (inContinuation and indentCheck):
         if ( "&" == remainingLine[addDummyLength]):
           remainingLine = addDummy+remainingLine[addDummyLength+1:].lstrip()

     insert_prefix = False
     while (len(remainingLine) > 0):
       (word, remainingLine, sep, wordType, leadingSpaces) = processNextWord(remainingLine)
 
       if (add_prefix):
           if (word.upper() == keyword):
             insert_prefix=True
           elif (insert_prefix and len(word) > 0  and wordType != is_fortran_keyword ):
             insert_prefix=False
             if (word[0:prefixLength] != prefix):
               if (word not in wordsList):
                 print 'adding ', word
                 wordsList[wordsListCounter] = word
                 wordsListReplace[wordsListCounter] = prefix+word
                 wordsListCounter=wordsListCounter+1

#       if (replaceWord and len(word) > 0 and (lineHasTypeKeyword or lineHasPublicKeyword)):
       if (replaceWord and len(word) > 0):
         for i in range(wordsListCounter):
           if (word == wordsList[i]):
             print word, wordsListReplace[i]
             word = wordsListReplace[i]
             break
       if (leadingSpaces > 0) :
         outLine += spaces[0:leadingSpaces]
       outLine += word + sep

     if (noWrite):
       continue
    
     #end while (len(remainingLine) > 0):
     # process line is done

     if (lineHasFlushKeyword):
       warning(outLine,"line with flush keyword")

     if (indentCheck):
       hasNewLevel = lineHasNewLevel
       hasEndLevel = lineHasEndLevel
       if (lineHasNoLevel):
         hasNewLevel = False
         hasEndLevel = False

       if (lineHasWhereKeyword):
        l = findClosedParenthesis(outLine)
        if (l < len(outLine)):
          hasNewLevel = False

       if (lineHasCaseKeyword):
         caseRepeats = caseRepeats + 1
         if (lineHasSelectKeyword):
            caseLevels = caseLevels+1
         elif(caseRepeats > 2):
           hasEndLevel = True

       if (lineHasElseKeyword or lineHasContainsKeyword):
         hasEndLevel = True

       if (hasEndLevel and not(inContinuation)):
         if (indentLevel < 1):
           error(inLine, "unmatched indentLevel")
         else:
           indentLevel =  indentLevel - 1

       if (hasPreprocessKeyword or lineHasCompileDirecitve):
         indent = ''
       elif(inContinuation and len(outLine) > 0):
         indent=spaces[0:(indentLevel+1)*noOfTabSpaces] + '&'+ spaces[0:noOfTabSpaces-1]
       else:
         indent=spaces[0:indentLevel*noOfTabSpaces]
     
       writeLine = indent + outLine[addDummyLength:] + spaces[0:commentSpaces]+commentLine+'\n'

     else:
       writeLine = outLine + spaces[0:commentSpaces]+commentLine+'\n'

     if (len(writeLine) > maxLineLength):
       if (len(outLine) > 0):
         warning(writeLine,"line too long")
       #(line1,line2) = splitLine(outLine, commentLine)
       #if (len(line1) > 0):
         #outfile.write(indent+line1+'\n')
       #if (len(line2) > 0):
         #outfile.write(indent+line2+'\n')
     #else:
     if (not(writeLine == inLine)):
       fileChanged = True
       if (verbose):
         print lineNumber,":\n","<",inLine,">",writeLine
       
     outfile.write(writeLine)
     
     #if (writeLine != inLine):
     #  print inLine
     #  print writeLine
     
     if (indentCheck):
     
      if (lineHasElseKeyword or lineHasContainsKeyword):
        hasNewLevel = True
        hasEndLevel = False
            
      if (lineHasTypeKeyword and not(hasEndLevel) and not(lineHasParenthesis)):
        hasNewLevel = True
      
      if (lineHasCaseKeyword and caseRepeats > 1):
        hasNewLevel = True
        hasEndLevel = False
      
      inContinuation = (sep == '&')
            
      if ( not(hasEndLevel) and not(inContinuation) and hasNewLevel):
        indentLevel =  indentLevel + 1
  #     outfile.write(str(hasNewLevel)+"-"+str(hasEndLevel)+"->"+str(inContinuation)+str(indentLevel)+"\n")
      
      if (hasEndLevel and lineHasSelectKeyword):
        caseLevels = caseLevels-1
        if (caseLevels < 1):
          caseRepeats = 0      
      
      if (len(stringContinousDelim) > 0 ):
        addDummy=stringContinousDelim
      else:
        addDummy=''
          
        
     #s = line.split()
     #for i in range(len(s)):
     #  print i, s[i]
  if (not noWrite):
    outfile.close()
    if (indentLevel > 0):
      error('end of program','Wrong indentation')
      
    if (fileChanged):
      print outName, ' is changed.','\n'
    else:
      print outName, ' is not changed.','\n'
      
  return 0


def usage():
  print sys.argv[0], "[-v] [-h] [-c --casecheck (default), --no_casecheck] [-i --indentcheck] [-s --sed <sed_script_file>] <inputfile> <outputfile> "

def processArgs():
  global verbose
  global sedScript, sedScriptFile
  global caseCheck, indentCheck
  global noWrite
  global add_prefix
  global prefix
  global keyword
  global printfileName
  global filelistName
  global replaceWord
  
  caseCheck = True
  verbose = False
  sedScript = False
  indentCheck = False
  add_prefix  = False
  sedScriptFile=''
  noWrite = False
  replaceWord = False
  prefix=""
  keyword=""
  printfileName=''
  filelistName=''
  
  
  flags="vshicnr"
  longFlags=["verbose","casecheck", "no_casecheck", "indentcheck",
  'no_write', "keyword=", 'printfile=', 'filelist=', "add_trail_prefix="]
  try:
    opts, args = getopt.getopt(sys.argv[1:], flags, longFlags)
  except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()
    sys.exit(2)
      
  #if (not (len(args) == 2)):
    #usage()
    #sys.exit(3)

  
  for o, a in opts:
    if o == "-v":
       verbose = True
    elif o in ("-h", "--help"):
       usage()
       sys.exit()
    elif o in ("-s", "--sed"):
       sedScriptFile = a
       sedScript = True
    elif o in ("-n", "--no_write"):
       noWrite = True
    elif o in ("-c", "--casecheck"):
       caseCheck = True
    elif o in ("--no_casecheck"):
       caseCheck = False
    elif o in ("-i", "--indentcheck"):
       indentCheck = True
    elif o in ("-r", "--replace"):
       replaceWord = True
    elif o in ( "-v", "--verbose" ):
       verbose = True
    elif o in ( "--keyword" ):
       add_prefix = True
       keyword = a
    elif o in ( "--printfile" ):
       printfileName = a
    elif o in ( "--filelist" ):
       filelistName = a
 #      print 'keyword=',keyword
    elif o in ( "--add_trail_prefix" ):
       add_prefix = True
       prefix = a
#       print 'prefix=',prefix
    else:
       assert False, "unhandled option"

  if (len(args) > 1):
    return args[0], args[1]
  else:
    return '', ''


#print "Python version", sys.version
def main_method():

  global wordsListCounter
  wordsListCounter=0
  
  infile, outfile = processArgs()

  if (replaceWord):
    printfile = open(printfileName, 'r')
    for line in printfile:
      arg_list = line.split()
      wordsList[wordsListCounter] = arg_list[0]
      wordsListReplace[wordsListCounter] = arg_list[1]
      wordsListCounter = wordsListCounter+1
    printfile.close()
  
  if (filelistName != ''):
    filelist = open(filelistName, 'r')
    for line in filelist:
      arg_list = line.split()
      print arg_list
      processSourceFortran(arg_list[0], arg_list[1])
    filelist.close()
  else:    
    processSourceFortran(infile, outfile)

  if (printfileName != ''):
    printfile = open(printfileName, 'a')
    for i in range(wordsListCounter):
      if (len(wordsListReplace[i]) > maxNameLength):
        printfile.write(wordsList[i]+' '+wordsListReplace[i]+'  !@! \n')
      else:
        printfile.write(wordsList[i]+' '+wordsListReplace[i]+'\n')
    printfile.close()

  return

main_method()
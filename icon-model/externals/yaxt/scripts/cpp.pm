# cpp.pm --- routines to emulate C preprocessor
#
# Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
#
# Version: 1.0
# Author: Thomas Jahns <jahns@dkrz.de>
# Keywords:
# Maintainer: Thomas Jahns <jahns@dkrz.de>
# URL: https://swprojects.dkrz.de/redmine/projects/scales-ppm
#
# Redistribution and use in source and binary forms, with or without
# modification, are  permitted provided that the following conditions are
# met:
#
# Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# Neither the name of the DKRZ GmbH nor the names of its contributors
# may be used to endorse or promote products derived from this software
# without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
package cpp;

use strict;
use warnings;

require Exporter;
our @ISA=('Exporter');
our @EXPORT_OK = qw(transliterateTrigraphs tokenize);

use Data::Dumper ();
use File::Temp ();

our ($debug, $operator_match, $integral_constant_match, $definition_match);

# defaults for LP64
my ($c_sizeof_int, $c_sizeof_long, $c_sizeof_long_long,
    $c_char_is_unsigned) = (4, 8, 8, 0);

BEGIN {
  if (exists($ENV{'DEBUG_CPP'}) && $ENV{'DEBUG_CPP'} =~ m{^\d+$}) {
    $debug = $ENV{'DEBUG_CPP'} + 0;
  } else {
    $debug = 0;
  }
};

BEGIN {

   my %trigraphs = ('=' => '#', '/' => '\\', "'" => '^', '(' => '[',
                    ')' => ']', '!' => '|', '<' => '{', '>' => '}', '-' => '~');
   my $trigraphsyms = quotemeta(join('', keys %trigraphs));
   my $trigraphre = qr{\?\?([$trigraphsyms])}xo;

   sub transliterateTrigraphs
   {
      local $_;
      my ($text) = @_;
      print STDERR 'transliterating "', $text, '"', $trigraphre, "\n"
          if ($debug > 4);
      while ($text =~ s{$trigraphre}{$trigraphs{$1}})
      {
      }
      return $text;
   }

   sub validPreProcSym
   {
      return scalar ($_[0] =~ m{[a-zA-Z_]\w*});
   }

   $operator_match = qr{(?:==
                        |!=
                        |[><]=
                        |[-+/&^%|]=
                        |(?:\|\||\&\&|>>|<<)=
                        |(?:\|\||\&\&|>>|<<)
                        |\(
                        |\)
                        |\#\#
                        |[-,+\*/%&|^!=><])}xo;

   $integral_constant_match
       = qr{^(\d+|0[0-7]+|0[xX][0-9a-fA-F]+)([ul]{1,3})?$}xio;

   $definition_match
       = qr{^\s*\#+(?:define)
               \s+([A-Za-z_]\w*)
               (\(\s*[A-Za-z_]\w*(?:\s*,[A-Za-z_]\w*)*\s*\))?
               \s+(.*?)\s*$}xo;

   sub tokenize
   {
      print(STDERR 'In ', whoami(), ', line ', __LINE__, "\n") if $debug > 2;
      local $_;
      my $expression = $_[0];
      my @tokens;
      if ($expression =~ m{^([^"']*)(".*?(?<=[^\\])"|'.*?(?<=[^\\])')}go)
      {
         my ($prefix, $stringOrCharConst,
             $stringOrCharConstPos, $stringOrCharConstEnd)
             = ($1, $2,  $-[2], $+[2]);
         die('invalid character constant in expression >', $expression, '<')
             if ($stringOrCharConst eq "''");
         $expression =~ m{\G.*}go;
         my ($suffix, $suffixEnd) = ($&, pos($expression));
         print(STDERR
               Data::Dumper->Dump([\$prefix, \$stringOrCharConst, \$suffix],
                                  ['prefix', 'stringOrCharConst', 'suffix']),
               "\n") if ($debug > 1);
         my ($prefixTokenList, $prefixTokenPosMap) = tokenize($prefix);
         my ($suffixTokenList, $suffixTokenPosMap) = tokenize($suffix);
         @tokens = ([@$prefixTokenList, $stringOrCharConst, @$suffixTokenList],
                    [@$prefixTokenPosMap, [$stringOrCharConstPos, $stringOrCharConstEnd],
                     map { [$_->[0] + $stringOrCharConstEnd + 1,
                            $_->[1] + $stringOrCharConstEnd + 1] }
                     @$suffixTokenPosMap]);
         print(STDERR Data::Dumper->Dump([\@tokens], ['tokens']), "\n")
             if ($debug > 1);
      } elsif ($expression =~ m{["']}) {
         die('Unbalanced quote " in expression: ', $expression);
      } elsif ($expression ne '') {
         my (@tokenList, @tokenMappings);
         while ($expression =~ m{(.*?)
                                 (?=\s*$operator_match
                                 |\s+
                                 |$)}xgco) {
            if (defined($1) and $1 ne '')
            {
               my $tokenEndPos = pos($expression);
               push(@tokenList, $1);
               push(@tokenMappings, [$tokenEndPos - length($1), $tokenEndPos - 1]);
               print(STDERR 'Found ', $1, ' at ', $tokenMappings[-1][0],
                     '..', $tokenMappings[-1][1], ': `',
                     substr($expression, $tokenMappings[-1][0],
                            $tokenMappings[-1][1] - $tokenMappings[-1][0] + 1),
                     "'\n") if ($debug > 2);
            }
            $expression =~ m/\s*($operator_match)|\s+|$/xgco;
            if (defined($1) and $1 ne '') {
               my $tokenEndPos = pos($expression);
               push(@tokenList, $1);
               push(@tokenMappings, [$tokenEndPos - length($1), $tokenEndPos]);
            }
         }
         @tokens = (\@tokenList, \@tokenMappings);
         print(STDERR 'Tokenized \'', $expression, '\' as ',
               join("\n", @tokenList), "\n") if $debug > 1;
      } else {
         @tokens = ([], []);
      }
      return (@tokens);
   }

}

sub parseDefinition
{
   my ($argString, $definition, $file, $line) = @_;
   my (@tokenizedDef);
   if (defined($argString))
   {
      $argString =~ s{^\s+}{};
      $argString =~ s{\s+$}{};
   }
   else { $argString = ''; }
   my @tokens = tokenize($definition);
   if (@tokens)
   {
      die(defined($file)?($file, ':', $line):('cmd arg ', $line),
          'error: \'##\' must not appear at either end or start of a macro',
          ' expansion')
          if ($tokens[0] eq '##' or $tokens[-1] eq '##');
   }
   return {
           'parameters' => [split(/\s*,\s*/, $argString)],
           'definition' => \@tokens
          };
}

sub evaluatePreProcExpression
{
   local $_;
   my ($expression, $defines, $cppKeyOut, $dontEvaluate) = @_;
   my @ast = parsePreProcExpression($expression, $defines, $cppKeyOut);
   print(STDERR 'DontEvaluate: \'', (!defined($dontEvaluate))?
         'undef':$dontEvaluate, "'\n") if $debug > 2;
   return (defined($dontEvaluate) && $dontEvaluate) ?
       1 : cpp_true(cpp_evaluate(@ast));
}

sub parsePreProcExpression
{
   print(STDERR 'In ', whoami(), ', line ', __LINE__, "\n") if $debug > 1;
   local $_;
   my ($expression, $defines, $cppKeyOut) = @_;
   my @tokenMapping = tokenize($expression);
   my @tokens = @{$tokenMapping[0]};
   @tokenMapping = @{$tokenMapping[1]};
   my $didSubst = 1;
   my $i;
   for ($i = 0; $i < @tokens; ++$i)
   {
      print(STDERR 'In ', whoami(), ', line ', __LINE__, "\n") if $debug > 1;
      my $token = $tokens[$i];
      print(STDERR 'Looking at token \'', $token, "'\n") if $debug > 2;
      if ($token eq 'defined')
      {
         # the predicate defined causes the next symbol to be
         # swallowed and evaluates to 1 or 0 depending on whether the
         # symbol is defined or not
         my ($symreflen, $symoffset);
         if ($tokens[$i + 1] eq '(')
         {
           my $leadingParens = 1;
           while ($tokens[$i + $leadingParens + 1] eq '(')
           {
             ++$leadingParens;
           }
           $symreflen = 1 + $leadingParens * 2;
           $symoffset = $leadingParens + 1;
           for (my $j = 1; $j <= $leadingParens; ++$j)
           {
             #         defined +   number of (  + symbol + )
             my $ofs =    $i   + $leadingParens +    1   + $j;
             die('EOL in defined predicate') if (@tokens <= $ofs);
             die('Expected closing parenthesis, got ', $tokens[$ofs])
                 if ($tokens[$ofs] ne ')');
           }
         }
         else
         {
            $symreflen = 1;
            $symoffset = 1;
         }
         die('EOL in defined predicate') if ($i + $symreflen > @tokens);
         die('Invalid preprocessor symbol ', $tokens[$i + $symoffset])
             unless(isSymbol($tokens[$i + $symoffset]));
         insertDefinitionOnDemand($defines, $tokens[$i + $symoffset])
             if (!exists($defines->{$tokens[$i + $symoffset]}));
         splice(@tokens, $i, 1 + $symreflen,
                (defined($defines->{$tokens[$i + $symoffset]})?'1L':'0L'));
         next;
      }
      elsif (validPreProcSym($token))
      {
         print(STDERR 'In ', whoami(), ', line ', __LINE__, "\n") if $debug > 2;
         print(STDERR
               Data::Dumper->Dump([\$defines], ['defines']),
               "\n") if $debug > 4;

         push(@$cppKeyOut, $token) if (defined($cppKeyOut));
         insertDefinitionOnDemand($defines, $token);
         print(STDERR 'token ', $token,
               (defined($defines->{$token}) ?
                    ( ': found definition ',
                      join(' ', @{$defines->{$token}{'definition'}[0]}))
                    : ' undefined'), "\n") if $debug > 3;
         if (defined($defines->{$token}))
         {
            my @definition = @{$defines->{$token}{'definition'}};
            print(STDERR 'Replacing token \'', $token, '\' with ',
                  Data::Dumper->Dump([$definition[0]], [qw(definition)]))
                if $debug > 2;
            my @parameters = @{$defines->{$token}{'parameters'}};
            if (!@parameters)
            {
               splice(@tokens, $i, 1, @{$definition[0]});
               splice(@tokenMapping, $i, 1, @{$definition[1]});
               if ($i < @tokens) { redo; } else { last; }
            }
            elsif ($tokens[$i + 1] eq '(')
            {
               my @args;
               my (@argTokens);
               my $parencount = 1;
               my $j = $i + 2;
               while($j < @tokens and $parencount)
               {
                  if ($tokens[$j] eq ')')
                  {
                     last if !--$parencount;
                  }
                  elsif ($tokens[$j] eq '(')
                  {
                     ++$parencount;
                  }
                  elsif ($tokens[$j] eq ',' and $parencount == 1)
                  {
                     push(@args, [ @argTokens ]);
                     @argTokens = ();
                  }
                  push(@argTokens, $tokens[$j]);
                  ++$j;
               }
               die ('EOL in gathering macro arguments: ', @tokens[$i..$#tokens])
                   if ($parencount);
               die('Argument count mismatch: ', @parameters, ' expected, ',
                   @args, ' given.', "\n")
                   if @args != @parameters;
               my %argDefs = ( 'tokens' => map { $parameters[$_] => $args[$_]
                                              }
                               (0..$#args),
                               'tokenStream' => $expression,
                               'tokenMapping' => \@tokenMapping
                             );
               my @expandedDefinition
                   = expandMacro(\%argDefs, @tokens);
               splice(@tokens, $i, $j - $i + 1, @expandedDefinition);
               redo;
            }
         }
      }
   }
   # A12.5 macro replacement of undefined symbols
   @tokens = map { isSymbol($_)?'0L':$_ } @tokens;
   {
      print STDERR Data::Dumper->Dump([\@tokens, \@tokenMapping],
                                      ['tokens', 'tokenMapping'])
          if $debug > 1;
      my %lex = makeGetsym(@tokens);
      my @ast = conditional_expression(\%lex, 1);
      if (defined($lex{'currentsym'}()))
      {
         print(STDERR join(', ', $lex{'currentsym'}(), $lex{'rest'}()),
               "\n", Data::Dumper->Dump([\@ast], ['ast']), "\n")
             if $debug > 1;
         die('Parse of conditional incomplete: ', join($lex{'rest'}()));
      }
      print(STDERR Data::Dumper->Dump([\@ast], ['ast']))
          if $debug > 1;
      return @ast;
   }
}

sub makeGetsym
{
   my (@tokens) = @_;
   my $i = 0;
   return ('currentsym' => sub { return $i<@tokens?$tokens[$i]:undef; },
           'advancesym' => sub { return $i<@tokens?$#tokens - $i++:0; },
           'rest' => sub { return @tokens[$i..$#tokens] },
           'parse_state' => sub { return $i },
           'set_parse_state' => sub { $i = $_[0] },
           'consumed' => sub { my $state = $_[0]; return @tokens[$state..$i] });
}

sub expect
{
   my ($lex, $force, $expected) = @_;
   my $currentsym = $lex->{'currentsym'}();
   die('expect: unexpected symbol ', $currentsym,
       ' vs. expectation ', $expected)
       if ($currentsym ne $expected);
   return $lex->{'advancesym'}();
}

sub expandMacro
{
   my ($argDefs, @tokenMapping) = @_;
   print(STDERR 'In ', whoami(), ', line ', __LINE__, "\n") if $debug > 1;
   my @tokens = @{$tokenMapping[0]};
   @tokenMapping = @{$tokenMapping[1]};
   my $i;
   for ($i = 1; $i < @tokens - 1; ++$i)
   {
      if ($tokens[$i] eq '#'
          and exists($argDefs->{$tokens[$i+1]}))
      {
         splice(@tokens, $i, 2, '"' . $argDefs->{$tokens[$i+1]} . '"');
      }
      elsif ($tokens[$i] eq '##')
      {
         splice(@tokens, $i - 1, 3, @tokens[$i-1,$i+1]);
         redo;
      }
   }
   return (@tokens, @tokenMapping);
}

BEGIN {

  our (@preprocCmd, $tempDir);

   sub preprocResult
   {
      local $_;
      my ($inputSym) = @_;
      if (!@preprocCmd)
      {
        die('The preprocessing command is not set, but symbol "', $inputSym,
            ' requires', "\n",
            'invoking the Fortran preprocessor, please set the FPP',
            ' environment', "\n", 'variable to retry.');
      }
      print(STDERR 'In ', whoami(), ', line ', __LINE__, "\n") if $debug > 1;
      # FIXME: mode depends on file type used
      if (!defined($tempDir))
      {
        $tempDir = File::Temp::tempdir('CLEANUP' => 1);
      }
      my ($tempFileFh, $tempFileName)
          = File::Temp::tempfile('DIR' => $tempDir, 'SUFFIX' => '.f90');
      print($tempFileFh $inputSym, "\n")
          or die('Write to temporary file ', $tempFileName, ' failed!');
      close($tempFileFh)
          or die('Close of temporary file ', $tempFileName, ' failed!');
      my $preprocOutputFh;
      print(STDERR 'Opening ', join(' ', @preprocCmd, $tempFileName, '|'), "\n")
            if $debug > 1;
      my $fppin;
      my $child_pid = IPC::Open2::open2($preprocOutputFh, $fppin,
                                        @preprocCmd, $tempFileName);
      close($fppin) or die('Closing input to FPP child failed!');
      my $preprocOutput = '';
      my $p = '';
      while (<$preprocOutputFh>)
      {
         $p .= $_;
         $preprocOutput .= $_ if (rindex($_, '#', 0) && $_ ne "\n");
      }
      while (chomp($preprocOutput))
      {
      }
      my $pidmatch = waitpid($child_pid, 0);
      if ($pidmatch != $child_pid or $? != 0 or not close($preprocOutputFh))
      {
         die('Close of pipe to preprocessor failed!');
      }
      print(STDERR $tempFileName, "\n", $p, "\n")
          if $debug > 2;
      system("cat $tempFileName >&2")
          if $debug > 2;
      unlink($tempFileName)
          or die('Unlinking temporary file ', $tempFileName, ' failed!');
      print(STDERR "symbol '", $inputSym, "' replaced with '",
            $preprocOutput, "'\n") if $debug > 2;
      return $preprocOutput eq $inputSym ? undef : $preprocOutput;
   }

   sub insertDefinitionOnDemand
   {
      my ($defines, $token) = @_;
      print(STDERR 'In ', whoami(), ', line ', __LINE__, "\n") if $debug > 1;
      print(STDERR Data::Dumper->Dump([\$defines, \$token],
                                      ['defines', 'token']), "\n")
          if $debug > 4;
      if (exists($defines->{$token}))
      {
        print(STDERR 'Returning ',
              defined($defines->{$token}) ?
                  @{$defines->{$token}{'definition'}[0]}
                  : 'undef',
              ' for token ', $token, "\n") if $debug > 3;
        return $defines->{$token};
      }
      else
      {
         my $definition = preprocResult($token);
         if (defined($definition))
         {
           $defines->{$token} = $defines->{$token}
               = parseDefinition(undef, $definition, '', '');
         }
         else
         {
           $defines->{$token} = $defines->{$token} = undef;
         }
      }
   }

   sub parse_backoff
   {
      my ($lex, $save, $force, @msg) = @_;
      $lex->{'set_parse_state'}($save);
      if ($force)
      {
         die(@msg);
      }
      else
      {
         return ();
      }
   }

   my %evaluations;

   sub expression
   {
      local $_;
      my ($lex, $force) = @_;
      my $save = $lex->{'parse_state'}();
      print(STDERR 'In ', whoami(), ', left over: ',
            join(', ', $lex->{'rest'}()), "\n")
          if $debug > 1;
      my @expression = ('expression');
      my @assignment_expression = assignment_expression($lex, $force);
      if (@assignment_expression
          and (!defined($_ = $lex->{'currentsym'}())
               or $_ ne ','))
      {
         return @assignment_expression;
      }
      else
      {
         my @sub_expression = @assignment_expression
             ? @assignment_expression : expression($lex, $force);
         my $currentsym;
         if (@sub_expression
             and defined ($currentsym = $lex->{'currentsym'}())
             and $currentsym eq ','
             and $lex->{'advancesym'}())
         {
            print(STDERR 'In ', whoami(), ', line ', __LINE__, ', left over: ',
                  join(', ', $lex->{'rest'}()), "\n") if $debug > 1;
            if (@assignment_expression = assignment_expression($lex, $force))
            {
               return ('expression', \@sub_expression, \@assignment_expression);
            }
         }
      }
   }

   my @assignment_operators = ('=', '*=', '/=', '%=', '+=', '-=',
                               '<<=', '>>=', '&=', '^=', '|=');

   sub assignment_expression
   {
      local $_;
      my ($lex, $force) = @_;
      my $save = $lex->{'parse_state'}();
      print(STDERR 'In ', whoami(), ', left over: ',
            join(', ', $lex->{'rest'}()), "\n")
          if $debug > 1;
      my @conditional_expression = conditional_expression($lex, $force);
      if (@conditional_expression)
      {
         return @conditional_expression;
      }
      else
      {
         my @unary_expression = unary_expression($lex, $force);
         my ($currentsym, @assignment_expression);
         if (@unary_expression
             and defined ($currentsym = $lex->{'getcurrentsym'}())
             and (grep { $currentsym eq $_ } @assignment_operators)
             and $lex->{'advancesym'}()
             and (@assignment_expression = assignment_expression($lex,$force)))
         {
            print(STDERR 'In ', whoami(), ', line ', __LINE__, ', left over: ',
                  join(', ', $lex->{'rest'}()), "\n") if $debug > 1;
            return ('assignment_expression', $currentsym,
                    \@unary_expression, \@assignment_expression);
         }
      }
   }

   $evaluations{'conditional_expression'}
       = sub { return cpp_true(cpp_evaluate(@{$_[0]}))
                   ? (cpp_evaluate(@{$_[1]}))
                       : (cpp_evaluate(@{$_[2]})); };

   sub conditional_expression
   {
      local $_;
      my ($lex, $force) = @_;
      my $save = $lex->{'parse_state'}();
      print(STDERR 'In ', whoami(), ', left over: ',
            join(', ', $lex->{'rest'}()), "\n")
          if $debug > 1;
      my @logical_or_expression
          = logical_or_expression($lex, $force);
      my $currentsym;
      if (defined($currentsym = $lex->{'currentsym'}())
          and $currentsym eq '?')
      {
         my @elseexpression;
         my @ifexpression;
         if ($lex->{'advancesym'}()
             and (@ifexpression = expression($lex, $force))
             and defined($currentsym = $lex->{'currentsym'}())
             and $currentsym eq ':'
             and $lex->{'advancesym'}()
             and (@elseexpression = conditional_expression($lex, $force)))
         {
            return ('conditional_expression', @logical_or_expression,
                    @ifexpression, @elseexpression);
         }
         else
         {
            return parse_backoff($lex, $save, $force,
                                 'ternary operator \':\' not' .
                                 ' followed by else expression!');
         }
      }
      return @logical_or_expression;
   }

   $evaluations{'logical_or_expression'}
       = sub { my $accum = 0;
               foreach my $subexpression (@_) {
                  if (cpp_true(cpp_evaluate(@$subexpression))) {
                     $accum = 1;
                     last;
                  }
               }
               return ('constant', 'int', $accum);
            };

   sub logical_or_expression
   {
      local $_;
      my ($lex, $force) = @_;
      my $save = $lex->{'parse_state'}();
      print(STDERR 'In ', whoami(), ', left over: ',
            join(', ', $lex->{'rest'}()), "\n")
          if $debug > 1;
      my @logical_and_expression
          = logical_and_expression($lex, 0);
      if (@logical_and_expression
          and (!defined($lex->{'currentsym'}())
               or $lex->{'currentsym'}() ne '||'))
      {
         return @logical_and_expression;
      }
      elsif (@logical_and_expression)
      {
         my @logical_or_expression = ('logical_or_expression',
                                      [ @logical_and_expression ]);
         my $currentsym;
         while (defined($currentsym = $lex->{'currentsym'}())
                and $lex->{'currentsym'}() eq '||'
                and $lex->{'advancesym'}())
         {
            print(STDERR 'In ', whoami(), ', left over: ',
                  join(', ', $lex->{'rest'}()), "\n") if $debug > 1;
            if (@logical_and_expression = logical_and_expression($lex, $force))
            {
               push(@logical_or_expression, [ @logical_and_expression ]);
            }
            else
            {
               return parse_backoff($lex, $save, $force,
                                    'Invalid logical_or_expression ',
                                    $lex->{'consumed'}($save));
            }
         }
         return @logical_or_expression;
      }
      return parse_backoff($lex, $save, $force,
                           'Invalid logical_or_expression ',
                           $lex->{'consumed'}($save));
   }

   $evaluations{'logical_and_expression'}
       = sub { my $accum = 1;
               foreach my $subexpression (@_) {
                  print(STDERR 'In ', whoami(), ', result: ',
                        cpp_true(cpp_evaluate(@$subexpression)),
                        "\n") if $debug > 1;
                  if (!cpp_true(cpp_evaluate(@$subexpression))) {
                     $accum = 0;
                     last;
                  }
               }
               return ('constant', 'int', $accum);
            };

   sub logical_and_expression
   {
      local $_;
      my ($lex, $force) = @_;
      my $save = $lex->{'parse_state'}();
      print(STDERR 'In ', whoami(), ', left over: ',
            join(', ', $lex->{'rest'}()), "\n")
          if $debug > 1;
      my @inclusive_or_expression
          = inclusive_or_expression($lex, 0);
      if (@inclusive_or_expression
          and (!defined($lex->{'currentsym'}())
               or $lex->{'currentsym'}() ne '&&'))
      {
         return @inclusive_or_expression;
      }
      elsif (@inclusive_or_expression)
      {
         my @logical_and_expression = ('logical_and_expression',
                                       [ @inclusive_or_expression ]);
         my $currentsym;
         while (defined($currentsym = $lex->{'currentsym'}())
                and $currentsym eq '&&'
                and $lex->{'advancesym'}())
         {
            print(STDERR 'In ', whoami(), ', left over: ',
                  join(', ', $lex->{'rest'}()), "\n") if $debug > 1;
            if (@inclusive_or_expression
                = inclusive_or_expression($lex, $force))
            {
               push(@logical_and_expression, [ @inclusive_or_expression ] );
            }
            else
            {
               return parse_backoff($lex, $save, $force,
                                    'Invalid logical_and_expression ',
                                    $lex->{'consumed'}($save));
            }
         }
         return @logical_and_expression;
      }
      return parse_backoff($lex, $save, $force,
                           'Invalid logical_and_expression ',
                           $lex->{'consumed'}($save));
   }

   $evaluations{'inclusive_or_expression'}
       = sub
   {
      my @accum = cpp_evaluate(@{$_[0]});
      die('Invalid computation on non-constant: ', join(':', @accum))
          if ($accum[0] ne 'constant');
      foreach my $subexpression (@_[1..$#_]) {
         my @operand = cpp_evaluate(@$subexpression);
         my $promotion = cpp_arithmetic_conversion(\@accum, \@operand);
         die('Invalid operation | (OR) performed on float type ', $promotion)
             if (grep { $_ eq $promotion } ('float', 'double', 'long-double'));
         die('Invalid computation on non-constant: ',
             join(':', @operand))
             if ($operand[0] ne 'constant');
         $accum[1] = $promotion;
         $accum[2] |= $operand[2];
      }
      return @accum;
   };

   sub inclusive_or_expression
   {
      local $_;
      my ($lex, $force) = @_;
      my $save = $lex->{'parse_state'}();
      print(STDERR 'In ', whoami(), ', left over: ',
            join(', ', $lex->{'rest'}()), "\n")
          if $debug > 1;
      my @exclusive_or_expression
          = exclusive_or_expression($lex, $force);
      if (@exclusive_or_expression
          and (!defined($lex->{'currentsym'}())
               or $lex->{'currentsym'}() ne '|'))
      {
         return @exclusive_or_expression;
      }
      elsif (@exclusive_or_expression)
      {
         my @inclusive_or_expression = ('inclusive_or_expression',
                                        [ @exclusive_or_expression ] );
         my $currentsym;
         while (defined($currentsym = $lex->{'currentsym'}())
                and $currentsym eq '|'
                and $lex->{'advancesym'}())
         {
            print(STDERR 'In ', whoami(), ', left over: ',
                  join(', ', $lex->{'rest'}()), "\n") if $debug > 1;
            if (@exclusive_or_expression
                = exclusive_or_expression($lex, $force))
            {
               push(@inclusive_or_expression, [ @exclusive_or_expression ] );
            }
            else
            {
               return parse_backoff($lex, $save, $force,
                                    'Invalid inclusive_or_expression ',
                                    $lex->{'consumed'}($save));
            }
         }
         return @inclusive_or_expression;
      }
      return parse_backoff($lex, $save, $force,
                           'Invalid inclusive_or_expression ',
                           $lex->{'consumed'}($save));
   }

   $evaluations{'exclusive_or_expression'}
       = sub
   {
      my @accum = cpp_evaluate(@{$_[0]});
      die('Invalid computation on non-constant: ', join(':', @accum))
          if ($accum[0] ne 'constant');
      foreach my $subexpression (@_[1..$#_]) {
         my @operand = cpp_evaluate(@$subexpression);
         my $promotion = cpp_arithmetic_conversion(\@accum, \@operand);
         die('Invalid operation ^ (XOR) performed on float type ', $promotion,
             ': ', join(':', @accum), ' ^ ', join(':', @operand))
             if (grep { $_ eq $promotion } ('float', 'double', 'long-double'));
         die('Invalid computation on non-constant: ',
             join(':', @operand))
             if ($operand[0] ne 'constant');
         $accum[1] = $promotion;
         $accum[2] ^= $operand[2];
      }
      return @accum;
   };

   sub exclusive_or_expression
   {
      local $_;
      my ($lex, $force) = @_;
      my $save = $lex->{'parse_state'}();
      print(STDERR 'In ', whoami(), ', left over: ',
            join(', ', $lex->{'rest'}()), "\n")
          if $debug > 1;
      my @and_expression
          = and_expression($lex, $force);
      if (@and_expression
          and (!defined($lex->{'currentsym'}())
               or $lex->{'currentsym'}() ne '^'))
      {
         return @and_expression;
      }
      elsif (@and_expression)
      {
         my @exclusive_or_expression = ('exclusive_or_expression',
                                        [ @and_expression ] );
         my $currentsym;
         while (defined($currentsym = $lex->{'currentsym'}())
                and $currentsym eq '^'
                and $lex->{'advancesym'}())
         {
            print(STDERR 'In ', whoami(), ', left over: ',
                  join(', ', $lex->{'rest'}()), "\n") if $debug > 1;
            if (@and_expression = and_expression($lex, $force))
            {
               push(@exclusive_or_expression, [ @and_expression ] );
            }
            else
            {
               return parse_backoff($lex, $save, $force,
                                    'Invalid exclusive_or_expression ',
                                    $lex->{'consumed'}($save));
            }
         }
         return @exclusive_or_expression;
      }
      return parse_backoff($lex, $save, $force,
                           'Invalid exclusive_or_expression ',
                           $lex->{'consumed'}($save));
   }

   $evaluations{'and_expression'}
       = sub
   {
      my @accum = cpp_evaluate(@{$_[0]});
      die('Invalid computation on non-constant: ', join(':', @accum))
          if ($accum[0] ne 'constant');
      foreach my $subexpression (@{$_}[1..$#{$_}]) {
         my @operand = cpp_evaluate(@$subexpression);
         my $promotion = cpp_arithmetic_conversion(\@accum, \@operand);
         die('Invalid operation & performed on float type ', $promotion,
             ': ', join(':', @accum), ' ^ ', join(':', @operand))
             if (grep { $_ eq $promotion } ('float', 'double', 'long-double'));
         die('Invalid computation on non-constant: ',
             join(':', @operand))
             if ($operand[0] ne 'constant');
         $accum[1] = $promotion;
         $accum[2] &= $operand[2];
      }
      return @accum;
   };

   sub and_expression
   {
      local $_;
      my ($lex, $force) = @_;
      my $save = $lex->{'parse_state'}();
      print(STDERR 'In ', whoami(), ', left over: ',
            join(', ', $lex->{'rest'}()), "\n")
          if $debug > 1;
      my @equality_expression
          = equality_expression($lex, $force);
      if (@equality_expression
          and (!defined($lex->{'currentsym'}())
               or $lex->{'currentsym'}() ne '&'))
      {
         return @equality_expression;
      }
      elsif (@equality_expression)
      {
         my @and_expression = ('and_expression', [ @equality_expression ]);
         my $currentsym;
         while (defined($currentsym = $lex->{'currentsym'}())
                and $currentsym eq '&'
                and $lex->{'advancesym'}())
         {
            print(STDERR 'In ', whoami(), ', left over: ',
                  join(', ', $lex->{'rest'}()), "\n") if $debug > 1;
            if (@equality_expression = equality_expression($lex, $force))
            {
               push(@and_expression, [ @equality_expression ]);
            }
            else
            {
               return parse_backoff($lex, $save, $force,
                                    'Invalid and_expression ',
                                    $lex->{'consumed'}($save));
            }
         }
         return @and_expression;
      }
      return parse_backoff($lex, $save, $force,
                           'Invalid and_expression ',
                           $lex->{'consumed'}($save));
   }

   $evaluations{'equality_expression'}
       = sub
   {
      my @val = ([cpp_evaluate(@{$_[1]})], [cpp_evaluate(@{$_[2]})]);
      my $result = $val[0][2] == $val[1][2];
      if ($_[0] eq '!=')
      {
         $result = !$result;
      }
      return ('constant', 'int', $result);
   };

   sub equality_expression
   {
      local $_;
      my ($lex, $force) = @_;
      my $save = $lex->{'parse_state'}();
      print(STDERR 'In ', whoami(), ', left over: ',
            join(', ', $lex->{'rest'}()), "\n")
          if $debug > 1;
      my @relational_expression
          = relational_expression($lex, $force);
      if (@relational_expression
          and (!defined($lex->{'currentsym'}())
               or ($lex->{'currentsym'}() ne '=='
                   and $lex->{'currentsym'}() ne '!=')))
      {
         return @relational_expression;
      }
      elsif (@relational_expression)
      {
         my @equality_expression;
         my @equality_subexpression = @relational_expression;
         my $op;
         while (defined($op = $lex->{'currentsym'}())
                and ($op eq '==' || $op eq '!=')
                and $lex->{'advancesym'}())
         {
            print(STDERR 'In ', whoami(), ', left over: ',
                  join(', ', $lex->{'rest'}()), "\n") if $debug > 1;
            if (@relational_expression = relational_expression($lex, $force))
            {
               @equality_subexpression
                   = @equality_expression
                       = ('equality_expression', $op,
                          [ @equality_subexpression ],
                          [ @relational_expression ]);
            }
            else
            {
               return parse_backoff($lex, $save, $force,
                                    'Invalid equality_expression ',
                                    $lex->{'consumed'}($save));
            }
         }
         return @equality_expression;
      }
      return parse_backoff($lex, $save, $force,
                           'Invalid equality_expression ',
                           $lex->{'consumed'}($save));
   }

   my @relationalOps = ('>', '>=', '<', '<=');

   $evaluations{'relational_expression'}
       = sub
   {
      my @val = ([cpp_evaluate(@{$_[1]})], [cpp_evaluate(@{$_[2]})]);
      my $result;
      if ($_[0] eq '>')
      {
         $result = $val[0][2] > $val[1][2];;
      }
      elsif ($_[0] eq '!=')
      {
         $result = $val[0][2] != $val[1][2];;
      }
      elsif ($_[0] eq '>=')
      {
         $result = $val[0][2] >= $val[1][2];;
      }
      elsif ($_[0] eq '<')
      {
         $result = $val[0][2] < $val[1][2];;
      }
      elsif ($_[0] eq '<=')
      {
         $result = $val[0][2] <= $val[1][2];;
      }
      return ('constant', 'int', $result);
   };

   sub relational_expression
   {
      local $_;
      my ($lex, $force) = @_;
      my $save = $lex->{'parse_state'}();
      print(STDERR 'In ', whoami(), ', left over: ',
            join(', ', $lex->{'rest'}()), "\n")
          if $debug > 1;
      my @shift_expression
          = shift_expression($lex, $force);
      if (@shift_expression
          and (!defined($lex->{'currentsym'}())
               or !grep { $lex->{'currentsym'}() eq $_} @relationalOps))
      {
         return @shift_expression;
      }
      elsif (@shift_expression)
      {
         my @relational_expression;
         my $op;
         my @relational_subexpression = @shift_expression;
         while (defined($op = $lex->{'currentsym'}())
                and (grep { $_ eq $op } @relationalOps)
                and $lex->{'advancesym'}())
         {
            print(STDERR 'In ', whoami(), ', left over: ',
                  join(', ', $lex->{'rest'}()), "\n") if $debug > 1;
            if (@shift_expression = shift_expression($lex, $force))
            {
               @relational_subexpression =
                   @relational_expression = ('relational_expression', $op,
                                             [ @relational_subexpression ],
                                             [ @shift_expression ]);
            }
            else
            {
               return parse_backoff($lex, $save, $force,
                                    'Invalid relational_expression ',
                                    $lex->{'consumed'}($save));
            }
         }
         return @relational_expression;
      }
      return parse_backoff($lex, $save, $force,
                           'Invalid relational_expression ',
                           $lex->{'consumed'}($save));
   }

   $evaluations{'shift_expression'}
       = sub
   {
      local $_;
      my ($op, @val) = ($_[0],
                        [cpp_evaluate(@{$_[1]})], [cpp_evaluate(@{$_[2]})]);
      my $promotion = cpp_arithmetic_conversion($val[0], $val[1]);
      die('Invalid operation ', $op, ' performed on float type ', $promotion)
          if (grep { $_ eq $promotion } ('float', 'double', 'long-double'));
      die('Undefined operation ', $op, ' performed on with negative shift ',
          $val[1][2]) if ($val[1][2] < 0);
      if ($op eq '<<')
      {
         return ('constant', $val[0][1], $val[0][2] << $val[1][2]);
      }
      else
      {
         return ('constant', $val[0][1], $val[0][2] >> $val[1][2]);
      }
   };

   sub shift_expression
   {
      local $_;
      my ($lex, $force) = @_;
      my $save = $lex->{'parse_state'}();
      print(STDERR 'In ', whoami(), ', left over: ',
            join(', ', $lex->{'rest'}()), "\n")
          if $debug > 1;
      my @additive_expression
          = additive_expression($lex, $force);
      if (@additive_expression
          and (!defined($_ = $lex->{'currentsym'}())
               or ($_ ne '<<' and $_ ne '>>')))
      {
         return @additive_expression;
      }
      elsif (@additive_expression)
      {
         my @shift_expression;
         my @shift_subexpression
             = @additive_expression;
         my $op;
         while (defined($op = $lex->{'currentsym'}())
                and ($op eq '<<' || $op eq '>>')
                and $lex->{'advancesym'}())
         {
            print(STDERR 'In ', whoami(), ', left over: ',
                  join(', ', $lex->{'rest'}()), "\n") if $debug > 1;
            if (@additive_expression = additive_expression($lex, $force))
            {
               @shift_subexpression
                   = @shift_expression
                       = ('shift_expression', $op,
                          [ @shift_subexpression ],
                          [ @additive_expression ]);
            }
            else
            {
               return parse_backoff($lex, $save, $force,
                                    'Invalid shift_expression ',
                                    $lex->{'consumed'}($save));
            }
         }
         return @shift_expression;
      }
      return parse_backoff($lex, $save, $force,
                           'Invalid shift_expression ',
                           $lex->{'consumed'}($save));
   }

   $evaluations{'additive_expression'}
       = sub
   {
      local $_;
      my ($op, @val) = ($_[0], [cpp_evaluate(@{$_[1]})],
                        [cpp_evaluate(@{$_[2]})]);
      my $promotion = cpp_arithmetic_conversion($val[0], $val[1]);
      if ($op eq '+')
      {
         return ('constant', $promotion, $val[0][2] + $val[1][2]);
      }
      else
      {
         return ('constant', $promotion, $val[0][2] - $val[1][2]);
      }
   };

   sub additive_expression
   {
      local $_;
      my ($lex, $force) = @_;
      my $save = $lex->{'parse_state'}();
      print(STDERR 'In ', whoami(), ', left over: ',
            join(', ', $lex->{'rest'}()), "\n")
          if $debug > 1;
      my @multiplicative_expression
          = multiplicative_expression($lex, $force);
      if (@multiplicative_expression
          and (!defined($_ = $lex->{'currentsym'}())
               or ($_ ne '-' and $_ ne '+')))
      {
         return @multiplicative_expression;
      }
      elsif (@multiplicative_expression)
      {
         my @additive_expression;
         my @additive_subexpression = @multiplicative_expression;
         my $op;
         while (defined($op = $lex->{'currentsym'}())
             and ($op eq '+' || $op eq '-')
             and $lex->{'advancesym'}())
         {
            print(STDERR 'In ', whoami(), ', left over: ',
                  join(', ', $lex->{'rest'}()), "\n") if $debug > 1;
            if (@multiplicative_expression
                = multiplicative_expression($lex, $force))
            {
               @additive_subexpression
                   = @additive_expression
                       = ('additive_expression', $op,
                          [ @additive_subexpression ],
                          [ @multiplicative_expression ]);
            }
            else
            {
               return parse_backoff($lex, $save, $force,
                                    'Invalid additive_expression ',
                                    $lex->{'consumed'}($save));
            }
         }
         return @additive_expression;
      }
      return parse_backoff($lex, $save, $force,
                           'Invalid additive_expression ',
                           $lex->{'consumed'}($save));
   }

   my @multiplicativeOps = ('*', '/', '%');

   $evaluations{'multiplicative_expression'}
       = sub
   {
      local $_;
      my ($op, @val) = ($_[0], [cpp_evaluate(@{$_[1]})],
                        [cpp_evaluate(@{$_[2]})]);
      my $promotion = cpp_arithmetic_conversion($val[0], $val[1]);
      if ($op eq '*')
      {
         return ('constant', $promotion, $val[0][2] * $val[1][2]);
      }
      elsif ($op eq '/')
      {
         return ('constant', $promotion, $val[0][2] / $val[1][2]);
      }
      else
      {
         return ('constant', $promotion, $val[0][2] % $val[1][2]);
      }
   };

   sub multiplicative_expression
   {
      local $_;
      my ($lex, $force) = @_;
      my $save = $lex->{'parse_state'}();
      print(STDERR 'In ', whoami(), ', left over: ',
            join(', ', $lex->{'rest'}()), "\n")
          if $debug > 1;
      my @cast_expression
          = cast_expression($lex, $force);
      if (@cast_expression
          and (!defined($lex->{'currentsym'}())
               or !grep { $_ eq $lex->{'currentsym'}() } @multiplicativeOps))
      {
         return @cast_expression;
      }
      elsif (@cast_expression)
      {
         my @multiplicative_expression;
         my $op;
         my @multiplicative_subexpression = @cast_expression;
         if (defined($op = $lex->{'currentsym'}())
             and (grep { $_ eq $op } @multiplicativeOps)
             and $lex->{'advancesym'}())
         {
            print(STDERR 'In ', whoami(), ', left over: ',
                  join(', ', $lex->{'rest'}()), "\n") if $debug > 1;
            if (@cast_expression = cast_expression($lex, $force))
            {
               @multiplicative_subexpression
                   = @multiplicative_expression
                       = ('multiplicative_expression', $op,
                          [ @multiplicative_subexpression ],
                          [ @cast_expression ]);
            }
            else
            {
               return parse_backoff($lex, $save, $force,
                                    'Invalid multiplicative_expression ',
                                    $lex->{'consumed'}($save));
            }
         }
         return @multiplicative_expression;
      }
      return parse_backoff($lex, $save, $force,
                           'Invalid multiplicative_expression ',
                           $lex->{'consumed'}($save));
   }

   $evaluations{'cast_expression'}
       = sub
   {
      local $_;
      my ($conversion, @val) = ($_[0], cpp_evaluate(@{$_[1]}));
      return ('constant', $conversion, $val[2]);
   };

   sub cast_expression
   {
      local $_;
      my ($lex, $force) = @_;
      my $save = $lex->{'parse_state'}();
      print(STDERR 'In ', whoami(), ', left over: ',
            join(', ', $lex->{'rest'}()), "\n")
          if $debug > 1;
      my @unary_expression
          = unary_expression($lex, $force);
      if (@unary_expression)
      {
         return @unary_expression;
      }
      elsif ($lex->{'currentsym'}() eq '('
             and $lex->{'advancesym'}())
      {
         my @cast_expression;
         my @typename
             = typename($lex, $force);
         if (@typename
             and $lex->{'currentsym'}() eq ')'
             and $lex->{'advancesym'}())
         {
            print(STDERR 'In ', whoami(), ', left over: ',
                  join(', ', $lex->{'rest'}()), "\n") if $debug > 1;
            @unary_expression
                = unary_expression($lex, $force);
            if (@unary_expression)
            {
               @cast_expression = ('cast_expression',
                                   \@typename, \@unary_expression);
               return @cast_expression;
            }
         }
      }
      return parse_backoff($lex, $save, $force,
                           'Invalid cast_expression ',
                           $lex->{'consumed'}($save));
   }

   my @unaryOps = ('&', '*', '+', '-', '~', '!');

   $evaluations{'unary_expression'}
       = sub
   {
      local $_;
      my ($op, @val) = ($_[0], cpp_evaluate(@{$_[1]}));
      if ($op eq '-')
      {
         return ('constant', $val[1], -$val[2]);
      }
      elsif ($op eq '!')
      {
         return ('constant', $val[1], 0+!$val[2]);
      }
      elsif ($op eq '+')
      {
         return ('constant', $val[1], +$val[2]);
      }
      elsif ($op eq '~')
      {
         return ('constant', $val[1], ~$val[2]);
      }
      die('Unhandled unary expression: (', join(', ', @val), ")\n")
   };

   sub unary_expression
   {
      local $_;
      my ($lex, $force) = @_;
      my $save = $lex->{'parse_state'}();
      print(STDERR 'In ', whoami(), ', left over: ',
            join(', ', $lex->{'rest'}()), "\n")
          if $debug > 1;
      my @postfix_expression;
      my $op;
      my @unary_expression;
      if (defined($op = $lex->{'currentsym'}())
          and (grep { $_ eq $op } @unaryOps)
          and $lex->{'advancesym'}())
      {
         print(STDERR 'In ', whoami(), ', line ', __LINE__, ', left over: ',
               join(', ', $lex->{'rest'}()), "\n") if $debug > 1;
         my @cast_expression
             = cast_expression($lex, $force);
         if (@cast_expression)
         {
            @unary_expression = ('unary_expression', $op,
                                 \@cast_expression);
            return @unary_expression;
         }
      }
      elsif (@postfix_expression = postfix_expression($lex, $force))
      {
         return @postfix_expression;
      }

      return parse_backoff($lex, $save, $force,
                           'Invalid unary_expression ',
                           $lex->{'consumed'}($save));
   }

   my @postfixOps = ('[', '(', '.', '->', '++', '--');
   my %postfix_opdescs
       = ('[' => 'array subscript ([])',
          '(' => 'function call (())',
          '.' => 'member access (.)',
          '->' => 'indirect member access (->)',
          '++' => 'postfix increment (++)',
          '--' => 'postfix decrement (--)'
         );

   $evaluations{'postfix_expression'}
       = sub
   {
      die('Error: postfix expressions (type: ',
          ') are not allowed in CPP expressions.',
          "\n");
   };

   sub postfix_expression
   {
      local $_;
      my ($lex, $force) = @_;
      my $save = $lex->{'parse_state'}();
      print(STDERR 'In ', whoami(), ', left over: ',
            join(', ', $lex->{'rest'}()), "\n")
          if $debug > 1;
      my @primary_expression
          = primary_expression($lex, $force);
      if (@primary_expression
          and (!defined($lex->{'currentsym'}())
               or !grep { $_ eq $lex->{'currentsym'}() } @postfixOps))
      {
         return @primary_expression;
      }
      elsif (@primary_expression)
      {
         # we have a primary expression extensible to a postfix expression
         my (@postfix_expression) = @primary_expression;
         $save = $lex->{'parse_state'}();
         my $op;
         while (defined($op = $lex->{'currentsym'}())
                and (grep { $_ eq $op } @postfixOps)
                and $lex->{'advancesym'}())
         {
            print(STDERR 'In ', whoami(), ', left over: ',
                  join(', ', $lex->{'rest'}()), "\n") if $debug > 1;
            if ($op eq '[')
            {
               if (my @expression = expression($lex, $force)
                   and defined($op = $lex->{'currentsym'}())
                   and $op eq ']')
               {
                  @postfix_expression = ('postfix_expression', '[',
                                         \@postfix_expression, \@expression);
                  $save = $lex->{'parse_state'}();
               }
               else
               {
                  last;
               }
            }
            elsif ($op eq '(')
            {
               # argument expression list
               my @arguments;
               while (my @assignment_expression
                      = assignment_expression($lex, $force))
               {
                  my $token;
                  push(@arguments, \@assignment_expression);
                  if (defined($token = $lex->{'currentsym'}()))
                  {
                     last if ($token eq ')');
                     if ($token eq ',')
                     {
                        $lex->{'advancesym'}();
                        next;
                     }
                  }
                  else
                  {
                     last;
                  }
               }
               {
                  my $token;
                  if (defined($token = $lex->{'currentsym'}())
                      and $token eq ')'
                      and $lex->{'advancesym'}())
                  {
                     @postfix_expression = ('postfix_expression', '(',
                                            \@postfix_expression, \@arguments);
                     $save = $lex->{'parse_state'}();
                  }
                  else
                  {
                     last;
                  }
               }
            }
            elsif ($op eq '.')
            {
               if (my @identifier = identifier($lex, $force))
               {
                  @postfix_expression = ('postfix_expression', '.',
                                         \@postfix_expression, \@identifier);
                  $save = $lex->{'parse_state'}();
               }
               else
               {
                  last;
               }
            }
            elsif ($op eq '->')
            {
               if (my @identifier = identifier($lex, $force))
               {
                  @postfix_expression = ('postfix_expression', '->',
                                         \@postfix_expression, \@identifier);
                  $save = $lex->{'parse_state'}();
               }
               else
               {
                  last;
               }
            }
            elsif ($op eq '++')
            {
               @postfix_expression = ('postfix_expression', '++',
                                      \@postfix_expression);
               $save = $lex->{'parse_state'}();
            }
            elsif ($op eq '--')
            {
               @postfix_expression = ('postfix_expression', '--',
                                      \@postfix_expression);
               $save = $lex->{'parse_state'}();
            }
         }
         $lex->{'set_parse_state'}($save);
         return @postfix_expression;
      }
      return parse_backoff($lex, $save, $force,
                           'Invalid postfix_expression ',
                           $lex->{'consumed'}($save));
   }

   sub primary_expression
   {
      local $_;
      my ($lex, $force) = @_;
      my $save = $lex->{'parse_state'}();
      print(STDERR 'In ', whoami(), ', left over: ',
            join(', ', $lex->{'rest'}()), "\n")
          if $debug > 1;
      my (@identifier, @constant, @string);
      if (@identifier = identifier($lex, 0))
      {
         return @identifier;
      }
      elsif (@constant = constant($lex, 0))
      {
         return @constant;
      }
      elsif (@string = string($lex, 0))
      {
         return @string;
      }
      else
      {
         my $currentsym;
         if (defined($currentsym = $lex->{'currentsym'}())
             and ($currentsym eq '(')
             and $lex->{'advancesym'}())
         {
            print(STDERR 'In ', whoami(), ', left over: ',
                  join(', ', $lex->{'rest'}()), "\n") if $debug > 1;
            my @expression
                = expression($lex, 0);
            if (@expression
                and defined($currentsym = $lex->{'currentsym'}())
                and $currentsym eq ')')
            {
               $lex->{'advancesym'}();
               return @expression;
            }
         }
      }
      return parse_backoff($lex, $save, $force,
                           'Invalid primary_expression ',
                           $lex->{'consumed'}($save));
   }

   sub identifier
   {
      local $_;
      my ($lex, $force) = @_;
      print(STDERR 'In ', whoami(), ', left over: ',
            join(', ', $lex->{'rest'}()), "\n")
          if $debug > 1;
      my $currentsym;
      my $save = $lex->{'parse_state'}();
      if (isSymbol($currentsym = $lex->{'currentsym'}()))
      {
         $lex->{'advancesym'}();
         return ('identifier', $currentsym);
      }
      return parse_backoff($lex, $save, $force,
                           'Invalid identifier ',
                           $lex->{'consumed'}($save));
   }

   $evaluations{'constant'}
       = sub
   {
      return ('constant', @_);
   };

   sub constant
   {
      local $_;
      my ($lex, $force) = @_;
      my $save = $lex->{'parse_state'}();
      print(STDERR 'In ', whoami(), ', left over: ',
            join(', ', $lex->{'rest'}()), "\n")
          if $debug > 1;
      my $token = $lex->{'currentsym'}();
      my @const_def;
      if (@const_def = integer_constant($lex, 0))
      {
         return @const_def;
      }
      elsif (@const_def = float_constant($lex, 0))
      {
         return @const_def;
      }
      elsif (@const_def = character_constant($lex, 0))
      {
         return @const_def;
      }
      elsif (@const_def = enumeration_constant($lex, 0))
      {
         return @const_def;
      }
      else
      {
         return parse_backoff($lex, $save, $force,
                              'Invalid primary_expression ',
                              $lex->{'consumed'}($save));
      }
   }

   sub integer_constant
   {
      my ($lex, $force) = @_;
      my $save = $lex->{'parse_state'}();
      print(STDERR 'In ', whoami(), ', left over: ',
            join(', ', $lex->{'rest'}()), "\n")
          if $debug > 1;
      local $_ = $lex->{'currentsym'}();
      if (m{$integral_constant_match}xi)
      {
         my ($digits, $integer_suffix) = ($1, $2);
         my $inttype = 'int';
         print(STDERR 'In ', whoami(), ', digits: ',
               $digits, ", suffix: ",
               defined($integer_suffix)?$integer_suffix:'none', "\n")
             if $debug > 1;
         if (defined($integer_suffix))
         {
            if ($integer_suffix =~ m{(u)?(l){0,2}(u)?}xi
                and !(defined($1) && defined($3)))
            {
               print(STDERR 'In ', whoami(), ', line ', __LINE__,
                     ', digits: ', $digits, ", suffix: ", $integer_suffix, "\n")
                   if $debug > 1;
               my ($unsigned_flag, $length_modifier)
                   = (lc($1 || $3 || ''), lc($2 || ''));
               $inttype = ($unsigned_flag eq 'u'?'unsigned-':'')
                   . ($length_modifier eq 'll'?
                      'long_long-':($length_modifier eq 'l'?'long-':''))
                       . 'int';
            }
            else
            {
               return parse_backoff($lex, $save, $force,
                                    'Invalid integer_constant ',
                                    $lex->{'consumed'}($save));
            }
         }
         $lex->{'advancesym'}();
         print(STDERR 'In ', whoami(), ', line ', __LINE__,
               ', left over: ', join(', ', $lex->{'rest'}()), "\n")
                   if $debug > 1;
         return ('constant', $inttype, 0+$digits);
      }
      return parse_backoff($lex, $save, $force,
                           'Invalid integer_constant ',
                           $lex->{'consumed'}($save));
   }

   sub float_constant
   {
      my ($lex, $force) = @_;
      my $save = $lex->{'parse_state'}();
      print(STDERR 'In ', whoami(), ', left over: ',
            join(', ', $lex->{'rest'}()), "\n")
          if $debug > 1;
      local $_ = $lex->{'currentsym'}();
      my $float_token = $_;
      if (m{^(\d*\.\d+|\d+\.)(?:([e])(\d*))?([fl]?)$}xi
          or m{^(\d+)([e])(\d*)([fl]?)$}xi)
      {
         my ($fractional_constant, $exponent_indicator,
             $exponent, $float_suffix) = ($1, $2, $3, $4);
         my $tokens_left = $lex->{'advancesym'}();
         my ($exponent_token, $exponent_sign_token) = (undef, '');
         $float_token =~ s{[fl]$}{}xi if defined($float_suffix);
         if (defined($exponent_indicator) && !defined($exponent))
         {
            if (!defined($float_suffix)
                and ((($exponent_token = $lex->{'currentsym'}()) eq '-'
                      or $exponent_token eq '+')
                     ? (($exponent_sign_token = $exponent_token)
                        and $lex->{'advancesym'}()
                        and (($exponent_token = $lex->{'currentsym'}())
                             =~ m{\d+[fl]?}xi))
                     : $exponent_token =~ m{\d+}x))
            {
               if ($exponent_token =~ m{([fl])$}xi)
               {
                  $float_suffix = $1;
                  $exponent_token =~ s{[fl]$}{}xi;
               }
               $float_token .= $exponent_sign_token . $exponent_token;
            }
            else
            {
               return parse_backoff($lex, $save, $force,
                                    'Invalid float_constant ',
                                    $lex->{'consumed'}($save));
            }
         }
         my $float_type = 'double';
         if (defined($float_suffix))
         {
            $float_suffix = lc($float_suffix);
            if ($float_suffix eq 'l')
            {
               $float_type = 'long-double';
            }
            elsif ($float_suffix eq 'f')
            {
               $float_type = 'float';
            }
            else
            {
               die('Unexpected error in float constant recognition.');
            }
         }
         return ('constant', $float_type, 0.0+$float_token);
      }
      else
      {
         return parse_backoff($lex, $save, $force,
                              'Invalid float_constant ',
                              $lex->{'consumed'}($save));
      }
   }

   sub character_constant
   {
      my ($lex, $force) = @_;
      my $save = $lex->{'parse_state'}();
      print(STDERR 'In ', whoami(), ', left over: ',
            join(', ', $lex->{'rest'}()), "\n")
          if $debug > 1;
      local $_ = $lex->{'currentsym'}();
      if (m{^'(.|\\[ntvbrfa\\?'"]|\\[0-7]{1,3}|\\x[0-9a-fA-F]+)'$}x)
      {
         return ('constant', 'char', $1);
      }
      else
      {
         return parse_backoff($lex, $save, $force,
                              'Invalid character_constant ',
                              $lex->{'consumed'}($save));
      }
   }

   sub string
   {
      my ($lex, $force) = @_;
      my $save = $lex->{'parse_state'}();
      print(STDERR 'In ', whoami(), ', left over: ',
            join(', ', $lex->{'rest'}()), "\n")
          if $debug > 1;
      local $_ = $lex->{'currentsym'}();
      if (m{^"((?:.|\\[ntvbrfa\\?'"]|\\[0-7]{1,3}|\\x[0-9a-fA-F]+)+)"$}x)
      {
         return ('constant', 'string', $1);
      }
      else
      {
         return parse_backoff($lex, $save, $force,
                              'Invalid string constant ',
                              $lex->{'consumed'}($save));
      }
   }

   sub enumeration_constant
   {
      return isSymbol($_[0]);
   }

   sub isSymbol
   {
      return $_[0] =~ m{^[a-zA-Z_]\w+$};
   }

   sub cpp_evaluate
   {
      local $_;
      my @ast = @_;
      my $toptype = $ast[0];
      print(STDERR 'Evaluating ', Data::Dumper->Dump([\@ast], ['ast']))
          if $debug > 1;
      if (exists($evaluations{$toptype}))
      {
         return $evaluations{$toptype}(@ast[1..$#ast]);
      }
      else
      {
         die('Unknown expression type: ', $toptype);
      }
   }

   sub cpp_true
   {
      my @value = @_;
      print(STDERR 'Returning truth value of ',
            Data::Dumper->Dump([\@value], ['value']))
          if $debug > 1;
      if ($value[0] ne 'constant')
      {
         die('Cannot evaluate non-constant expressions!');
      }
      if ($value[1] =~ m{(?:unsigned)?(?:long_long-|long-)?int})
      {
         return $value[2] != 0;
      }
      die('Unhandled truth expression: ', join(':', @_));
   }

   sub cpp_arithmetic_conversion
   {
      my @val = @_;
      if ($val[0][1] eq 'long-double' or $val[1][1] eq 'long-double')
      {
         return 'long-double';
      }
      elsif ($val[0][1] eq 'double' or $val[1][1] eq 'double')
      {
         return 'double';
      }
      elsif ($val[0][1] eq 'float' or $val[1][1] eq 'float')
      {
         return 'float';
      }
      else
      {
         @val = ([integral_promotion(@{$val[0]})],
                 [integral_promotion(@{$val[1]})]);
         if ($val[0][1] eq 'unsigned-long_long-int'
             or $val[1][1] eq 'unsigned-long_long-int')
         {
            return 'unsigned-long_long-int';
         }
         if ($c_sizeof_long_long > $c_sizeof_long)
         {
            if ($val[0][1] eq 'long_long-int'
                or $val[1][1] eq 'long_long-int')
            {
               return 'long_long-int';
            }
         }
         else
         {
            # long long can't represent every value in unsigned long
            if (($val[0][1] eq 'long_long-int'
                 and $val[1][1] eq 'unsigned-long-int')
                or ($val[1][1] eq 'long_long-int'
                    and $val[0][1] eq 'unsigned-long-int'))
            {
               return 'unsigned-long_long-int';
            }
         }
         if ($val[0][1] eq 'unsigned-long_long-int'
             or $val[1][1] eq 'unsigned-long_long-int')
         {
            return 'unsigned-long_long-int';
         }
         if ($c_sizeof_long > $c_sizeof_int)
         {
            if ($val[0][1] eq 'long-int'
                or $val[1][1] eq 'long-int')
            {
               return 'long-int';
            }
         }
         else
         {
            # long can't represent every value in unsigned
            if (($val[0][1] eq 'long-int'
                 and $val[1][1] eq 'unsigned-int')
                or ($val[1][1] eq 'long-int'
                    and $val[0][1] eq 'unsigned-int'))
            {
               return 'unsigned-long-int';
            }
         }
         if ($val[0][1] eq 'long-int'
             or $val[1][1] eq 'long-int')
         {
            return 'long-int';
         }
         elsif ($val[0][1] eq 'unsigned-int'
                or $val[1][1] eq 'unsigned-int')
         {
            return 'unsigned-int';
         }
         else
         {
            return 'int';
         }
      }
   }

   sub integral_promotion
   {
      my @val = @_;
      if ($val[1] eq 'char')
      {
         return ($val[0], 'int', unpack($c_char_is_unsigned?'C':'c', $val[2]));
      }
      else
      {
         return @val;
      }
   }

  sub init(%)
  {
    my (%conf) = %{$_[0]};
    if (exists($conf{'FPP'}))
    {
      @preprocCmd = split /\s+/, $conf{'FPP'};
    }
    if (exists($conf{'FPPFLAGS'}))
    {
      push(@preprocCmd, @{$conf{'FPPFLAGS'}});
    }
    $c_sizeof_int = $conf{'c_sizeof_int'}
        if (exists($conf{'c_sizeof_int'}));
    $c_sizeof_long = $conf{'c_sizeof_long'}
        if (exists($conf{'c_sizeof_long'}));
    $c_sizeof_long_long = $conf{'c_sizeof_long_long'}
        if (exists($conf{'c_sizeof_long_long'}));
    $c_char_is_unsigned = $conf{'c_char_is_unsigned'}
        if (exists($conf{'c_char_is_unsigned'}));
    print STDERR Data::Dumper->Dump([\%conf],['cond']),
        Data::Dumper->Dump([\@preprocCmd],['preprocCmd'])
            if ($debug);
  }

}


sub whoami  { (caller(1))[3] }

1;

# Local Variables:
# mode: cperl
# cperl-indent-level: 2
# license-project-url: "https://swprojects.dkrz.de/redmine/projects/scales-ppm"
# license-default: "bsd"
# End:

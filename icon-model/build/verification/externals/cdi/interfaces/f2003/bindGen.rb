#!/usr/bin/env ruby
$VERBOSE=true

# This script generates a fortran source file that uses the ISO_C_BINDINGS to interface to the functions defined in the given C header file.
# The basic approach is, that every C function is wrapped in a fortran function/subroutine, which internally uses a bind(c) interface to the C code.
# This wrapper based approach has the advantage that the wrapper is free to provide a true fortran interface
# that enables full type checking of its arguments; the pure bind(c) interface would not be able to distinguish
# between different opaque pointer types, for instance, nor would it be able to infer the size of a static string returned by a C function.
#
# Within this header file, the following constructs are recognized:
#
#   * #define FOO 123
#   * typedef struct foo foo;
#   * typedef struct foo { ... } foo;
#   * ... foo(...);
#
# These constructs are used to divide a source line into parts that are recognizable by the templates defined below.
# A function definition, for instance, is divided into a return type, a function name, and a number of argument definitions,
# the return type and argument descriptions are matched against templates which define the translation of these parts into fortran code.
# Note that all these constructs must be one-liners since processing in this script is line based.
#
# Every template is a hash that contains an entry :regex, which is used to match it against the corresponding C declaration.
# There are a couple of placeholders that may be used within these regex strings, they are expanded by matchTemplate() before a Regexp object is constructed from the string in :regex.
# These placeholders are:
#	<integerTypes>	matches the C integer types that can be used within Fortran by prefixing 'c_' to the type
#	<floatTypes>	matches the C floating point types that can be used within Fortran by prefixing 'c_' to the type
#	<opaqueTypes>	matches all the opaque types defined within the header
#	<publicTypes>	matches all the public types defined within the header
#
# In the case of argument and type templates, this :regex may contain one or more named subexpressions /(?<name>...)/,
# which can be included in the other fields by means of a corresponding placeholder "<name>".
# The names of the subexpressions that are to be substituted in this way need to be listed in the :placeholders key.
# This is usually used to capture the variable name, and then use "<name>_foo" to derive fortran variable names from the argument name,
# but it may also be used to capture the size of an array declaration.
# Since fortran uses so many keywords that can easily conflict with C argument names, it is a good idea not to use a naked "<name>";
# always append something to it as in "<name>_dummy"
#
# Argument templates must provide the following fields:
#	:regex	A regex that matches the whole definition of a C argument. Make sure it only matches the cases that the template can actually handle!
#	:placeholders	An array of the name of the named subexpressions used in the regex. For the :regex => /(?<foo>.),(?<bar>.)/ you would use :placeholders => %w[foo bar]
#	:dummyName	The name of the fortran dummy argument. Both the wrapper function and the `bind(c)` interface use the same name.
#	:acceptAs	The declaration of the dummy argument in the fortran wrapper.
#	:helperVars	Declarations of additional variables needed to provide the desired functionality in the wrapper function.
#	:precallStatements	Code that needs to be executed before the C function is called.
#	:callExpression	The actual argument that the wrapper passes to the C function.
#	:passAs	The declaration of the dummy argument in the `bind(c)` interface.
#	:postcallStatements	Code that needs to be executed after the C function returns.
#
#
#
# Return type templates are similar to argument templates, but they have to deal with the fact that fortran differentiates between subroutines and functions. Because of this, return type templates add the :isVoid key which is only true if the C function returns `void`.
#
# Return type templates must provide the following fields:
#	:regex	A regex that matches the whole definition of a C return type. Make sure it only matches the cases that the template can actually handle!
#	:isVoid	Always false, except for the template for `void`.
#	:returnAs	The type of the fortran wrapper function.
#	:helperVars	Declarations of additional variables needed to provide the desired functionality in the wrapper function.
#	:precallStatements	Code that needs to be executed before the C function is called.
#	:receiveAs	The type of the `bind(c)` interface function.
#	:assignVariable	The expression that the result of the C function is assigned to.
#	:postcallStatements	Code that needs to be executed after the C function returns.
#
#
#
# Type templates are used for the variables in public `struct` definitions. These are much simpler as they only have to translate a C variable declaration into an interoperable fortran variable declaration.
#
# Type templates must provide the following fields:
#	:regex	A regex that matches the whole C variable definition. Make sure it only matches the cases that the template can actually handle!
#	:placeholders	An array of the name of the named subexpressions used in the regex. Same semantics as in an argument template.
#	:declareAs	The declaration of the corresponding fortran derived type member.
#
#
#
# The wrapper that is generated for a non-void C function looks like this:
#
#	function fname(:dummyName...) result(result)
#		:returnAs :: result
#		:acceptAs...
#		:helperVars...
#		interface
#			:receiveAs function lib_fname(:dummyName...) bind(c, name = 'fname')
#				import <importConstants>
#				:passAs...
#			end function lib_fname
#		end interface
#		:precallStatements
#		:assignVariable = lib_fname(:callExpression)
#		:postcallStatements
#	end function fname
#
#
#
# The wrapper that is generated for a void C function looks like this:
#
#	subroutine fname(:dummyName...)
#		:acceptAs...
#		:helperVars...
#		interface
#			subroutine lib_fname(:dummyName...) bind(c, name = 'fname')
#				import <importConstants>
#				:passAs...
#			end subroutine lib_fname
#		end interface
#		:precallStatements
#		call lib_fname(:callExpression)
#		:postcallStatements
#	end subroutine fname

$debug = 0
$wrapperResultVarName = 'f_result'

####################################################################################################
# Template definitions #############################################################################
####################################################################################################

$argumentTemplates =
	[
	 {	#Dummy for declarations using foo(void).
		 :regex => '^\s*void\s*$',
		 :placeholders => %w[],
		 :dummyName => '',
		 :acceptAs => '',
		 :helperVars => '',
		 :precallStatements => '',
		 :callExpression => '',
		 :passAs => '',
		 :postcallStatements => ''
	 }, {	#<integerTypes>
		 :regex => '^\s*(?<type><integerTypes>)\s+(?<name>\w+)\s*$',
		 :placeholders => %w[name type],
		 :dummyName => '<name>_dummy',
		 :acceptAs => 'integer(c_<type>), value :: <name>_dummy',
		 :helperVars => '',
		 :precallStatements => '',
		 :callExpression => '<name>_dummy',
		 :passAs => 'integer(c_<type>), value :: <name>_dummy',
		 :postcallStatements => ''
	 }, {	#logical
		 :regex => '^\s*(?<type>bool)\s+(?<name>\w+)\s*$',
		 :placeholders => %w[name type],
		 :dummyName => '<name>_dummy',
		 :acceptAs => 'logical(c_<type>), value :: <name>_dummy',
		 :helperVars => '',
		 :precallStatements => '',
		 :callExpression => '<name>_dummy',
		 :passAs => 'logical(c_<type>), value :: <name>_dummy',
		 :postcallStatements => ''
	 }, {	#<floatTypes>
		 :regex => '^\s*(?<type><floatTypes>)\s+(?<name>\w+)\s*$',
		 :placeholders => %w[name type],
		 :dummyName => '<name>_dummy',
		 :acceptAs => 'real(c_<type>), value :: <name>_dummy',
		 :helperVars => '',
		 :precallStatements => '',
		 :callExpression => '<name>_dummy',
		 :passAs => 'real(c_<type>), value :: <name>_dummy',
		 :postcallStatements => ''
	 },
	 #Optional pointer arguments. These are wrapped by optional
	 #arguments which have to be named in calling code, which is why we
	 #don't use the _dummy suffix for them.
	 {	#<integerTypes>*
		 :regex => '^\s*(?<type><integerTypes>)\s*\*\s*(?<name>\w+)_optional\s*$',
		 :placeholders => %w[name type],
		 :dummyName => '<name>',
		 :acceptAs => 'integer(c_<type>), optional, target, intent(inout) :: <name>',
		 :helperVars => "type(c_ptr) :: <name>_cptr",
		 :precallStatements => "<name>_cptr = c_null_ptr\nif(present(<name>)) <name>_cptr = c_loc(<name>)",
		 :callExpression => '<name>_cptr',
		 :passAs => 'type(c_ptr), value :: <name>',
		 :postcallStatements => ""
	 },
	 {	#<floatTypes>*
		 :regex => '^\s*(?<type><floatTypes>)\s*\*\s*(?<name>\w+)_optional\s*$',
		 :placeholders => %w[name type],
		 :dummyName => '<name>',
		 :acceptAs => 'real(c_<type>), target, optional, intent(inout) :: <name>',
		 :helperVars => "type(c_ptr) :: <name>_cptr",
		 :precallStatements => "<name>_cptr = c_null_ptr\nif(present(<name>)) <name>_cptr = c_loc(<name>)",
		 :callExpression => '<name>_cptr',
		 :passAs => 'type(c_ptr), value :: <name>',
		 :postcallStatements => ""
	 },
	 {	#unsigned char (*<name>)[<size>]
		 :regex => '^\s*unsigned\s+char\s*\(\s*\*\s*(?<name>\w+)_optional\s*\)\s*\[\s*(?<size>[^\]]+)\s*\]\s*$',
		 :placeholders => %w[name size],
		 :dummyName => '<name>',
		 :acceptAs => 'integer(kind = c_signed_char), target, optional, intent(inout) :: <name>(<size>)',
		 :helperVars => "type(c_ptr) :: <name>_cptr",
		 :precallStatements => "<name>_cptr = c_null_ptr\nif(present(<name>)) <name>_cptr = c_loc(<name>)",
		 :callExpression => '<name>_cptr',
		 :passAs => 'type(c_ptr), value :: <name>',
		 :postcallStatements => ""
	 },
	 {	#<integerTypes> []
		 :regex => '^\s*(?<type><integerTypes>)\s*(?<name>\w+)_optional\[\]\s*$',
		 :placeholders => %w[name type],
		 :dummyName => '<name>',
		 :acceptAs => 'integer(c_<type>), target, optional, intent(inout) :: <name>(*)',
		 :helperVars => "type(c_ptr) :: <name>_cptr",
		 :precallStatements => "<name>_cptr = c_null_ptr\nif(present(<name>)) <name>_cptr = c_loc(<name>)",
		 :callExpression => '<name>_cptr',
		 :passAs => 'type(c_ptr), value :: <name>',
		 :postcallStatements => ""
	 }, {	#<floatTypes> []
		 :regex => '^\s*(?<type><floatTypes>)\s*(?<name>\w+)_optional\[\]\s*$',
		 :placeholders => %w[name type],
		 :dummyName => '<name>',
		 :acceptAs => 'real(c_<type>), target, optional, intent(inout) :: <name>(*)',
		 :helperVars => "type(c_ptr) :: <name>_cptr",
		 :precallStatements => "<name>_cptr = c_null_ptr\nif(present(<name>)) <name>_cptr = c_loc(<name>)",
		 :callExpression => '<name>_cptr',
		 :passAs => 'type(c_ptr), value :: <name>',
		 :postcallStatements => ""
	 }, {	#unsigned char <name>[<size>]
		 :regex => '^\s*unsigned\s+char\s*(?<name>\w+)_optional\s*\[\s*(?<size>[^\]]+)\s*\]\s*$',
		 :placeholders => %w[name size],
		 :dummyName => '<name>',
		 :acceptAs => 'integer(kind = c_signed_char), target, optional, intent(inout) :: <name>(<size>)',
		 :helperVars => "type(c_ptr) :: <name>_cptr",
		 :precallStatements => "<name>_cptr = c_null_ptr\nif(present(<name>)) <name>_cptr = c_loc(<name>)",
		 :callExpression => '<name>_cptr',
		 :passAs => 'type(c_ptr), value :: <name>',
		 :postcallStatements => ""
	 },
	 {	#optional char *
		 :regex => '^\s*char\s*(?:\*\s*(?<name>\w+)_optional|(?<name>\w+)_optional\[\])\s*$',
		 :placeholders => %w[name],
		 :dummyName => '<name>',
		 :acceptAs => 'character(kind = c_char, len = *), optional, intent(inout) :: <name>',
		 :helperVars => "character(kind = c_char), allocatable, target :: <name>_temp(:)\n" +
		 "integer :: <name>_i\n" +
		 "type(c_ptr) :: <name>_cptr",
		 :precallStatements => "if (present(<name>)) then\n" +
		 "\tallocate(<name>_temp(len(<name>) + 1))\n" +
		 "\t<name>_temp(len(<name>) + 1) = c_null_char\n" +
		 "\tdo <name>_i = len(<name>), 1, -1\n" +
		 "\t\tif(<name>(<name>_i:<name>_i) /= ' ') exit\n" +
		 "\t\t<name>_temp(<name>_i) = c_null_char\n" +
		 "\tend do\n" +
		 "\tdo <name>_i = <name>_i, 1, -1\n" +
		 "\t\t\t<name>_temp(<name>_i) = <name>(<name>_i:<name>_i)\n" +
		 "\tend do\n" +
		 "\t<name>_cptr = c_loc(<name>_temp)\n" +
		 "else\n" +
		 "\t<name>_cptr = c_null_ptr\n" +
		 "end if",
		 :callExpression => '<name>_cptr',
		 :passAs => 'type(c_ptr), value :: <name>',
		 :postcallStatements => "if (present(<name>)) then\n" +
		 "\tdo <name>_i = 1, len(<name>)\n" +
		 "\t\tif(<name>_temp(<name>_i) == c_null_char) exit\n" +
		 "\t\t<name>(<name>_i:<name>_i) = <name>_temp(<name>_i)\n" +
		 "\tend do\n" +
		 "\tdo <name>_i = <name>_i, len(<name>)\n" +
		 "\t\t<name>(<name>_i:<name>_i) = ' '\n" +
		 "\tend do\n" +
		 "end if"
	 },
	 {	#optional const char *
		 :regex => '^\s*const\s+char\s*(?:\*\s*(?<name>\w+)_optional|(?<name>\w+)_optional\[\])\s*$',
		 :placeholders => %w[name],
		 :dummyName => '<name>',
		 :acceptAs => 'character(kind = c_char, len = *), optional, intent(in) :: <name>',
		 :helperVars => "character(kind = c_char), allocatable, target :: <name>_temp(:)\n" +
		 "integer :: <name>_i\n" +
		 "type(c_ptr) :: <name>_cptr",
		 :precallStatements => "if (present(<name>)) then\n" +
		 "\tallocate(<name>_temp(len(<name>) + 1))\n" +
		 "\t<name>_temp(len(<name>) + 1) = c_null_char\n" +
		 "\tdo <name>_i = len(<name>), 1, -1\n" +
		 "\t\tif(<name>(<name>_i:<name>_i) /= ' ') exit\n" +
		 "\t\t<name>_temp(<name>_i) = c_null_char\n" +
		 "\tend do\n" +
		 "\tdo <name>_i = <name>_i, 1, -1\n" +
		 "\t\t\t<name>_temp(<name>_i) = <name>(<name>_i:<name>_i)\n" +
		 "\tend do\n" +
		 "\t<name>_cptr = c_loc(<name>_temp)\n" +
		 "else\n" +
		 "\t<name>_cptr = c_null_ptr\n" +
		 "end if",
		 :callExpression => '<name>_cptr',
		 :passAs => 'type(c_ptr), value :: <name>',
		 :postcallStatements => ""
	 },
	 #Array arguments. These are marked by a `_vec` suffix by convention.
	 #Since it's near impossible to write regexs that only match names that do *not* end in a given suffix,
	 #these templates must precede the more general templates for pointer arguments.
	 #That way, we can override the more general template with the more special one if both match.
	 {	#<integerTypes>* <name>_vec
		 :regex => '^\s*(?<type><integerTypes>)\s*(?:\*\s*(?<name>\w+_vec)|(?<name>\w+)\[\])\s*$',
		 :placeholders => %w[name type],
		 :dummyName => '<name>_dummy',
		 :acceptAs => 'integer(c_<type>), intent(inout) :: <name>_dummy(*)',
		 :helperVars => "",
		 :precallStatements => "",
		 :callExpression => '<name>_dummy',
		 :passAs => 'integer(c_<type>), intent(inout) :: <name>_dummy(*)',
		 :postcallStatements => ""
	 }, {	#<floatTypes>* <name>_vec
		 :regex => '^\s*(?<type><floatTypes>)\s*(?:\*\s*(?<name>\w+_vec)|(?<name>\w+)\[\])\s*$',
		 :placeholders => %w[name type],
		 :dummyName => '<name>_dummy',
		 :acceptAs => 'real(c_<type>), intent(inout) :: <name>_dummy(*)',
		 :helperVars => "",
		 :precallStatements => "",
		 :callExpression => '<name>_dummy',
		 :passAs => 'real(c_<type>), intent(inout) :: <name>_dummy(*)',
		 :postcallStatements => ""
	 }, {	#unsigned char <name>[<size>]
		 :regex => '^\s*unsigned\s+char\s+(?<name>\w+)\s*\[\s*(?<size>[^\]]+)\s*\]\s*$',
		 :placeholders => %w[name size],
		 :dummyName => '<name>_dummy',
		 :acceptAs => 'integer(kind = c_signed_char), intent(inout) :: <name>_dummy(<size>)',
		 :helperVars => "",
		 :precallStatements => "",
		 :callExpression => '<name>_dummy',
		 :passAs => 'integer(kind = c_signed_char), intent(inout) :: <name>_dummy(<size>)',
		 :postcallStatements => ""
	 }, {	#unsigned char <name>[]
		 :regex => '^\s*unsigned\s+char\s+(?<name>\w+)\s*\[\]\s*$',
		 :placeholders => %w[name],
		 :dummyName => '<name>_dummy',
		 :acceptAs => 'integer(kind = c_signed_char), intent(inout) :: <name>_dummy(*)',
		 :helperVars => "",
		 :precallStatements => "",
		 :callExpression => '<name>_dummy',
		 :passAs => 'integer(kind = c_signed_char), intent(inout) :: <name>_dummy(*)',
		 :postcallStatements => ""
	 }, {	#const <integerTypes>* <name>_vec
		 :regex => '^\s*const\s+(?<type><integerTypes>)\s*(?:\*\s*(?<name>\w+_vec)|(?<name>\w+)\[\])\s*$',
		 :placeholders => %w[name type],
		 :dummyName => '<name>_dummy',
		 :acceptAs => 'integer(c_<type>), intent(in) :: <name>_dummy(*)',
		 :helperVars => "",
		 :precallStatements => "",
		 :callExpression => '<name>_dummy',
		 :passAs => 'integer(c_<type>), intent(in) :: <name>_dummy(*)',
		 :postcallStatements => ""
	 }, {	#const <floatTypes>* <name>_vec
		 :regex => '^\s*const\s+(?<type><floatTypes>)\s*(?:\*\s*(?<name>\w+_vec)|(?<name>\w+)\[\])\s*$',
		 :placeholders => %w[name type],
		 :dummyName => '<name>_dummy',
		 :acceptAs => 'real(c_<type>), intent(in) :: <name>_dummy(*)',
		 :helperVars => "",
		 :precallStatements => "",
		 :callExpression => '<name>_dummy',
		 :passAs => 'real(c_<type>), intent(in) :: <name>_dummy(*)',
		 :postcallStatements => ""
	 }, {	#const unsigned char <name>[]
		 :regex => '^\s*(const\s+unsigned\s+char|unsigned\s+char\s+const)\s+(?<name>\w+)\s*\[\]\s*$',
		 :placeholders => %w[name],
		 :dummyName => '<name>_dummy',
		 :acceptAs => 'integer(kind = c_signed_char), intent(in) :: <name>_dummy(*)',
		 :helperVars => "",
		 :precallStatements => "",
		 :callExpression => '<name>_dummy',
		 :passAs => 'integer(kind = c_signed_char), intent(in) :: <name>_dummy(*)',
		 :postcallStatements => ""
	 }, {	#const unsigned char <name>[<size>]
		 :regex => '^\s*(const\s+unsigned\s+char|unsigned\s+char\s+const)\s+(?<name>\w+)\s*\[\s*(?<size>[^\]]+)\s*\]\s*$',
		 :placeholders => %w[name size],
		 :dummyName => '<name>_dummy',
		 :acceptAs => 'integer(kind = c_signed_char), intent(in) :: <name>_dummy(<size>)',
		 :helperVars => "",
		 :precallStatements => "",
		 :callExpression => '<name>_dummy',
		 :passAs => 'integer(kind = c_signed_char), intent(in) :: <name>_dummy(<size>)',
		 :postcallStatements => ""
	 }, {	#const <integerTypes> <name>[<lineCount>][<lineSize>]
		 :regex => '^\s*const\s+(?<type><integerTypes>)\s+(?<name>\w+)\s*\[\s*(?<lineCount>[^\]]+)\s*\]\s*\[\s*(?<lineSize>[^\]]+)\s*\]\s*$',
		 :placeholders => %w[name type lineCount lineSize],
		 :dummyName => '<name>_dummy',
		 :acceptAs => 'integer(c_<type>), intent(in) :: <name>_dummy(<lineSize>, <lineCount>)',
		 :helperVars => "",
		 :precallStatements => "",
		 :callExpression => '<name>_dummy',
		 :passAs => 'integer(c_<type>), intent(in) :: <name>_dummy(<lineSize>, <lineCount>)',
		 :postcallStatements => ""
	 }, {	#const <integerTypes> <name>[][<lineSize>]
		 :regex => '^\s*const\s+(?<type><integerTypes>)\s+(?<name>\w+)\s*\[\s*\]\s*\[\s*(?<lineSize>[^\]]+)\s*\]\s*$',
		 :placeholders => %w[name type lineSize],
		 :dummyName => '<name>_dummy',
		 :acceptAs => 'integer(c_<type>), intent(in) :: <name>_dummy(<lineSize>, *)',
		 :helperVars => "",
		 :precallStatements => "",
		 :callExpression => '<name>_dummy',
		 :passAs => 'integer(c_<type>), intent(in) :: <name>_dummy(<lineSize>, *)',
		 :postcallStatements => ""
	 },
	 #Non-optional pointer arguments. These match both pointers and
	 #arrays, so they must appear after the more special array templates.
	 {	#<integerTypes>*
		 :regex => '^\s*(?<type><integerTypes>)\s*\*\s*(?<name>\w+)\s*$',
		 :placeholders => %w[name type],
		 :dummyName => '<name>_dummy',
		 :acceptAs => 'integer(c_<type>), intent(inout) :: <name>_dummy',
		 :helperVars => "",
		 :precallStatements => "",
		 :callExpression => '<name>_dummy',
		 :passAs => 'integer(c_<type>), intent(inout) :: <name>_dummy',
		 :postcallStatements => ""
	 }, {	#<floatTypes>*
		 :regex => '^\s*(?<type><floatTypes>)\s*\*\s*(?<name>\w+)\s*$',
		 :placeholders => %w[name type],
		 :dummyName => '<name>_dummy',
		 :acceptAs => 'real(c_<type>), intent(inout) :: <name>_dummy',
		 :helperVars => "",
		 :precallStatements => "",
		 :callExpression => '<name>_dummy',
		 :passAs => 'real(c_<type>), intent(inout) :: <name>_dummy',
		 :postcallStatements => ""
	 }, {	#unsigned char (*<name>)[<size>]
		 :regex => '^\s*unsigned\s+char\s*\(\s*\*\s*(?<name>\w+)\s*\)\s*\[\s*(?<size>[^\]]+)\s*\]\s*$',
		 :placeholders => %w[name size],
		 :dummyName => '<name>_dummy',
		 :acceptAs => 'integer(kind = c_signed_char), intent(inout) :: <name>(<size>)',
		 :helperVars => "",
		 :precallStatements => "",
		 :callExpression => '<name>_dummy',
		 :passAs => 'integer(kind = c_signed_char), intent(inout) :: <name>(<size>)',
		 :postcallStatements => ""
	 },
	 #String arguments.
	 {	#char*	Unsafe buffer passing
		 :regex => '^\s*char\s*\*\s*(?<name>\w+)\s*$',
		 :placeholders => %w[name],
		 :dummyName => '<name>_dummy',
		 :acceptAs => 'character(kind = c_char, len = *), intent(inout) :: <name>_dummy',
		 :helperVars => "character(kind = c_char) :: <name>_temp(len(<name>_dummy) + 1)\n" +
		 "integer :: <name>_i",
		 :precallStatements => "<name>_temp(len(<name>_dummy) + 1) = c_null_char\n" +
		 "do <name>_i = len(<name>_dummy), 1, -1\n" +
		 "\tif(<name>_dummy(<name>_i:<name>_i) /= ' ') exit\n" +
		 "\t<name>_temp(<name>_i) = c_null_char\n" +
		 "end do\n" +
		 "do <name>_i = <name>_i, 1, -1\n" +
		 "\t\t<name>_temp(<name>_i) = <name>_dummy(<name>_i:<name>_i)\n" +
		 "end do",
		 :callExpression => '<name>_temp',
		 :passAs => 'character(kind = c_char) :: <name>_dummy(*)',
		 :postcallStatements => "do <name>_i = 1, len(<name>_dummy)\n" +
		 "\tif(<name>_temp(<name>_i) == c_null_char) exit\n" +
		 "\t<name>_dummy(<name>_i:<name>_i) = <name>_temp(<name>_i)\n" +
		 "end do\n" +
		 "do <name>_i = <name>_i, len(<name>_dummy)\n" +
		 "\t<name>_dummy(<name>_i:<name>_i) = ' '\n" +
		 "end do"
	 }, {	#const char*	Safe passing of an input string.
		 :regex => '^\s*(const\s+char|char\s+const)\s*\*\s*(?<name>\w+)\s*$',
		 :placeholders => %w[name],
		 :dummyName => '<name>_dummy',
		 :acceptAs => 'character(kind = c_char, len = *), intent(in) :: <name>_dummy',
		 :helperVars => "character(kind = c_char) :: <name>_temp(len(<name>_dummy) + 1)\ninteger :: <name>_i",
		 :precallStatements => "do <name>_i = 1, len(<name>_dummy)\n\t<name>_temp(<name>_i) = <name>_dummy(<name>_i:<name>_i)\nend do\n<name>_temp(len(<name>_dummy) + 1) = c_null_char",
		 :callExpression => '<name>_temp',
		 :passAs => 'character(kind = c_char) :: <name>_dummy(*)',
		 :postcallStatements => ''
	 }, {	#char**	Safe returning of an output string.
		 :regex => '^\s*char\s*\*\s*\*\s*(?<name>\w+)_optional\s*$',
		 :placeholders => %w[name],
		 :dummyName => '<name>',
		 :acceptAs => 'character(kind = c_char), pointer, optional, intent(inout) :: <name>(:)',
		 :helperVars => "type(c_ptr), target :: <name>_cptr\n" +
		 "type(c_ptr) :: <name>_handle\n" +
		 "integer :: <name>_shape(1)\n" +
		 "character(kind = c_char), pointer :: <name>_fptr(:)",
		 :precallStatements => "<name>_handle = c_null_ptr\n" +
		 "if(present(<name>)) <name>_handle = c_loc(<name>_cptr)",
		 :callExpression => '<name>_handle',
		 :passAs => 'type(c_ptr), value :: <name>',
		 :postcallStatements => "if(present(<name>)) then\n" +
		 "\tif(c_associated(<name>_cptr)) then\n" +
		 "\t\t<name>_shape(1) = int(lib_strlen(<name>_cptr))\n" +
		 "\t\tcall c_f_pointer(<name>_cptr, <name>_fptr, <name>_shape)\n" +
		 "\t\tallocate(<name>(<name>_shape(1)))\n" +
		 "\t\t<name> = <name>_fptr\n" +
		 "\t\tcall lib_free(<name>_cptr)\n" +
		 "\telse\n" +
		 "\t\t<name> => null()\n" +
		 "\tend if\n" +
		 "end if"
	 },
	 #Public and opaque types
	 {	#[const] <opaqueTypes>*
		 :regex => '^\s*(const\s+|)(?<type><opaqueTypes>)(\s+const|)\s*\*\s*(?<name>\w+)\s*$',
		 :placeholders => %w[name type],
		 :dummyName => '<name>_dummy',
		 :acceptAs => 'type(t_<type>), intent(in) :: <name>_dummy',
		 :helperVars => '',
		 :precallStatements => '',
		 :callExpression => '<name>_dummy%ptr',
		 :passAs => 'type(c_ptr), value :: <name>_dummy',
		 :postcallStatements => ''
	 }
	]

$returnTypeTemplates = [
	{	#void
		:regex => '^\s*void\s*$',
		:placeholders => %w[],
		:isVoid => true
	}, {	#<integerTypes>
		:regex => '^\s*(?<type><integerTypes>)\s*$',
		:placeholders => %w[type],
		:isVoid => false,
		:returnAs => 'integer(c_<type>)',
		:helperVars => '',
		:precallStatements => '',
		:receiveAs => 'integer(c_<type>)',
		:assignVariable => $wrapperResultVarName,
		:postcallStatements => ''
	}, {	#bool
		:regex => '^\s*(?<type>bool)\s*$',
		:placeholders => %w[type],
		:isVoid => false,
		:returnAs => 'logical(c_<type>)',
		:helperVars => '',
		:precallStatements => '',
		:receiveAs => 'logical(c_<type>)',
		:assignVariable => $wrapperResultVarName,
		:postcallStatements => ''
	}, {	#<floatTypes>
		:regex => '^\s*(?<type><floatTypes>)\s*$',
		:placeholders => %w[type],
		:isVoid => false,
		:returnAs => 'real(c_<type>)',
		:helperVars => '',
		:precallStatements => '',
		:receiveAs => 'real(c_<type>)',
		:assignVariable => $wrapperResultVarName,
		:postcallStatements => ''
	}, {	#char*
		:regex => '^\s*char\s*\*\s*$',
		:placeholders => %w[],
		:isVoid => false,
		:returnAs => 'character(kind = c_char), dimension(:), pointer',
		:helperVars => "type(c_ptr) :: cString\n" +
		               "integer :: rv_shape(1)\n" +
		               "character(kind = c_char), dimension(:), pointer :: temp",
		:precallStatements => '',
		:receiveAs => 'type(c_ptr)',
		:assignVariable => 'cString',
		:postcallStatements => "if(c_associated(cString)) then\n" +
		                       "\trv_shape(1) = int(lib_strlen(cString))\n" +
		                       "\tcall c_f_pointer(cString, temp, rv_shape)\n" +
		                       "\tallocate(#{$wrapperResultVarName}(rv_shape(1)))\n" +
		                       "\t#{$wrapperResultVarName} = temp\n" +
		                       "\tcall lib_free(cString)\n" +
		                       "else\n" +
		                       "\t#{$wrapperResultVarName} => null()\n" +
		                       "end if"
	}, {	#const char*
		:regex => '^\s*const\s+char\s*\*\s*$',
		:placeholders => %w[],
		:isVoid => false,
		:returnAs => 'character(kind = c_char), dimension(:), pointer',
		:helperVars => "type(c_ptr) :: ptr\ninteger :: rv_shape(1)",
		:precallStatements => $wrapperResultVarName + ' => null()',
		:receiveAs => 'type(c_ptr)',
		:assignVariable => 'ptr',
		:postcallStatements => "if(c_associated(ptr)) then\n" +
		                       "\trv_shape(1) = int(lib_strlen(ptr))\n" +
		                       "\tcall c_f_pointer(ptr, #{$wrapperResultVarName}, rv_shape)\n" +
		                       "end if"
	}, {	#const int*	This returns the naked pointer because we can't know the length of the returned array within the wrapper. The user has to call c_f_pointer() himself.
		:regex => '^\s*const\s+(?<type><integerTypes>)\s*\*\s*$',
		:placeholders => %w[type],
		:isVoid => false,
		:returnAs => 'type(c_ptr)',
		:helperVars => '',
		:precallStatements => '',
		:receiveAs => 'type(c_ptr)',
		:assignVariable => $wrapperResultVarName,
		:postcallStatements => ''
	}, {	#const double*	This returns the naked pointer because we can't know the length of the returned array within the wrapper. The user has to call c_f_pointer() himself.
		:regex => '^\s*const\s+(?<type><floatTypes>)\s*\*\s*$',
		:placeholders => %w[type],
		:isVoid => false,
		:returnAs => 'type(c_ptr)',
		:helperVars => '',
		:precallStatements => '',
		:receiveAs => 'type(c_ptr)',
		:assignVariable => $wrapperResultVarName,
		:postcallStatements => ''
	},
	#Public and opaque types.
	{	#<publicTypes>
		:regex => '^\s*(?<type><publicTypes>)\s+$',
		:placeholders => %w[type],
		:isVoid => false,
		:returnAs => 'type(t_<type>)',
		:helperVars => '',
		:precallStatements => '',
		:receiveAs => 'type(t_<type>)',
		:assignVariable => $wrapperResultVarName,
		:postcallStatements => ''
	}, {	#<opaqueTypes>*
		:regex => '^\s*(?<type><opaqueTypes>)\s*\*\s*$',
		:placeholders => %w[type],
		:isVoid => false,
		:returnAs => 'type(t_<type>)',
		:helperVars => '',
		:precallStatements => '',
		:receiveAs => 'type(c_ptr)',
		:assignVariable => $wrapperResultVarName + '%ptr',
		:postcallStatements => ''
	}
]

$typeTemplates = [
	{	#<integerTypes>
		:regex => '^\s*(?<type><integerTypes>)\s+(?<name>\w+)\s*;$',
		:placeholders => %w[name type],
		:declareAs => "integer(c_<type>) :: <name>"
	}, {	#<floatTypes>
		:regex => '^\s*(?<type><floatTypes>)\s+(?<name>\w+)\s*;$',
		:placeholders => %w[name type],
		:declareAs => "real(c_<type>) :: <name>"
	}
]

####################################################################################################
# Verbatim Fortran Code ############################################################################
####################################################################################################

$verbatimDeclarations = '
	public ctrim
	public c_len

	interface
		integer(c_size_t) function lib_strlen(charPtr) bind(c, name = "strlen")
			import c_size_t, c_ptr
			type(c_ptr), value :: charPtr
		end function lib_strlen

		subroutine lib_free(ptr) bind(c, name = "free")
			import c_ptr
			type(c_ptr), value, intent(in) :: ptr
		end subroutine lib_free
	end interface
'

$verbatimDefinitions = "
	subroutine ctrim(str)
		character(kind = c_char, len = *), intent(inout) :: str
		integer :: i

		do i=1,len(str)
			if (str(i:i) == c_null_char) then
				str(i:len(str)) = ' '
				exit
			end if
		end do
	end subroutine ctrim

	function c_len(s) result(i)
		character(kind = c_char, len = *), intent(in) :: s
		integer :: i

		do i = 1, len(s)
			if (s(i:i) == c_null_char) exit
		end do
		i = i - 1
	end function
"

####################################################################################################
# Code to interpret the templates ##################################################################
####################################################################################################

$declarationLines = []
$definitionLines = []
$interfaceLines = []
$opaqueTypes = []
$publicTypes = []

def rubyVersionOk()
       version = RUBY_VERSION.split(".")
       if version[0].to_i > 1
               return true
       elsif version[0].to_i == 1
               return version[1].to_i >= 9
       else
               return false
       end
end

#This substitutes the placeholders <opaqueTypes> and <publicTypes> in the regexString prior to constructing a Regexp out of it.
def matchTemplate(regexString, matchString)
	opaqueTypesString = '(' + $opaqueTypes.collect{ |type| type	}.join('|') + ')'
	regexString = regexString.gsub("<opaqueTypes>", opaqueTypesString)
	publicTypesString = '(' + $publicTypes.collect{ |type| type }.join('|') + ')'
	regexString = regexString.gsub('<publicTypes>', publicTypesString)
	regexString = regexString.gsub('<integerTypes>', '\b(short|int|long|size_t|intmax_t|int_(least|fast)(8|16|32|64)_t|int(8|16|32|64)_t)\b')
	regexString = regexString.gsub('<floatTypes>',  '(float|double)')
	return Regexp.new(regexString).match(matchString)
end

class TemplateInstanciation
	def initialize(argumentString, template)
		@template = template
		@matchData = matchTemplate(template[:regex], argumentString)
		@placeholders = []
		template[:placeholders].each { |placeholder|
			@placeholders.push({ :name => placeholder, :regex => Regexp.new("<#{placeholder}>") })
		}
	end

	def expandTemplate(templateKey)
		result = @template[templateKey]
		#Replace all placeholders with their expansion.
		@placeholders.each { |current|
			result = result.gsub(current[:regex], @matchData[current[:name]])
		}
		return result
	end
	def nonEmptyKey?(templateKey)
		result = @template[templateKey] != ''
		return result
	end
	attr_reader :template
end

def formatLines(lineArray, indentation, string)
	$stderr.puts("Formatting '" + string + "'") if $debug > 3
	if string == "" && indentation == 0
		lineArray.push("")	#split() does not return anything if the string is empty, killing our empty lines
	end
	string.split("\n").each { |line|
		lineArray.push("\t"*indentation + line)
	}
	tail = ''
	if lineArray.length > 1
		tail = lineArray[-2..-1].join("\n")
	else
		tail = lineArray[-1]
	end
	$stderr.puts("Tail '" + tail + "'") if $debug > 3
end

def haveTemplateKey(templates, templateKey)
	result = false
	templates.count{ |template|
		result ||= template.nonEmptyKey?(templateKey)
		break if result
	}
	return result
end


def dumpStatements(indentation, argumentArray, templateKey, outputArray)
	argumentArray.each{ |argument|
		formatLines(outputArray, indentation, argument.expandTemplate(templateKey))
	}
end

def defineConstant(name, value)
	if /^(\+|-|)\d+$/.match(value)
		formatLines($declarationLines, 1, "integer(c_int), public, parameter :: #{name} = " + value)
	else
		$stderr.puts("Error: value '#{value}' of constant '#{name}' is not an integer literal")
	end
end

def defineOpaqueType(name)
	formatLines($declarationLines, 0, "")
	formatLines($declarationLines, 1, "public :: t_#{name}")
	formatLines($declarationLines, 1, "type :: t_#{name}")
	formatLines($declarationLines, 2, "type(c_ptr) :: ptr")
	formatLines($declarationLines, 1, "end type t_#{name}")
	$opaqueTypes.push(name)
end

def findTemplate(string, templateArray)
	templateArray.each do |template|
		if matchTemplate(template[:regex], string)
			return template
		end
	end
	return nil
end

def definePublicType(name, body)
	formatLines($declarationLines, 0, "")
	formatLines($declarationLines, 1, "public t_#{name}")
	formatLines($declarationLines, 1, "type, bind(c) :: t_#{name}")
	body.gsub(/[^;]+;/) do |variableDeclaration|
		if template = findTemplate(variableDeclaration, $typeTemplates)
			variable = TemplateInstanciation.new(variableDeclaration, template)
			formatLines($declarationLines, 2, variable.expandTemplate(:declareAs))
		else
			$stderr.puts("Error: Can't translate the declaration '" + variableDeclaration + "'")
		end
	end
	formatLines($declarationLines, 1, 'end type t_' + name)
	$publicTypes.push(name)
end

def collectImportConstants(importConstantsArray, typeString)
	$stderr.puts('Considering \'' + typeString + "' for import\n") if $debug > 2
	if importConstant = typeString[/\b[ct]_\w+\b/]
		importConstantsArray.push(importConstant)
	end
end

def collectArrayBoundImports(importConstantsArray, argPassString)
	if array_bounds_match =
			/(?<=\()\s*(?:[a-zA-Z_]\w*|:|\d+)\s*(?:,\s*(?:[a-zA-Z_]\w*|:|\d+)\s*?)*(?=\s*\)$)/.match(argPassString)
		$stderr.puts('Found array bounds: \'' +
								 array_bounds_match.to_a.join('\', \'') + "'",
								 'for \'' + argPassString + "'") if $debug > 3
		array_bounds_match = array_bounds_match[0]
		array_bounds_match.lstrip!
		array_bounds_match = array_bounds_match.split(/\s*,\s*/)
		array_bounds_match.select! do |bound|
			bound != ':' && bound !~ /^\d+$/
		end
		$stderr.puts('Possible import symbols in array bounds: \'' +
								 array_bounds_match.join('\', \'') + "'") if $debug > 3
		importConstantsArray.push(*array_bounds_match)
	end
end

#Collect the c_* and t_* constants/types from the arguments and the return type and build the corresponding `import` statement from them.
def importStatement(returnType, argumentArray)
	importConstants = []
	collectImportConstants(importConstants, returnType)
	argumentArray.each { |arg|
		argStr = arg.expandTemplate(:passAs)
		collectArrayBoundImports(importConstants, argStr)
		collectImportConstants(importConstants, argStr)
	}
	$stderr.puts('Returning \'' + importConstants.sort.uniq.join(', ') +
							 "' for import\n") if $debug > 2
	return (importConstants.length != 0) ? ('import ' + importConstants.sort.uniq.join(', ')) : ''
end

def defineFunction(name, arguments, returnType)
	#Find the relevant templates.
	if returnTemplate = findTemplate(returnType, $returnTypeTemplates)
		returnData = TemplateInstanciation.new(returnType, returnTemplate)
		argArray = []
		arguments.gsub(/[^,]+/) do |argument|
			if template = findTemplate(argument, $argumentTemplates)
				argArray.push(TemplateInstanciation.new(argument, template))
			else
				$stderr.puts("Error: type of argument '#{argument}' to function #{name}() is not supported")
				return
			end
		end
	else
		$stderr.puts("Error: Can't translate return type '#{returnType}' of function #{name}()")
		return
	end

	needArgPrologue = haveTemplateKey(argArray, :precallStatements)
	needRetValPrologue = !returnTemplate[:isVoid] &&
		returnData.nonEmptyKey?(:precallStatements)
	needArgEpilogue = haveTemplateKey(argArray, :postcallStatements)
	needRetValEpilogue = !returnTemplate[:isVoid] &&
		returnData.nonEmptyKey?(:postcallStatements)
	needRetValConversion = !returnTemplate[:isVoid] &&
		returnTemplate[:returnAs] != returnTemplate[:receiveAs]
	needArgConversion = false
	argArray.each do |arg|
		needArgConversion = (arg.template[:passAs] != arg.template[:acceptAs])
		break if needArgConversion
	end

	needWrapper = needArgPrologue || needRetValPrologue ||
		needArgEpilogue || needRetValEpilogue ||
		needRetValConversion || needArgConversion
	baseIndent = 0
	formatLines($declarationLines, 1, 'public :: ' + name)
	$stderr.puts(name + "\n" +
							 [ "\t"+'needArgPrologue => ' + needArgPrologue.to_s,
								 "\t"+'needRetValPrologue => ' + needRetValPrologue.to_s,
								 "\t"+'needArgEpilogue => ' + needArgEpilogue.to_s,
								 "\t"+'needRetValEpilogue => ' + needRetValEpilogue.to_s,
								 "\t"+'needRetValConversion => ' + needRetValConversion.to_s,
								 "\t"+'needArgConversion => ' + needArgConversion.to_s,
								 "\t"+'needWrapper => ' +
		needWrapper.to_s].join("\n")) if $debug > 1
	subprogramtype = returnTemplate[:isVoid] ? 'subroutine' : 'function'
	dummyArguments = argArray.collect{ |arg|
		arg.expandTemplate(:dummyName)
	}.join(", ")
	if needWrapper
		#Generate the wrapper function.
		line = subprogramtype + ' ' + name + '(' + dummyArguments + ')' +
			(returnTemplate[:isVoid] ? '' : " result(#{$wrapperResultVarName})")
		formatLines($definitionLines, 1, line)
		if !returnTemplate[:isVoid]
			formatLines($definitionLines, 2,
									returnData.expandTemplate(:returnAs) + ' :: ' +
									$wrapperResultVarName)
		end
		dumpStatements(               2, argArray, :acceptAs, $definitionLines)
		dumpStatements(               2, argArray, :helperVars, $definitionLines)
		if !returnTemplate[:isVoid]
			formatLines($definitionLines, 2, returnData.expandTemplate(:helperVars))
		end
		formatLines($definitionLines, 2, 'interface')
		line=subprogramtype + ' lib_' + name + '(' + dummyArguments + ') ' +
			'bind(c, name = \'' + name + '\')' +
			(returnTemplate[:isVoid] ? '' : ' result(c_result)')
		formatLines($definitionLines, 3, line)
		formatLines($definitionLines, 4,
								importStatement(returnTemplate[:isVoid] ? '' :
																returnData.expandTemplate(:receiveAs), argArray))
		if !returnTemplate[:isVoid]
			formatLines($definitionLines, 4,
									returnData.expandTemplate(:receiveAs) + ' :: c_result')
		end
		dumpStatements(               4, argArray, :passAs, $definitionLines)
		formatLines($definitionLines, 3, 'end ' + subprogramtype + ' lib_' + name)
		formatLines($definitionLines, 2, 'end interface')
		if needArgPrologue
			dumpStatements(               2, argArray, :precallStatements, $definitionLines)
		end
		if !returnTemplate[:isVoid]
			formatLines($definitionLines, 2,
									returnData.expandTemplate(:precallStatements))
		end
		formatLines($definitionLines, 2,
								(returnTemplate[:isVoid] ? 'call' :
								 (returnData.expandTemplate(:assignVariable) + ' =')) +
								' lib_' + name + '(' + argArray.collect{ |arg|
									arg.expandTemplate(:callExpression)
								}.join(', ') + ')')
		if (needArgEpilogue)
			dumpStatements(               2, argArray, :postcallStatements, $definitionLines)
		end
		if needRetValEpilogue
			formatLines($definitionLines, 2,
									returnData.expandTemplate(:postcallStatements))
		end
		formatLines($definitionLines, 1, 'end ' + subprogramtype + ' ' + name)
		formatLines($definitionLines, 0, '')
	else
		#Generate interface block only
		baseIndent = 1
		line = subprogramtype + ' ' + name + '(' + dummyArguments +
			') bind(c, name = \'' + name + '\')' +
			(returnTemplate[:isVoid] ? '' : " result(#{$wrapperResultVarName})")
		formatLines($interfaceLines, baseIndent + 1, line)
		line = importStatement(returnTemplate[:isVoid] ? '' :
													 returnData.expandTemplate(:returnAs),
													 argArray)
		$stderr.puts('Emitting \'' + line + "'") if $debug > 2
		formatLines($interfaceLines, baseIndent + 2, line)

		dumpStatements(baseIndent + 2, argArray, :passAs, $interfaceLines)
		if !returnTemplate[:isVoid]
			formatLines($interfaceLines, baseIndent + 2,
									returnData.expandTemplate(:returnAs) + ' :: ' +
									$wrapperResultVarName)
		end
		formatLines($interfaceLines, baseIndent + 1, 'end ' + subprogramtype + ' ' + name)
		formatLines($interfaceLines, 0, '')
	end
end

#Scan the given header and collect the interface information in the global variables.
def scanHeader(headerPath)
	#Scan the given header.
	headerFile = File.open(headerPath, 'r')
	commentTerm = %r{\*/}
	commentInit = %r{/[*/]}
	while line = headerFile.gets
		line.chomp!

		line.gsub!(/SizeType/,'int')
		line.gsub!(/DateType/,'int')

		commentBeginLineNo = headerFile.lineno
		while foundComment = commentInit.match(line)
			commentBegin = foundComment.begin(0)
			typeChar = line[commentBegin+1]
			if typeChar == "/"
				line = line[0,commentBegin]
				line.rstrip!
				foundComment = false
			elsif typeChar == "*"
				commentSubStr=line[commentBegin..-1]
				if commentTerm.match(commentSubStr)
					newLine =	line.gsub(%r{/\*.*?\*/},'')
					line = newLine.rstrip
				else
					while newLine = headerFile.gets and not foundCommentEnd = commentTerm.match(newLine)
					end
					if not newLine
						$stderr.puts("Unterminated comment started at line ",
												 commentBeginLineNo, "'", line, "'")
						exit(1)
					end
					newLine.chomp!
					line = line[0,commentBegin] +
								 newLine[foundCommentEnd.begin(0)+2..-1]
				end
			end
		end
		if /^\s*$/.match(line)
			#Empty lines are ignored.

		#Preprocessor stuff
		elsif matchedLine = /^\s*#\s*define\s+(?<symbol>\w+)\s+(?<value>.+)$/.match(line)
			defineConstant(matchedLine['symbol'], matchedLine['value'])
		elsif /^\s*#/.match(line)
			#All other preprocessor directives are ignored.

		#User defined types
		elsif matchedLine = /^\s*typedef\s+struct\s+(?<typeName>\w+)\s+\k<typeName>\s*;\s*$/.match(line)
			defineOpaqueType(matchedLine['typeName'])
		elsif matchedLine = /^\s*typedef\s+struct\s+(?<typeName>\w+)\s*{(?<body>.*)}\s*\k<typeName>\s*;\s*$/.match(line)
			definePublicType(matchedLine['typeName'], matchedLine['body'])

		#Function declarations
		elsif matchedLine = /^\s*(?<returnType>[^()]+)\b(?<functionName>\w+)\s*\((?<arguments>.*)\)\s*;\s*$/.match(line)
			defineFunction(matchedLine['functionName'], matchedLine['arguments'], matchedLine['returnType'])

		else
			STDERR.puts("Warning: Unrecognized line '#{line}'")
		end
	end
end

#Prints the line if it does not consist only of indentation, adding continuation lines as necessary.
def fortranLine(file, line)
	unless /^\t+$/.match(line)	#Intentionally empty lines don't contain indentation, so we preserve totally empty lines while throwing away the ones with leading tabs.
		# replace tabs with spaces first
		spacesPerTab = 2
		mline = line.gsub("\t", ' '*spacesPerTab)
		charsPerLine=79
		indentation = /^ */.match(mline)[0]
		if not %r{^ *!}.match(mline)
			while mline.length > charsPerLine
				# last position of space preceding line break
				tspos = mline[0..charsPerLine].rindex(' ') || charsPerLine
				file.puts(mline[0...tspos] + "&")
				mline = indentation + "&" + mline[tspos...mline.length]
			end
		end
		file.puts(mline)
	end
end

#Output the interface information in the global variables to a fortran file.
def writeFortranModule(scriptPath, headerPath, modulePath, moduleName)
	file = File.new(modulePath, "w")
	fortranLine(file, "! >>> Warning: This is a generated file. If you modify it, you get what you deserve. <<<")
	fortranLine(file, "!")
	fortranLine(file, "! Generated by \"#{scriptPath}\" from input file \"#{headerPath}\".")
	fortranLine(file, "");

	fortranLine(file, 'module ' + moduleName)
	fortranLine(file, "\tuse iso_c_binding")
	fortranLine(file, "\timplicit none")
	fortranLine(file, "\tprivate")

	$verbatimDeclarations.each_line do |line|
		fortranLine(file, line)
	end
	fortranLine(file, '')
	$declarationLines.each do |line|
		fortranLine(file, line)
	end
	fortranLine(file, '')
	unless $interfaceLines.empty?
		fortranLine(file, "\tinterface")
		$interfaceLines.each do |line|
			fortranLine(file, line)
		end
		fortranLine(file, "\tend interface")
		fortranLine(file, '')
	end
	fortranLine(file, "contains")
	$verbatimDefinitions.each_line do |line|
		fortranLine(file, line)
	end
	fortranLine(file, '')
	$definitionLines.each do |line|
		fortranLine(file, line)
	end

	fortranLine(file, 'end module ' + moduleName)
end

def main
	printUsage = false
	ARGV.each { |argument|
		if argument == "-h" || argument == "--help"
			printUsage = true
		end
	}
	unless printUsage
		case ARGV.length
			when 0
				$stderr.puts("Error: no input file given")
				printUsage = true
			when 1
				$stderr.puts("Error: no output file given")
				printUsage = true
			when 2
				moduleName = /(?<basename>[^.\/]+)\.[^\/]+/.match(ARGV[1])['basename']
			when 3
				moduleName = ARGV[2]
			else
				$stderr.puts("Error: too many arguments")
				printUsage = true
		end
	end
	unless printUsage
		headerPath = ARGV[0]
		outputPath = ARGV[1]
		scanHeader(headerPath)
		writeFortranModule($0, headerPath, outputPath, moduleName)
	else
		puts("Usage:")
		puts($0 + ' cHeader outputPath [ moduleName ]')
		puts($0 + ' ( -h | --help )')
		puts("")
		puts("\tcHeader:    input C header file")
		puts("\toutputPath: output fortran file name")
		puts("\tmoduleName: name of the resulting fortran module, defaults to the basename of outputPath")
	end
end

if rubyVersionOk()
	main()
else
	$stderr.puts("Error: Ruby version #{RUBY_VERSION} is too old (version 1.9 is required). Skipping fortran interface generation.")
end
#
# Local Variables:
# mode: ruby
# tab-always-indent: nil
# tab-width: 2
# ruby-indent-tabs-mode: t
# indent-tabs-mode: t
# End:
#

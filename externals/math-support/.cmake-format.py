# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
#
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: CC0-1.0
# ---------------------------------------------------------------

with section("parse"):
	additional_commands = {
	  "check_macro_defined": {
	  	"pargs": 2,
	  	"flags": ["QUIET"],
	  	"kwargs": {
	  		"LANG": 1
	  	}
	  },
	  "list_sources": {
	  	"pargs": 1,
	  	"flags": ["EXCLUDE_GENERATED"],
	  	"kwargs": {
	  		"DIRECTORY": 1,
	  		"INCLUDE_REGEX": 1
	  	}
	  },
	  "FortUTF_Find_Tests": {}
	}

with section("format"):
	line_width = 80
	autosort = True
	keyword_case = 'upper'

with section("markup"):
	first_comment_is_literal = True

with section("lint"):
	disabled_codes = ['C0301']
	function_pattern = '[0-9a-z_]+'
	macro_pattern = '[0-9a-z_]+'
	local_var_pattern = '[a-zA-Z][0-9a-zA-z_]+'
	private_var_pattern = '[a-z][a-z0-9_]+'
	public_var_pattern = '[A-Z][0-9a-zA-Z_]+'

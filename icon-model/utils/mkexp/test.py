# -*- coding: utf-8 -*-

import io
import locale
import os
import re
import subprocess
import sys
import unittest

import expconfig

from os.path import join


_preferred_encoding = locale.getpreferredencoding()

def align(string):
    return re.sub(r'\n\s*', '\n', string.lstrip())

def script(string, set_encoding=False):
    text = u"""
        set -e
        unset CDPATH
        cd test
        PATH=..:.:$PATH
        MKEXP_PATH=
    """
    if set_encoding and sys.getdefaultencoding() == 'ascii':
        text += u"""
            PYTHONIOENCODING={0}
            export PYTHONIOENCODING
        """.format(_preferred_encoding)
    text += string
    return align(text)

def output(command):
    p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (result, ignore) = p.communicate()
    try: result = str(result, _preferred_encoding)
    except TypeError: pass
    ### print("DEBUG:\n" + result) ###
    return result

def readfile(filename):
    with io.open(filename) as stream:
        return stream.read()

def writefile(filename, string):
    with io.open(filename, 'w') as stream:
        stream.write(string)

def writeconfig(exp_id, string):
    writefile(join("test", exp_id+".config"), string)

def writetemplate(exp_id, job_id, string):
    writefile(join("test", exp_id+"."+job_id+".tmpl"), string)



class MkexpTestCase(unittest.TestCase):

    script_clean = script(u"""
        rm -rf experiments
        rm -f test_* SETUP.config
    """)

    script_run = script(u"""
        mkexp test0001.config
    """)

    def setUp(self):
        os.system(self.script_clean)

    @classmethod
    def tearDownClass(cls):
        os.system(cls.script_clean)

class MkexpSimpleTestCase(MkexpTestCase):

    def setUp(self):
        self.exp_id = self.id().split('.')[-1]
        self.job_id = 'job'
        MkexpTestCase.setUp(self)

    def run_test(self, template, expected, additional='', epilog=''):
        writeconfig(self.exp_id, u"""
            EXP_TYPE =
            """+additional+u"""
            [jobs]
              [["""+self.job_id+u"""]]
            """+epilog+u"""
        """)
        writetemplate(self.exp_id, self.job_id, template)
        expected = align(expected)
        ignore = output(script("mkexp "+self.exp_id+".config"))
        result = readfile(join("test", "experiments", self.exp_id,
                               self.exp_id+"."+self.job_id))
        result = align(result)
        self.assertMultiLineEqual(expected, result)

    def run_no_template(self, result_path, expected, additional=''):
        writeconfig(self.exp_id, u"""
            EXP_TYPE =
            """+additional+u"""
            [jobs]
              [["""+self.job_id+u"""]]
        """)
        expected = align(expected)
        ignore = output(script("mkexp "+self.exp_id+".config"))
        result = readfile(join("test", "experiments", result_path))
        result = align(result)
        self.assertMultiLineEqual(expected, result)



class RunningTestCase(MkexpTestCase):

    def test_missing_config_file(self):
        result = output(script('mkexp'))
        self.assertTrue('error: too few arguments' in result or
            'error: the following arguments are required' in result)

    def test_clean_run(self):
        expected = align(u"""
            Script directory: 'experiments/test0001'
            Data directory: 'experiments/test0001' (already exists)
            Work directory: 'experiments/test0001' (already exists)
        """)
        result = output(self.script_run)
        self.assertMultiLineEqual(expected, result)

    def test_backup_run(self):
        expected = "Note: script directory already exists, "+\
                   "moving existing scripts to backup"
        ignore = output(self.script_run)
        result = output(self.script_run)
        self.assertIn(expected, result)

    def test_additional_dirs_run(self):
        expected = align(u"""
            Script directory: 'experiments/test0001'
            Data directory: 'experiments/test0001' (already exists)
            Work directory: 'experiments/test0001' (already exists)
            Log directory: 'experiments/test0001/log'
        """)
        result = output(script("mkexp test0001.config "
            r'LOG_DIR=\$DATA_DIR/log EXP_DIR_NAMES=LOG_DIR,WORK_DIR,LOG_DIR'))
        self.assertMultiLineEqual(expected, result)

class CommandLineTestCase(MkexpTestCase):

    def test_pass_section_variables(self):
        script_section = script(u"""
            mkexp test0001.config \
                namelists.namelist..echam.runctl.dt_start=2345,01,23,12,34,56 \
                namelists.namelist..echam.runctl.some_file=abcdefgh.ijk
        """)
        expecteds = ["dt_start = 2345, 1, 23, 12, 34, 56",
                     "some_file = 'abcdefgh.ijk'"]
        ignore = output(script_section)
        result = readfile('test/experiments/test0001/test0001.run')
        for expected in expecteds:
            self.assertIn(expected, result)

    def test_pass_new_job(self):
        output(script("mkexp test0001.config jobs.dummy...extends=run"))
        readfile('test/experiments/test0001/test0001.dummy')
        # Should exist, otherwise exception is thrown

    def test_options(self):
        script_option = script(u"""
            mkexp test0001.config EXP_OPTIONS=option1
        """)
        expected = "default_output = .false."
        ignore = output(script_option)
        result = readfile('test/experiments/test0001/test0001.run')
        self.assertIn(expected, result)

    def test_getexp_vv(self):
        script_getexp = script(u"""
            mkexp test0001.config
            mv experiments/test0001 experiments/test0001.orig
            getexp -vv test0001.config MODEL_DIR=. > test_getexp.dump
            mkexp --getexp test_getexp.dump
        """, set_encoding=True)
        ignore = output(script_getexp)
        expected = readfile('test/experiments/test0001.orig/test0001.run')
        result = readfile('test/experiments/test0001/test0001.run')
        self.assertMultiLineEqual(expected, result)

    def test_getexp_k(self):
        result = output(script('getexp -k VAR1 test0001.config'))
        expected = align(u"""
            Note: data for experiment 'test0001' does not exist
            value1
        """)
        self.assertMultiLineEqual(expected, result)

    def test_getexp_k_k(self):
        result = output(script('getexp -k VAR1 -k VAR2 test0001.config'))
        expected = align(u"""
            Note: data for experiment 'test0001' does not exist
            value1
            value1
        """)
        self.assertMultiLineEqual(expected, result)

    def test_getexp_s(self):
        result = output(script('getexp -s -k VAR1 test0001.config'))
        expected = align(u"""
            Note: data for experiment 'test0001' does not exist
            VAR1='value1'
        """)
        self.assertMultiLineEqual(expected, result)

    def test_getconfig(self):
        ignore = output(script(u"""
            mkexp test0001.config VAR4=value4 jobs.run.time_limit=12:34:56
        """))
        result = output(script('getconfig experiments/test0001/update',
                               set_encoding=True))
        self.assertIn('VAR4 = value4', result)
        self.assertIn('time_limit = 12:34:56', result)

    def test_setconfig(self):
        result = output(script(u"""
            setconfig test0001.config VAR4=value4 jobs.run.time_limit=12:34:56
        """, set_encoding=True))
        self.assertIn('VAR4 = value4', result)
        self.assertIn('time_limit = 12:34:56', result)

    def test_selconfig(self):
        result = output(script(u"""
            selconfig VAR1 test0001.config
        """))
        self.assertIn('VAR1 = value1', result)
        self.assertNotIn('VAR2', result)

    def test_rmexp(self):
        script_getexp = script(u"""
            mkexp test0001.config
            (echo n; echo y) | rmexp test0001.config
        """)
        result = output(script_getexp)
        self.assertIn("Script directory: 'experiments/test0001'", result)
        self.assertIn("Data directory: 'experiments/test0001'", result)
        self.assertIn("Info for test0001's working directory:", result)
        self.assertIn("rmexp: remove test0001's working directory and its contents [no]? Info for test0001's script directory:", result)
        self.assertNotIn("rmexp: remove test0001's working directory and its contents [no]? rmexp: removed 'experiments/test0001'", result)
        self.assertIn("rmexp: remove test0001's script directory and its contents [no]? rmexp: removed 'experiments/test0001'", result)
        self.assertNotIn("rmexp: remove test0001's script directory and its contents [no]? Info for test0001's data directory:", result)
        self.assertNotIn("rmexp: remove test0001's data directory and its contents [no]?", result)

class ContentTestCase(MkexpSimpleTestCase):

    def test_job_override(self):
        exp_id = "test_job_override"
        writeconfig(exp_id, u"""
            EXP_TYPE =
            [jobs]
              key1 = global
              key2 = global
              [[job1]]
                key1 = local
        """)
        writetemplate(exp_id, "job1", u"""
            key1 = %{JOB.key1}
            key2 = %{JOB.key2}
        """)
        ignore = output(script("mkexp "+exp_id+".config"))
        result = readfile(join("test", "experiments", exp_id, exp_id+".job1"))
        self.assertIn("key1 = local", result)
        self.assertIn("key2 = global", result)

    def test_var_list_in_context(self):
        exp_id = "test_var_list_in_context"
        job_id = "job"
        writeconfig(exp_id, u"""
            EXP_TYPE =
            VAR1 = value1
            GLOBAL1 = $${VAR1} # Initialized
            GLOBAL2 = $${VAR2} # Uninitialized
            GLOBAL3 = $${VAR1} # Used twice, may only be listed once
            GLOBAL${FOUR} = 4  # (Uninitialized) Variable in key
            [jobs]
              [["""+job_id+u"""]]
        """)
        writetemplate(exp_id, job_id, u"""
            #% for var in VARIABLES_|sort:
            %{var}=%{context(var)}
            #% endfor
        """)
        expected = align(u"""
            FOUR=
            VAR1=value1
            VAR2=
        """)
        ignore = output(script("mkexp "+exp_id+".config"))
        result = readfile(join("test", "experiments", exp_id, exp_id+"."+job_id))
        result = align(result)
        self.assertMultiLineEqual(expected, result)

    def test_split_date(self):
        exp_id = 'test_split_date'
        job_id = 'job'
        writeconfig(exp_id, u"""
            EXP_TYPE =
            DATE_ISO = 1234-05-06
            DATE_RAW = 12340506
            DATE_LIST_ISO = split_date($DATE_ISO)
            DATE_LIST_RAW = split_date($DATE_RAW)
            [jobs]
              [["""+job_id+u"""]]
        """)
        writetemplate(exp_id, job_id, u"""
            %{DATE_LIST_ISO|join(',')}
            %{DATE_LIST_RAW|join(',')}
        """)
        expected = align(u"""
            1234,5,6,0,0,0
            1234,05,06,0,0,0
        """)
        ignore = output(script("mkexp "+exp_id+".config"))
        result = readfile(join("test", "experiments", exp_id, exp_id+"."+job_id))
        result = align(result)
        self.assertMultiLineEqual(expected, result)

    def test_add_years(self):
        exp_id = 'test_add_years'
        job_id = 'job'
        writeconfig(exp_id, u"""
            EXP_TYPE =
            DATE = 1234-05-06
            NEXT_DATE = 'add_years($DATE, 1)'
            PREVIOUS_DATE = 'add_years($DATE, -1)'
            NEGATIVE_DATE = 'add_years($DATE, -2000)'
            LONGYEAR_DATE = 'add_years($DATE, 10000)'
            [jobs]
              [["""+job_id+u"""]]
        """)
        writetemplate(exp_id, job_id, u"""
            %{NEXT_DATE}
            %{PREVIOUS_DATE}
            %{NEGATIVE_DATE}
            %{LONGYEAR_DATE}
        """)
        expected = align(u"""
            1235-05-06
            1233-05-06
            -0766-05-06
            11234-05-06
        """)
        ignore = output(script("mkexp "+exp_id+".config"))
        result = readfile(join("test", "experiments", exp_id, exp_id+"."+job_id))
        result = align(result)
        self.assertMultiLineEqual(expected, result)

    def test_add_days(self):
        exp_id = 'test_add_days'
        job_id = 'job'
        writeconfig(exp_id, u"""
            EXP_TYPE =
            DATE = 1234-05-06
            NEXT_DATE = 'add_days($DATE, 1)'
            PREVIOUS_DATE = 'add_days($DATE, -1)'
            NEGATIVE_DATE = 'add_days($DATE, -2000)'
            LONGYEAR_DATE = 'add_days($DATE, 10000)'
            LATE_DATE = 9999-12-31
            LATER_DATE = 'add_days($LATE_DATE, 1)'
            EARLY_DATE = 0000-01-01
            EARLIER_DATE = 'add_days($EARLY_DATE, -1)'
            [jobs]
              [["""+job_id+u"""]]
        """)
        writetemplate(exp_id, job_id, u"""
            %{NEXT_DATE}
            %{PREVIOUS_DATE}
            %{NEGATIVE_DATE}
            %{LONGYEAR_DATE}
            %{LATER_DATE}
            %{EARLIER_DATE}
        """)
        expected = align(u"""
            1234-05-07
            1234-05-05
            1228-11-13
            1261-09-21
            10000-01-01
            -0001-12-31
        """)
        ignore = output(script("mkexp "+exp_id+".config"))
        result = readfile(join("test", "experiments", exp_id, exp_id+"."+job_id))
        result = align(result)
        self.assertMultiLineEqual(expected, result)

    def test_eval(self):
        self.run_test(u"""
            %{VALUE}
        """, u"""
            42
        """, u"""
            VALUE = eval(5*8+2)
        """)

    def test_eval_time(self):
        self.run_test(u"""
            %{VALUE}
        """, u"""
            1970-01-01
        """, u"""
            VALUE = "eval(time.strftime('%Y-%m-%d', time.gmtime(0)))"
        """)

    def test_eval_is_set(self):
        self.run_test(u"""
            %{TRUE}
            %{FALSE}
        """, u"""
            True
            False
        """, u"""
            TRUE = eval(is_set('.true.'))
            FALSE = eval(is_set('F'))
        """)

    def test_initial_comment_boilerplate(self):
        writeconfig(self.exp_id, u"""
            ######
            #    #
            # 42
            #    #
            ######
            EXP_TYPE =
            [jobs]
              [["""+self.job_id+"""]]
        """)
        writetemplate(self.exp_id, self.job_id, u"""
            %{EXP_DESCRIPTION}
        """)
        expected = align("""
            42
        """)
        ignore = output(script("mkexp "+self.exp_id+".config"))
        result = readfile(join("test", "experiments", self.exp_id,
                               self.exp_id+"."+self.job_id))
        result = align(result)
        self.assertMultiLineEqual(expected, result)

class StaticAssignmentTestCase(MkexpSimpleTestCase):

    def test_add_to_list(self):
        self.run_test(u"""
            %{LIST1}
            %{LIST2}
            %{LIST3}
        """, u"""
            ['a']
            ['b', 'c']
            ['d', 'e', 'f', 'g']
        """, u"""
            LIST1 += a
            LIST2 = b
            LIST2 += c
            LIST3 = d, e
            LIST3 + = f, g
        """)

    def test_remove_from_list(self):
        self.run_test(u"""
            %{LIST1}
            %{LIST2}
            %{LIST3}
        """, u"""
            []
            ['b', 'c']
            ['d', 'f']
        """, u"""
            LIST1 = a
            LIST1 -= a
            LIST2 = b, c
            LIST2 -= z
            LIST3 = d, e, f, g
            LIST3 - = e, g
        """)

    def test_order_add_remove_list(self):
        self.run_test(u"""
            %{LIST1}
            %{LIST2}
        """, u"""
            ['a', 'b']
            ['a', 'b', 'c']
        """, u"""
            LIST1 = a, b
            LIST1 += c
            LIST1 -= c
            LIST2 = a, b
            LIST2 -= c
            LIST2 += c
        """)

    def test_rename_list_items(self):
        self.run:test(u"""
            %{LIST}
        """, u"""
            ['f', 'de', 'de']
        """, u"""
            LIST = a, ab, abc
            LIST >= ab>de, a>f, c>
        """)

    def test_invalid_rename(self):
        writeconfig(self.exp_id, u"""
            EXP_TYPE = DEFAULT
            LIST >= a
        """)
        result = output(script("mkexp "+self.exp_id+".config")) 
        self.assertEqual(result.rstrip(),
            "Oops: invalid rename 'a' while reading key 'LIST >'")
    

class NamelistTestCase(MkexpSimpleTestCase):

    def test_namelist_comments(self):
        self.run_test(u"""
            %{NAMELIST}
        """,u"""
            ! Comment group 1
            ! var_1c = 'test'
            &group_1
                ! Comment for var 1a
                var_1a = 42 ! Inline comment for var 1a
                ! var_1b = .true.
            /
            ! var_1d = 10.5
            &group_2
                ! Comment for var 2b
                var_2b = 21 ! Inline comment for var 2b
            /
        """,u"""
            [namelists]
              [[namelist]]
                # Comment group 1
                # var_1c = test
                [[[group_1]]]
                  # Comment for var 1a
                  var_1a = 42 # Inline comment for var 1a
                  # var_1b = true
                  .end =
                  # var_1d = 10.5
                [[[group_2]]]
                  # Comment for var 2b
                  var_2b = 21 # Inline comment for var 2b
        """)

    def test_var_in_namelist(self):
        self.run_test(u"""
            %{NAMELIST}
        """,u"""
            &group
                var_1 = $value_1
                var_2 = ${value_2}
                var_3 = 'a', $value_3, 'b'
                var_4 = 'a$value_4'
                var_5 = '${value_5}b'
            /
        """,u"""
            [namelists]
              [[namelist]]
                [[[group]]]
                  var_1 = raw($$value_1)
                  var_2 = raw($${value_2})
                  var_3 = a, raw($$value_3), b
                  var_4 = a$$value_4
                  var_5 = $${value_5}b
        """)

    def test_namelist_multi_groups(self):
        self.run_test(u"""
            %{NAMELIST}
        """, u"""
            &group ! '1'
            /
            &group ! ' 1'
            /
            &group ! '2'
            /
            &group ! 'i i i'
            /
        """, u"""
            [namelists]
              [[namelist]]
                [[[group 1]]]
                [[[group  1]]]
                [[[group 2]]]
                [[[group i i i]]]
        """)

    def test_namelist_case_twist(self):
        self.run_test(u"""
            %{NAMELIST}
        """, u"""
            &group
                value = 41
                value = 42
            /
        """, u"""
            [namelists]
              [[namelist]]
                [[[group]]]
                   value = 41
                   VALUE = 42
        """)

    def test_namelist_format(self):
        self.run_test(u"""
            %{format_namelist(namelists.namelist)}
            %{format_namelist(namelists.namelist, 'group2')}
            %{format_namelist(namelists.namelist, 'no such group')}
        """, u"""
            &group1
                value = 41
            /
            &group2
                value = 42
            /
            &group2
                value = 42
            /
        """, u"""
            [namelists]
              [[namelist]]
                .remove = group3
                [[[group1]]]
                   value = 41
                [[[group2]]]
                   value = 42
                [[[group3]]]
                   value = 43
        """)

class NamelistHiddenTestCase(MkexpSimpleTestCase):

    def test_namelist_hide(self):
        self.run_test(u"""
            %{NAMELIST}
        """, u"""
            &group1
                value = 41
            /
            &group3
                value = 43
            /
        """, u"""
            [namelists]
              [[namelist]]
                [[[group1]]]
                   value = 41
                [[[group2]]]
                   .hide = true
                   value = 42
                [[[group3]]]
                   .hide = dont_care_if_we_dont_start_with_t
                   value = 43
        """)

    def test_hidden_namelist_file(self):
        self.run_test(u"""
            %{NAMELIST}
        """, u"""
        """, u"""
            [namelists]
              [[namelist]]
                .hide = true
                [[[group]]]
                   value = 42
        """)

    def test_hidden_namelist_with_template(self):
        self.run_test(u"""
            %{NAMELIST}
        """, u"""
        """, u"""
            [namelists]
              [[namelist]]
                .hide = true
                .use_template = true
                [[[group]]]
                   value = 42
        """)

class NamelistDefaultValueTestCase(MkexpSimpleTestCase):

    def test_standard_default_value(self):
        self.run_test(u"""
            %{NAMELIST}
        """, u"""
            &group
            /
        """, u"""
            [namelists]
              [[namelist]]
                [[[group]]]
                  value =
        """)

    def test_custom_default_value(self):
        self.run_test(u"""
            %{NAMELIST}
        """, u"""
            &group
                value1 = ''
            /
        """, u"""
            [namelists]
              [[namelist]]
                .default = <DEFAULT>
                [[[group]]]
                  value1 =
                  value2 = ${.default}
        """)

    def test_namelist_default_value(self):
        self.run_test(u"""
            %{NAMELIST1}
            %{NAMELIST2}
        """, u"""
            &group
            /
            &group
                value = '<DEFAULT>'
            /
        """, u"""
            [namelists]
              [[namelist1]]
                .default = <DEFAULT>
                [[[group]]]
                  value = <DEFAULT>
              [[namelist2]]
                [[[group]]]
                  value = <DEFAULT>
        """)

    def test_global_default_value(self):
        self.run_test(u"""
            %{format_namelist(namelists.namelist, default_value='<DEFAULT>')}
        """, u"""
            &group
            /
        """, u"""
            [namelists]
              [[namelist]]
                [[[group]]]
                  value = <DEFAULT>
        """)

class NamelistInheritanceTestCase(MkexpSimpleTestCase):

    def test_basic_inheritance(self):
        self.run_test(u"""
            %{NAMELIST}
        """, u"""
            &group1
                value = 42
            /
            &group2
                value = 42
            /
        """, u"""
            [namelists]
              [[namelist]]
                [[[group1]]]
                  value = 42
                [[[group2]]]
                  .extends = group1
        """)

    def test_hidden_inheritance(self):
        self.run_test(u"""
            %{NAMELIST}
        """, u"""
            &group2
                value = 42
            /
        """, u"""
            [namelists]
              [[namelist]]
                .remove = group1
                [[[group1]]]
                  value = 42
                [[[group2]]]
                  .extends = group1
        """)

    def test_inheritance_by_name(self):
        self.run_test(u"""
            %{NAMELIST}
        """, u"""
            &group
                value = 42
            /
            &group ! 'clone'
                value = 42
            /
        """, u"""
            [namelists]
              [[namelist]]
                [[[group]]]
                  value = 42
                [[[group clone]]]
        """)

    def test_id_for_inheritance_by_name(self):
        self.run_test(u"""
            %{NAMELIST}
        """, u"""
            &group
                value = ''
                list = 1, "my id is ''", 3
            /
            &group ! 'clone'
                value = 'clone'
                list = 1, "my id is 'clone'", 3
            /
        """, u"""
            [namelists]
              [[namelist]]
                [[[group]]]
                  value = %{id}
                  list = 1, my id is '%{id}', 3
                [[[group clone]]]
        """)

class NamelistTemplateTestCase(MkexpSimpleTestCase):

    def test_use_template(self):
        self.run_test(u"""
            %{NAMELIST}
        """, u"""
            42
        """, u"""
            [namelists]
              [[namelist]]
                .use_template = true
                [[[group]]]
                   value = 42
        """)

    def test_named_template(self):
        writetemplate(self.exp_id, 'namelist_template', u"""
            %{group.value}
        """)
        self.run_test(u"""
            %{NAMELIST}
        """, u"""
            42
        """, u"""
            [namelists]
              [[namelist]]
                .use_template = namelist_template
                [[[group]]]
                   value = 42
        """)

    def test_uppercase_template_name(self):
        writetemplate(self.exp_id, 'UPPERCASE_namelist_template', u"""
            %{group.value}
        """)
        self.run_test(u"""
            %{NAMELIST}
        """, u"""
            42
        """, u"""
            [namelists]
              [[namelist]]
                .use_template = UPPERCASE_namelist_template
                [[[group]]]
                   value = 42
        """)

    def test_full_config_in_template(self):
        writetemplate(self.exp_id, 'namelist_template', u"""
            %{EXP_ID}
            %{value}
            %{_.value}
            %{namelists.namelist.value}
            %{namelists.namelist.group.value}
        """)
        self.run_test(u"""
            %{NAMELIST}
        """, u"""
            """+self.exp_id+"""
            namelist
            global
            namelist
            group
        """, u"""
            value = global
            [namelists]
              [[namelist]]
                .use_template = namelist_template
                .use_full_config = true
                value = namelist
                [[[group]]]
                   value = group 
        """)

    def test_normal_config_in_template(self):
        writetemplate(self.exp_id, 'namelist_template', u"""
#%  for key, value in context().get_all().items():
#%    if not key.startswith('.') and not value is callable:
%{key}
#%    endif
#%  endfor
        """)
        self.run_test(u"""
            %{NAMELIST}
        """, u"""
            value
            group
        """, u"""
            value = global
            [namelists]
              [[namelist]]
                .use_template = namelist_template
                value = namelist
                [[[group]]]
                   value = group 
        """)


class KeySplittingTestCase(MkexpSimpleTestCase):

    def test_job_splitting(self):
        self.run_test(u"""
            %{{jobs.{job_id}.seniority}}
            %{{jobs.{job_id}.common}}
            %{{jobs.job2.seniority}}
            %{{jobs.job2.common}}
        """.format(**self.__dict__), u"""
            elder
            parents
            younger
            parents
        """, epilog="""
              seniority = elder
            [[{job_id}, job2]]
              common = parents
            [[job2]]
              .extends = {job_id}
              seniority = younger
        """.format(**self.__dict__))

    def test_namelist_splitting(self):
        self.run_test(u"""
            %{NAMELIST_ATM}
            %{NAMELIST_OCE}
        """, u"""
            &group1
                value = 21
            /
            &group2
                value = 42
            /
            &group3
                value = 84
            /
            &group2
                value = 42
            /
        """, u"""
            [namelists]
              [[NAMELIST_atm]]
                [[[group1]]]
                  value = 21
              [[NAMELIST_atm, NAMELIST_oce]]
                [[[group2]]]
                  value = 42
              [[NAMELIST_oce]]
                [[[group3]]]
                  value = 84
        """)

    def test_job_namelist_splitting(self):
        self.run_test(u"""
            %{NAMELIST_ATM}
            %{NAMELIST_OCE}
        """, u"""
            &group1
                value = 21
            /
            &group2
                value = 42
            /
            &group3
                value = 84
            /
            &group2
                value = 42
            /
        """, u"""
            [namelists]
              [[NAMELIST_atm]]
                [[[group1]]]
                  value = 21
              [[NAMELIST_oce]]
                [[[group3]]]
                  value = 84
        """, u"""
            [[[namelists]]]
              [[[[NAMELIST_atm, NAMELIST_oce]]]]
                [[[[[group2]]]]]
                  value = 42
        """)


class JinjaTemplateTestCase(MkexpSimpleTestCase):

    def test_ignore_blocks(self):
        self.run_test(u"""
            {% set answer = 42 %}
        """, u"""
            {% set answer = 42 %}
        """)

    def test_ignore_comments(self):
        self.run_test(u"""
            {# no comment #}
            ${#ARRAY}
        """, u"""
            {# no comment #}
            ${#ARRAY}
        """)

class DefaultEnvironmentTestCase(MkexpSimpleTestCase):

    def test_basic(self):
       self.run_test(u"""
           %{ENVIRONMENT}
       """, u"""
           DEFAULT
       """)

    def test_explicit(self):
       self.run_test(u"""
           %{ENVIRONMENT}
       """, u"""
           green
       """, u"""
           ENVIRONMENT = green
       """)

    def test_setup(self):
       writeconfig('SETUP', u"""
           ENVIRONMENT = 
       """)
       self.run_test(u"""
           %{ENVIRONMENT}
       """, u"""
           DEFAULT
       """)

class SetupConfigTestCase(MkexpSimpleTestCase):

    def test_system_options(self):
       writeconfig('SETUP', u"""
           SETUP_OPTIONS = option1
       """)
       self.run_test(u"""
           %{NAMELIST_ECHAM}
       """, u"""
           &runctl
             default_output = .false.
           /
       """, u"""
           EXP_OPTIONS =
       """)
           
class MatchTestCase(MkexpSimpleTestCase):

    def test_basic(self):
       self.run_test(u"""
           %{'Douglas Adams'|match('Adam')}
       """, u"""
           Douglas Adams
       """)

    def test_no_match(self):
       self.run_test(u"""
           %{'Douglas Adams'|match('Eva')}
       """, u"""
           
       """)

    def test_with_default(self):
       self.run_test(u"""
           %{'Douglas Adams'|match('Abel', 'Kain')}
       """, u"""
           Kain
       """)

    def test_with_group(self):
       self.run_test(u"""
           %{'Douglas Adams'|match('l(.*)m')}
       """, u"""
           as Ada
       """)

class SplitTestCase(MkexpSimpleTestCase):

    def test_basic(self):
        self.run_test(u"""
            %{'Douglas Noel  Adams'|split(' ')}
            %{'Douglas Noel  Adams'|split('  ')}
        """,u"""
            ['Douglas', 'Noel', '', 'Adams']
            ['Douglas Noel', 'Adams']
        """)

    def test_max_split(self):
        self.run_test(u"""
            %{'Douglas Noel  Adams'|split(' ', 2)}
            %{'Douglas Noel  Adams'|split(' ', 1)}
        """,u"""
            ['Douglas', 'Noel', ' Adams']
            ['Douglas', 'Noel  Adams']
        """)

    def test_default_separator(self):
        self.run_test(u"""
            %{'Douglas Noel  Adams'|split()}
            %{'Douglas Noel  Adams'|split(none)}
        """,u"""
            ['Douglas', 'Noel', 'Adams']
            ['Douglas', 'Noel', 'Adams']
        """)

    def test_max_split_and_default_separator(self):
        self.run_test(u"""
            %{'Douglas Noel  Adams'|split(m=2)}
            %{'Douglas Noel  Adams'|split(m=1)}
        """,u"""
            ['Douglas', 'Noel', 'Adams']
            ['Douglas', 'Noel  Adams']
        """)

class FilterTestCase(MkexpSimpleTestCase):

    def test_basic(self):
        self.run_test(u"""
            %{['Douglas', 'Noel', '', 'Adams']|filter}
        """,u"""
            ['Douglas', 'Noel', 'Adams']
        """)

class WordwrapTestCase(MkexpSimpleTestCase):

    def test_basic(self):
        self.run_test(u"""
            %{'long-arbitrarilyhyphenated textlike-message'|wordwrap(15)}
        """, u"""
            long-arbitraril
            yhyphenated
            textlike-
            message
        """)

    def test_keep_long(self):
        self.run_test(u"""
            %{'long-arbitrarilyhyphenated textlike-message'|wordwrap(15,false)}
        """, u"""
            long-
            arbitrarilyhyphenated
            textlike-
            message
        """)

    def test_keep_long_hyphens(self):
        self.run_test(u"""
            %{'long-arbitrarilyhyphenated textlike-message'|wordwrap(15,false,false)}
        """, u"""
            long-arbitrarilyhyphenated
            textlike-message
        """)

    def test_keep_hyphens(self):
        self.run_test(u"""
            %{'long-arbitrarilyhyphenated textlike-message'|wordwrap(15,true,false)}
        """, u"""
            long-arbitraril
            yhyphenated tex
            tlike-message
        """)

class ListTestCase(MkexpSimpleTestCase):

    def test_list_on_string(self):
        self.run_test(u"""
            %{'first'|list}
        """, u"""
            ['first']
        """)

    def test_list_on_empty_string(self):
        self.run_test(u"""
            %{''|list}
        """, u"""
            []
        """)

    def test_list_keep_empty_string(self):
        self.run_test(u"""
            %{''|list(true)}
        """, u"""
            ['']
        """)

    def test_list_on_list(self):
        self.run_test(u"""
            %{['first', 'second', 'third']|list}
        """, u"""
            ['first', 'second', 'third']
        """)

    def test_list_on_int(self):
        self.run_test(u"""
            %{42|list}
        """, u"""
            [42]
        """)

    def test_list_on_tuple(self):
        self.run_test(u"""
            %{('first', 'second', 'third')|list}
        """, u"""
             ['first', 'second', 'third']
        """)

    def test_list_on_iterator(self):
        self.run_test(u"""
            %{('first', 'second', 'third')|reverse|list}
        """, u"""
             ['third', 'second', 'first']
        """)

class JoinTestCase(MkexpSimpleTestCase):

    def test_join_on_string(self):
        self.run_test(u"""
            %{'first'|join(', ')}
        """, u"""
            first
        """)

    def test_join_on_empty_string(self):
        self.run_test(u"""
            %{''|join}
        """, u"""
        """)

    def test_join_on_list(self):
        self.run_test(u"""
            %{['first', 'second', 'third']|join(', ')}
        """, u"""
            first, second, third
        """)

    def test_join_on_int(self):
        self.run_test(u"""
            %{42|join(', ')}
        """, u"""
            42
        """)

    def test_join_on_tuple(self):
        self.run_test(u"""
            %{('first', 'second', 'third')|join(', ')}
        """, u"""
             first, second, third
        """)

    def test_join_on_iterator(self):
        self.run_test(u"""
            %{('first', 'second', 'third')|reverse|join(', ')}
        """, u"""
             third, second, first
        """)

class IsSetTestCase(MkexpSimpleTestCase):

    def test_empty_string(self):
        self.run_test(u"""
            %{'' is set}
        """, u"""
            False
        """)

    def test_true(self):
        self.run_test(u"""
            %{'true' is set}
        """, u"""
            True
        """)

    def test_namelist_true(self):
        self.run_test(u"""
            %{'.true.' is set}
        """, u"""
            True
        """)

    def test_false(self):
        self.run_test(u"""
            %{'false' is set}
        """, u"""
            False
        """)

    def test_test(self):
        self.run_test(u"""
            %{'.test.' is set}
        """, u"""
            True
        """)

    def test_undefined(self):
        self.run_test(u"""
            %{undefined_variable_name is set}
        """, u"""
            False
        """)

    def test_undefined_with_true_default(self):
        self.run_test(u"""
            %{undefined_variable_name|d('t') is set}
        """, u"""
            True
        """)

class FilesTestCase(MkexpSimpleTestCase):

    def test_get_file_simple(self):
        self.run_test(u"""
            %{get_file(files, 'target.txt')}
            %{get_file(files, 'broken.txt')}
        """, u"""
            source.txt
            .
        """, u"""
            [files]
                target.txt = source.txt
                broken.txt = .
        """)

    def test_get_file_path(self):
        self.run_test(u"""
            %{get_file(files, 'target.txt')}
            %{get_file(files, 'path.txt')}
            %{get_file(files.subdir, 'target.txt')}
        """, u"""
            /path/to/source/source.txt
            /just/this/one/source.txt
            /path/to/source/subdir/source.txt
        """, u"""
            [files]
                .base_dir = /path/to/source
                target.txt = source.txt
                path.txt = /just/this/one/source.txt
                [[subdir]]
                    .sub_dir = subdir
                    target.txt = source.txt
        """)

    def test_get_file_variable(self):
        self.run_test(u"""
            %{JOB.id}
            %{get_file(files, 'target.txt')}
            %{get_file(files, 'broken.txt')}
            %{get_file(files, 'incomplete.txt')}
        """, u"""
            job
            source.txt
            $BASENAME.txt
            ${DOES_NOT_EXIST}.txt
        """, u"""
            BASENAME = source
            [files]
                target.txt = $${BASENAME}.txt
                broken.txt = $$BASENAME.txt
                incomplete.txt = $${DOES_NOT_EXIST}.txt
        """)

    def test_get_dir(self):
        self.run_test(u"""
            %{get_dir(files)}
            %{get_dir(files.subdir)}
        """, u"""
            /path/to/source
            /path/to/source/subdir
        """, u"""
            [files]
                .base_dir = /path/to/source
                [[subdir]]
                    .sub_dir = subdir
        """)

class GetTemplatesTestCase(MkexpSimpleTestCase):

    def test_by_config_file_name(self):
        other_exp_id = 'test_something_completely_different'
        writetemplate(self.exp_id, self.job_id, u"""
            selected by config file name
        """)
        self.run_no_template(join(other_exp_id, other_exp_id+'.'+self.job_id),
        u"""
            selected by config file name
        """, u"""
            EXP_ID = """+other_exp_id+u"""
        """)

    def test_by_exp_id(self):
        other_exp_id = 'test_something_completely_different'
        writetemplate(other_exp_id, self.job_id, u"""
            selected by EXP_ID
        """)
        self.run_no_template(join(other_exp_id, other_exp_id+'.'+self.job_id),
        u"""
            selected by EXP_ID
        """, u"""
            EXP_ID = """+other_exp_id+u"""
        """)

class DelimiterTestCase(MkexpSimpleTestCase):

    def test_statement(self):
        self.run_test(u"""
            {%__mkexp__
                set x = 'Hello, world!'
            %}
            %{x}
        """, u"""
            Hello, world!
        """)

    def test_comment(self):
        self.run_test(u"""
            {#__mkexp__
                Now you see me - now you don't
            #}
        """, u"""
        """)

class InheritanceTestCase(MkexpSimpleTestCase):

    def test_child_template(self):
        writeconfig(self.exp_id, u"""
            EXP_TYPE =
            [jobs]
              [[job1]]
              [[job2]]
                .extends = job1
        """)
        writetemplate(self.exp_id, 'job1', u"""
            %{JOB.id} as in job1
        """)
        writetemplate(self.exp_id, 'job2', u"""
            %{JOB.id} as in job2
        """)
        expected = align(u"""
            job2 as in job2
        """)
        ignore = output(script("mkexp "+self.exp_id+".config"))
        result = readfile(join("test", "experiments", self.exp_id,
                               self.exp_id+".job2"))
        result = align(result)
        self.assertMultiLineEqual(expected, result)

    def test_parent_template(self):
        writeconfig(self.exp_id, u"""
            EXP_TYPE =
            [jobs]
              [[job1]]
              [[job2]]
                .extends = job1
        """)
        writetemplate(self.exp_id, 'job1', u"""
            %{JOB.id} as in job1
        """)
        expected = align(u"""
            job2 as in job1
        """)
        ignore = output(script("mkexp "+self.exp_id+".config"))
        result = readfile(join("test", "experiments", self.exp_id,
                               self.exp_id+".job2"))
        result = align(result)
        self.assertMultiLineEqual(expected, result)

    def test_grandparent_template(self):
        writeconfig(self.exp_id, u"""
            EXP_TYPE =
            [jobs]
              [[job1]]
              [[job2]]
                .extends = job1
              [[job3]]
                .extends = job2
        """)
        writetemplate(self.exp_id, 'job1', u"""
            %{JOB.id} as in job1
        """)
        expected = align(u"""
            job3 as in job1
        """)
        ignore = output(script("mkexp "+self.exp_id+".config"))
        result = readfile(join("test", "experiments", self.exp_id,
                               self.exp_id+".job3"))
        result = align(result)
        self.assertMultiLineEqual(expected, result)

    def test_variable_ancestry(self):
        writeconfig(self.exp_id, u"""
            EXP_TYPE =
            [jobs]
              var_0 = from jobs
              [[job1]]
                var_1 = from job1
                var_2 = from job1
                var_3 = from job1
              [[job2]]
                .extends = job1
                var_2 = from job2
                var_3 = from job2
              [[job3]]
                .extends = job2
                var_0 = not needed
                var_3 = from job3
        """)
        writetemplate(self.exp_id, 'job1', u"""
            %{JOB.id}
            %{JOB.var_0}
            %{JOB.var_1}
            %{JOB.var_2}
            %{JOB.var_3}
        """)
        expecteds = list(map(align, (
            u"""
                job1
                from jobs
                from job1
                from job1
                from job1
            """,
            u"""
                job2
                from jobs
                from job1
                from job2
                from job2
            """,
            u"""
                job3
                not needed
                from job1
                from job2
                from job3
            """)))
        ignore = output(script("mkexp "+self.exp_id+".config"))
        for i in (1, 2, 3):
            result = readfile(join("test", "experiments", self.exp_id,
                                   self.exp_id+".job%d"%i))
            result = align(result)
            self.assertMultiLineEqual(expecteds[i-1], result)

    def test_namelist_override(self):
        writeconfig(self.exp_id, u"""
            EXP_TYPE =
            [namelists]
              [[namelist]]
                [[[group]]]
                  var = 999
            [jobs]
              [[job]]
              [[job1]]
                .extends = job
                [[[namelists]]]
                  [[[[namelist]]]]
                    [[[[[group]]]]]
                      var = 1
              [[job2]]
                .extends = job1
                [[[namelists]]]
                  [[[[namelist]]]]
                    [[[[[group]]]]]
                      var = 2
        """)
        writetemplate(self.exp_id, 'job', u"""
            %{NAMELIST}
        """)
        expecteds = {
            'job':  align(u"""
                        &group
                            var = 999
                        /
                    """),
            'job1': align(u"""
                        &group
                            var = 1
                        /
                    """),
            'job2': align(u"""
                        &group
                            var = 2
                        /
                    """)
        }
        ignore = output(script("mkexp "+self.exp_id+".config"))
        for i in ('job', 'job1', 'job2'):
            result = readfile(join("test", "experiments", self.exp_id,
                                   self.exp_id+"."+i))
            result = align(result)
            self.assertMultiLineEqual(expecteds[i], result)

    def test_job_dict_in_extend(self):
        writeconfig(self.exp_id, u"""
            EXP_TYPE =
            [jobs]
              var = 42
              [[job1]]
                .extends = job
              [[job]]
        """)
        writetemplate(self.exp_id, 'job', u"""
            %{JOB.var} # %{JOB.id}
        """)
        expecteds = {
            'job':  align(u"""
                        42 # job
                    """),
            'job1': align(u"""
                        42 # job1
                    """),
        }
        ignore = output(script("mkexp "+self.exp_id+".config"))
        for i in ('job', 'job1'):
            result = readfile(join("test", "experiments", self.exp_id,
                                   self.exp_id+"."+i))
            result = align(result)
            self.assertMultiLineEqual(expecteds[i], result)


class JobSiblingsTestCase(MkexpSimpleTestCase):

    def test_sibling_lookup(self):
        writeconfig(self.exp_id, u"""
            EXP_TYPE =
            [jobs]
              [[job1]]
                seniority = elder
              [[job2]]
                seniority = younger
        """)
        writetemplate(self.exp_id, 'job1', u"""
            %{JOB.id}: %{JOB.seniority}
            %{jobs.job1.id}: %{jobs.job1.seniority}
            %{jobs.job2.id}: %{jobs.job2.seniority}
        """)
        writetemplate(self.exp_id, 'job2', u"""
            %{JOB.id}: %{JOB.seniority}
            %{jobs.job2.id}: %{jobs.job2.seniority}
            %{jobs.job1.id}: %{jobs.job1.seniority}
        """)
        expected = align(u"""
            job1: elder
            job1: elder
            job2: younger
            job2: younger
            job2: younger
            job1: elder
        """)
        ignore = output(script("mkexp "+self.exp_id+".config"))
        result = readfile(join("test", "experiments", self.exp_id,
                               self.exp_id+".job1"))
        result += readfile(join("test", "experiments", self.exp_id,
                                self.exp_id+".job2"))
        result = align(result)
        self.assertMultiLineEqual(expected, result)

class NativeVariableTestCase(MkexpSimpleTestCase):

    def test_var_statement(self):
        exp_id = "test_var_statement"
        job_id = "job"
        writeconfig(exp_id, u"""
            EXP_TYPE = 
            GLOBAL1 = 123$${VAR1}456
            GLOBAL2 = $${VAR2}$${VAR3}
            GLOBAL3 = 1, $${VAR2}, 3
            GLOBAL${FOUR} = 4
            [namelists]
              [[namelist]]
                [[[group]]]
                  key = abc$${var}def
            [jobs]
              [["""+job_id+u"""]]
                .var_format = <<<%s>>>
        """)
        writetemplate(exp_id, job_id, u"""
            GLOBAL1=%{GLOBAL1}
            GLOBAL2=%{GLOBAL2}
            GLOBAL3='%{GLOBAL3|join(" ")}'
            GLOBAL4=%{context("GLOBAL<<<FOUR>>>")}
            %{NAMELIST}
        """)
        expected = align(u"""
            GLOBAL1=123<<<VAR1>>>456
            GLOBAL2=<<<VAR2>>><<<VAR3>>>
            GLOBAL3='1 <<<VAR2>>> 3'
            GLOBAL4=4
            &group
                key = 'abc<<<var>>>def'
            /
        """)
        ignore = output(script("mkexp "+exp_id+".config"))
        result = readfile(join("test", "experiments", exp_id, exp_id+"."+job_id))
        result = align(result)
        self.assertMultiLineEqual(expected, result)

    def test_var_separation(self):
        exp_id = "test_var_separation"
        writeconfig(exp_id, u"""
            EXP_TYPE = 
            [jobs]
              use_native_var = $${native_var}
              [[job1]]
                .var_format = <<<%s>>>
              [[job2]]
        """)
        writetemplate(exp_id, "job1", u"""
            %{JOB.use_native_var}
        """)
        writetemplate(exp_id, "job2", u"""
            %{JOB.use_native_var}
        """)
        expected1 = align(u"""
            <<<native_var>>>
        """)
        expected2 = align(u"""
            ${native_var}
        """)
        ignore = output(script("mkexp "+exp_id+".config"))
        result1 = readfile(join("test", "experiments", exp_id, exp_id+".job1"))
        result2 = readfile(join("test", "experiments", exp_id, exp_id+".job2"))
        result1 = align(result1)
        result2 = align(result2)
        self.assertMultiLineEqual(expected1, result1)
        self.assertMultiLineEqual(expected2, result2)

class UnicodeTestCase(MkexpSimpleTestCase):

    def test_value(self):
        self.run_test(u"""
            %{VAR}
        """, u"""
            ÄÖÜäöüß😉
        """, u"""
            VAR = ÄÖÜäöüß😉
        """)

class ExpConfigTestCase(unittest.TestCase):

    def setUp(self):
        text = u"""
            var1 = value1
            var2 = $var1
            [group]
              var1 = value2
              var3 = $var2
        """
        self.config = expconfig.ConfigObj(text.split('\n'))

    def test_default_interpolation_false(self):
        self.assertEqual(self.config['var2'], '$var1')

    def test_template_interpolation(self):
        default_interpolation = self.config.interpolation
        self.config.interpolation = 'template'
        self.assertEqual(self.config['var2'], 'value1')
        self.assertEqual(self.config['group']['var3'], 'value2') # sic!
        self.config.interpolation = default_interpolation

    def test_odict(self):
        odict = expconfig.odict(self.config)
        self.assertEqual(self.config.dict(), odict)

    def test_alias_by_self_merge(self):
        config = expconfig.ConfigObj(self.config)
        config['merged'] = {}
        config['merged'].merge(config)
        # Change in merged group also changes origin
        config['merged']['group']['var1'] = 'changed'
        self.assertEqual(config['group']['var1'], 'changed')
        # Change in merged scalar does not affect origin
        config['merged']['var1'] = 'changed'
        self.assertEqual(config['var1'], 'value1')

    def test_unalias_by_odict_merge(self):
        config = expconfig.ConfigObj(self.config)
        config['merged'] = {}
        config['merged'].merge(expconfig.odict(config))
        # Change in merged group does not affect origin
        config['merged']['group']['var1'] = 'changed'
        self.assertEqual(config['group']['var1'], 'value2')

if __name__ == '__main__':
    unittest.main()

'''
Generate an earth system model configuration from the given configuration file.

$Id$
'''

import collections
import io
import locale
import os
import re
import time # for 'eval' context only

from itertools import dropwhile

import configobj
from configobj import InterpolationError

import feedback

# Utilities used for config and templates

_preferred_encoding = locale.getpreferredencoding()
    
# - Check a namelist logical
def is_set(s):
    if not s:
        return False
    return s.strip('.').lower().startswith('t')

# - Provide merging of comments
def is_not_empty(arg):
    if arg is None:
        return None
    elif not isinstance(arg, (list, tuple)):
        return arg.rstrip()
    else:
        return [_f for _f in [x.rstrip() for x in arg] if _f]

def odict(self):
    '''Return a deepcopy of self as an ordered dictionary.

    >>> n = odict(a)
    >>> n == a
    1
    >>> n is a
    0
    '''
    newdict = collections.OrderedDict()
    for entry in self:
        this_entry = self[entry]
        if isinstance(this_entry, configobj.Section):
            this_entry = odict(this_entry)
        elif isinstance(this_entry, list):
            # create a copy rather than a reference
            this_entry = list(this_entry)
        elif isinstance(this_entry, tuple):
            # create a copy rather than a reference
            this_entry = tuple(this_entry)
        newdict[entry] = this_entry
    return newdict

def merge_comments(this, indict):
    '''Merge comments from indict into current configuration. 

    Requires indict to be merged before being called.
    '''
    if isinstance(indict, configobj.ConfigObj):
        if is_not_empty(indict.initial_comment):
            this.initial_comment = indict.initial_comment
        if is_not_empty(indict.final_comment):
            this.final_comment = indict.final_comment

    for key in indict.scalars:
        if is_not_empty(indict.comments[key]):
            this.comments[key] = indict.comments[key]
        if is_not_empty(indict.inline_comments[key]):
            this.inline_comments[key] = indict.inline_comments[key]
    
    for key in indict.sections:
        merge_comments(this[key], indict[key])

def section_key_info(key):
    key_info = key.split(' ', 1)
    if len(key_info) == 1:
        key_info.append('')
    return key_info

class ConfigObj(configobj.ConfigObj):

    def __init__(self, *args, **kwargs):
        default_args = {
            'encoding': _preferred_encoding,
            'default_encoding': _preferred_encoding,
            'interpolation': False,
        }
        for kw in default_args:
            if not kw in kwargs:
                kwargs[kw] = default_args[kw]
        configobj.ConfigObj.__init__(self, *args, **kwargs)

    def merge(self, indict):
        configobj.ConfigObj.merge(self, indict)
        if isinstance(indict, configobj.Section):
            merge_comments(self, indict)


class ExpConfigError(InterpolationError):
    def __init__(self, message, key):
        message = message.rstrip('.!')
        InterpolationError.__init__(self, 
            "{0} while reading key '{1}'".format(message, key))

class ExpConfig(ConfigObj):
    '''Read and store configuration info from input and experiments' library

    Store environment as default for control settings, then add config from files
    '''
    
    #
    # Basic settings
    #

    exp_lib_dir = 'standard_experiments'
    env_lib_dir = 'standard_environments'
    opt_lib_dir = 'standard_options'
    default_name = 'DEFAULT'
    id_name = 'EXP_ID'
    setup_config_name = 'SETUP.config'

    # Class constructor

    def __init__(self, experiment_config_name,
                 extra_dict={}, config_roots=[''], getexp=False):
        '''Read experiment config to get basic settings
        
        TODO: probably nicer if default experiment is given as argument
        '''

        # State variables
        self.version_info_missing = False

        #
        # Helper functions
        #

        def split_shared_sections(config):
            '''Process sections to expand shared entries as [[job1, job2]]

               Supports jobs, namelists, and job-specific namelists'''

            sep = re.compile(r'\s*,\s*')

            def split_key(current):
                for subkey_orig in current.sections:
                    subkeys = re.split(sep, subkey_orig)
                    if len(subkeys) > 1:
                        subconfig = current[subkey_orig]
                        for subkey in subkeys:
                            if subkey in current:
                                current[subkey].merge(subconfig.dict())
                            else:
                                current[subkey] = subconfig.dict()
                        del current[subkey_orig]

            if 'namelists' in config.sections:
                split_key(config['namelists'])
            if 'jobs' in config.sections:
                split_key(config['jobs'])
                for job in config['jobs'].sections:
                    if 'namelists' in config['jobs'][job].sections:
                        split_key(config['jobs'][job]['namelists'])

        def get_config_name(lib_name, base_name):
            '''Cycle through config path until a match is found.
               
               Return simple path otherwise'''
            config_name = os.path.join(lib_name, base_name)
            for config_root in config_roots:
                tentative_name = os.path.join(config_root, config_name)
                if os.path.exists(tentative_name):
                    config_name = tentative_name
                    break
            return config_name

        # Incremental list assignments (+=, -=)

        def add_to_list(section, key, base_key):
            if isinstance(section[key], (list, tuple)):
                section[base_key].extend(section[key])
            else:
                section[base_key].append(section[key])

        def remove_from_list(section, key, base_key):
            values = section[key]
            if not isinstance(values, (list, tuple)):
                values = [values]
            for value in values:
                try:
                    section[base_key].remove(value)
                except ValueError:
                    pass

        rename_spec_re = re.compile(r'(.*?)\s*>\s*(.*)$')

        def rename_list_items(section, key, base_key):
            values = section[key]
            if not isinstance(values, (list, tuple)):
                values = [values]
            for value in values:
                matches_rename_spec = re.match(rename_spec_re, value)
                if not matches_rename_spec:
                    raise ExpConfigError(
                        "invalid rename '{0}'".format(value), key)
                (old, new) = matches_rename_spec.groups()
                section[base_key] = [v.replace(old, new)
                                     for v in section[base_key]]
                

        list_assign_re = re.compile(r'(.*?)\s*([-+>])$')
        list_assign_op = {
            '+': add_to_list,
            '-': remove_from_list,
            '>': rename_list_items,
        }

        def list_assign(section):

            del_keys = []
            for key in section.scalars:
                is_list_assign = re.match(list_assign_re, key)
                if is_list_assign:
                    base_key = is_list_assign.group(1)
                    list_op = is_list_assign.group(2)

                    if not base_key in section:
                        section[base_key] = []
                    if not isinstance(section[base_key], (list, tuple)):
                        section[base_key] = [section[base_key]]

                    section.main.interpolation = 'template'
                    base_value = section[base_key]
                    section.main.interpolation = False
                    section[base_key] = base_value

                    list_assign_op[list_op](section, key, base_key)

                    del_keys.append(key)
            for key in del_keys:
                del section[key]

            for key in section.sections:
                list_assign(section[key])

        # Helper functions for value definitions

        def read_value(value):
            if os.path.exists(value):
                stream = open(value)
                result = stream.read().strip()
                stream.close()
            else:
                result = ''
            return result

        def sec2time(seconds):
            '''Create time string (HH:MM:SS) from second of day'''
            seconds = int(seconds)
            if seconds >= 86400:
                raise ValueError("invalid second of day '{0}'".format(seconds))
            minutes, s = divmod(seconds, 60)
            h, m = divmod(minutes, 60)
            return "{0:02}:{1:02}:{2:02}".format(h, m, s)

        def split_date(value):
            '''Re-format datetime string to list for use in namelists'''
            match = re.match(r'^0*(\d+)-0*(\d+)-0*(\d+)'
                             r'([T ]0*(\d+)(:0*(\d+)(:0*(\d+))?)?)?$', value)
            if match:
                return [match.groups('0')[i] for i in [0,1,2,4,6,8]]

            match = re.match(r'^0*(\d+?)(\d{2})(\d{2})'
                             r'([T ]0*(\d+)(:0*(\d+)(:0*(\d+))?)?)?$', value)
            if match:
                return [match.groups('0')[i] for i in [0,1,2,4,6,8]]
                
            raise ValueError("invalid date/time '{0}'".format(value))

        def add_years(value, years):
            '''Add specified number of years (possible negative) to date'''
            years = int(years)
            dt = list(map(int, split_date(value)))
            dt[0] += years
            return "{0:+05}-{1:02}-{2:02}".format(*dt).lstrip('+')

        def add_days(value, days):
            '''Add specified number of days (possible negative) to date'''
            def leap(year):
                return (not year % 4) and (not (not year % 100) or (not year % 400)) 
            def monlen(year, mon):
                monlens = (0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 0)
                return monlens[mon] + (mon == 2 and leap(year))
            def add_days_(year, mon, day, days):
                while True:
                    if mon == 0:
                        year -= 1
                        mon = 12
                        day = monlen(year, 12)
                        continue
                    if mon == 13:
                        year += 1
                        mon = 1
                        day = 1
                        continue
                    if day + days <= 0:
                        days += day
                        mon -= 1
                        day = monlen(year, mon)
                        continue
                    if day + days > monlen(year, mon):
                        days -= monlen(year, mon) - day + 1
                        mon += 1
                        day = 1
                        continue
                    day += days
                    break

                return (year, mon, day)

            days = int(days)
            dt = list(map(int, split_date(value)))
            dt = add_days_(dt[0], dt[1], dt[2], days)
            return "{0:+05}-{1:02}-{2:02}".format(*dt).lstrip('+')

        def eval_value(value):
            '''
                Evaluate key as python expression,
                return as string or sequence of strings.
            '''
            result = eval(value)
            if isinstance(result, (list, tuple)):
                result = list(map(str, result))
            else:
                result = str(result)
            return result

        def eval_value_string(value):
            '''
                Evaluate key as python expression,
                return as string or sequence of strings.
            '''
            result = eval_value(value)
            if isinstance(result, (list, tuple)):
                result = ", ".join(result)
            return result

        def eval_expression(value):
            '''
                Check if value is a supported expression.
                If so, evaluate and return result, otherwise just pass through.
            '''
            match = re.match(r'^eval\((.*)\)$', value, re.S)
            if match:
                return eval_value(match.group(1))

            match = re.match(r'^evals\((.*)\)$', value, re.S)
            if match:
                return eval_value_string(match.group(1))

            match = re.match(r'^add_(years|days)\(\s*([-\d]+([T ][\d:]+)?)\s*,\s*([-+]?\d+)\s*\)$', value, re.S)
            if match:
                if match.group(1) == 'days':
                    return add_days(match.group(2), match.group(4))
                return add_years(match.group(2), match.group(4))

            match = re.match(r'^split_date\((.*)\)$', value, re.S)
            if match:
                return split_date(match.group(1))

            match = re.match(r'^sec2time\((.*)\)$', value, re.S)
            if match:
                return sec2time(match.group(1))

            match = re.match(r'^read\((.*)\)$', value, re.S)
            if match:
                return read_value(match.group(1))

            return value

        # Interpolate and evaluate keys if they are an expression
        def eval_key(section, key):
            try:
                value = section[key]
                if isinstance(value, (list, tuple)):
                    value = list(map(eval_expression, value))
                elif not isinstance(value, dict):
                    value = eval_expression(value)
                if isinstance(value, (list, tuple)):
                    value = [v.replace('$', '$$') for v in value]
                elif not isinstance(value, dict):
                    value = value.replace('$', '$$')
            except (InterpolationError, ValueError) as error:
                raise ExpConfigError(str(error), key)
            section[key] = value

        # Undo remaining changes from walk with eval_key
        def uneval_key(section, key):
            try:
                value = section[key]
                if isinstance(value, (list, tuple)):
                    value = [v.replace('$$', '$') for v in value]
                elif not isinstance(value, dict):
                    value = value.replace('$$', '$')
            except (InterpolationError, ValueError) as error:
                raise ExpConfigError(str(error), key)
            section[key] = value

        # Move version info from local config to global list
        def register_version(pre_config, config_versions):
            if 'VERSION_' in pre_config:
                config_versions.append(pre_config['VERSION_'])
                del pre_config['VERSION_']
            else:
                self.version_info_missing = True

        #
        # Method body
        #

        # Pre-read basic experiment settings

        pre_config = None
        setup_config_name = get_config_name('', ExpConfig.setup_config_name)
        if os.path.exists(setup_config_name):
            pre_config = ConfigObj(setup_config_name)
        user_config = ConfigObj(experiment_config_name)
        if pre_config:
            pre_config.merge(user_config)
        else:
            pre_config = user_config

        experiment_type = extra_dict.get('EXP_TYPE', pre_config['EXP_TYPE'])
        # Empty environment should load default
        environment = extra_dict.get('ENVIRONMENT', 
                      pre_config.get('ENVIRONMENT',
                      ExpConfig.default_name))
        # Options should always be treated as a list
        setup_options = extra_dict.get('SETUP_OPTIONS',
                        pre_config.get('SETUP_OPTIONS',
                        ''))
        if not isinstance(setup_options, (list, tuple)):
            if setup_options:
                setup_options = [setup_options]
            else:
                setup_options = []
        exp_options = extra_dict.get('EXP_OPTIONS',
                      pre_config.get('EXP_OPTIONS',
                      ''))
        if not isinstance(exp_options, (list, tuple)):
            if exp_options:
                exp_options = [exp_options]
            else:
                exp_options = []
        options = setup_options + exp_options
        # Backwards compatibility ENVIRONMENT -> QUEUE_TYPE
        if environment == ExpConfig.default_name and 'QUEUE_TYPE' in pre_config:
            feedback.warning("found obsolete keyword 'QUEUE_TYPE'; "
                             "should be replaced by 'ENVIRONMENT'")
            environment = pre_config['QUEUE_TYPE']
        # Load default if environment was deliberately set to empty
        if not environment:
            environment = ExpConfig.default_name

        pre_config = None
        user_config = None

        # Start from empty configuration

        pre_config = ConfigObj()
        config_versions = []

        # Get default experiment id from file name
        pre_config[ExpConfig.id_name] = os.path.splitext(
            os.path.basename(experiment_config_name)
        )[0]

        # Read Environment

        env_dict = dict(os.environ)
        if not getexp:
            # Mask literal dollar characters
            for key, value in env_dict.items():
                env_dict[key] = value.replace('$', '$$')
        pre_config.merge({'DEFAULT': {}})
        for key, value in sorted(env_dict.items()):
            pre_config['DEFAULT'][key] = value

        # Read experiment settings from library (default and type specific)

        lib_config_name = get_config_name(ExpConfig.exp_lib_dir,
                                          ExpConfig.default_name+'.config')
        pre_config.merge(ConfigObj(lib_config_name))
        split_shared_sections(pre_config)
        register_version(pre_config, config_versions)

        if os.path.exists(setup_config_name):
            pre_config.merge(ConfigObj(setup_config_name))
            split_shared_sections(pre_config)
            list_assign(pre_config)
            register_version(pre_config, config_versions)

        lib_config_name = get_config_name(ExpConfig.exp_lib_dir, 
                                          experiment_type+'.config')
        if os.path.exists(lib_config_name):
            pre_config.merge(ConfigObj(lib_config_name))
            split_shared_sections(pre_config)
            list_assign(pre_config)
            register_version(pre_config, config_versions)
        else:
            feedback.warning("cannot find experiment config for '%s', "+
                             "using default only", experiment_type)

        for option in options:
            lib_config_name = get_config_name(ExpConfig.opt_lib_dir, 
                                              option+'.config')
            if os.path.exists(lib_config_name):
                pre_config.merge(ConfigObj(lib_config_name))
                split_shared_sections(pre_config)
                list_assign(pre_config)
                register_version(pre_config, config_versions)
            else:
                feedback.warning("cannot find config for option '%s', using "+
                                 "default/experiment type only", option)

        # Read host environment settings from library

        lib_config_name = get_config_name(ExpConfig.env_lib_dir,
                                          environment+'.config')

        if os.path.exists(lib_config_name):
            pre_config.merge(ConfigObj(lib_config_name))
            list_assign(pre_config)
            register_version(pre_config, config_versions)

        # Warn user if at least one config had no version info
        if self.version_info_missing:
            feedback.info("version info for standard config is incomplete")

        # Re-read config to allow overriding default settings
        # TODO: probably nicer if default experiment is given as argument
        experiment_config = ConfigObj(experiment_config_name)
        pre_config.merge(experiment_config)
        split_shared_sections(pre_config)
        list_assign(pre_config)

        # Add extra dictionary
        pre_config.merge(extra_dict)

        # Backwards compatibility ENVIRONMENT -> QUEUE_TYPE
        pre_config['ENVIRONMENT'] = environment

        # Add complete versioning info
        if not getexp:
            pre_config['VERSIONS_'] = config_versions

        # Extract experiment description from initial comment
        # if not set explicitly
        if 'EXP_DESCRIPTION' not in pre_config:
            is_empty = lambda s: re.match(r'^[\s#]*$', s)
            rm_comment = lambda s: re.sub(r'^\s*# ?', '', s)       
            pre_config['EXP_DESCRIPTION'] = "\n".join(
                reversed(list(
                    dropwhile(is_empty,
                        reversed(list(
                            dropwhile(is_empty,
                                list(map(rm_comment,
                                    experiment_config.initial_comment))
                            )
                        )) 
                    )
                ))
            )

        # Additional experiment info
        self.experiment_id = pre_config[ExpConfig.id_name]
        self.experiment_kind = re.sub(r'-\w+$', '', experiment_type)

        # Re-read merged config with interpolation set.
        # This works around incomprehensible inheritance of interpolation with
        # merge. Make sure that all values are interpolated

        config_lines = io.BytesIO()

        pre_config.write(config_lines)
        pre_config = None

        config_lines.seek(0)
        pre_config = ConfigObj(io.TextIOWrapper(config_lines),
                               interpolation=False if getexp else 'template')

        pre_config.walk(eval_key)

        # Re-read final config without interpolation.
        # This allows copying data without evaluation of version keywords.

        config_lines = io.BytesIO()
        
        pre_config.write(config_lines)
        pre_config = None

        config_lines.seek(0)
        ConfigObj.__init__(self, io.TextIOWrapper(config_lines))
        self.walk(uneval_key)


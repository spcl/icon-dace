'''
Utilities for handling 'files' section
'''

import os
import re

def _get_dir(section):
    '''
    Expand input files section into directory name
    '''

    # Create section list
    parent = section
    ancestry = [parent]
    while not parent.parent is parent:
        ancestry.append(parent.parent)
        parent = parent.parent
    ancestry.reverse()

    # Extract directory name   
    dir_name = ''
    for parent in ancestry:
        if '.base_dir' in parent:
            dir_name = parent['.base_dir']
        elif '.sub_dir' in parent:
            dir_name = os.path.join(dir_name, parent['.sub_dir'])

    return dir_name

def get_dir(section):
    return re.sub(r'\$\{(\w+)\}', lambda m: section.main[m.group(1)], _get_dir(section))


def get_file(section, name):
    '''
    Expand input files into full file names
    '''
   
    input_file = section.get(name, name)

    # Check input file for absolute path
    if not os.path.isabs(input_file):
        input_file = os.path.join(_get_dir(section), input_file)

    var_format = section.main['JOB'].get('.var_format', '${%s}')
    result = re.sub(r'\$\{(\w+)\}',
                    lambda m: section.main.get(m.group(1), var_format%m.group(1)),
                    input_file)

    return result

'''
Specialized argument parser for mkexp tool.

$Id$
'''

import argparse

import package_info

class MkexpArgumentParser(argparse.ArgumentParser):
    ''' Specialized argument parser for mkexp tool
    
    Also used for parsing the 'update' script provided with the run scripts
    '''
    def __init__(self):
        argparse.ArgumentParser.__init__(self, description=
            'Generate an experiment from templates and the given configuration file.')
        self.add_argument('config', help='experiment configuration file name')
        self.add_argument('assigns', metavar='key=value', nargs='*',
                           help='override configuration file settings')
        self.add_argument('-V', '--version', action='version',
                          version=package_info.version)
        self.add_argument('-p', '--path', 
                          help='search path for default config and templates')
        self.add_argument('-m', '--no-make-dirs',
                          action='store_false', dest='make_dirs',
                          help='do not create work and data directories')
        self.add_argument('-q', '--quiet',
                          action='store_true', dest='quiet',
                          help='suppress informative messages')
        self.add_argument('-g', '--getexp', action='store_true',
                          help='load flat config (from getexp -vv)') 

def get_key_chain(key):
    sections = []
    pos = key.rfind('.')
    while(pos >= 0):
        if key[pos-1] == '.':
            key = key[0:pos-1]+key[pos:]
            pos = key.rfind('.', 0, pos-1)
        else:
            sections.append(key[pos+1:])
            key = key[0:pos]
            pos = key.rfind('.')
    if key:
        sections.append(key)
    return sections

def assigns_to_dicts(args):
    
    def value_to_list(value):
        result = value.split(',')
        if len(result) == 1:
            result = result[0]
        return result

    def assign_to_dict(assign):
        (key, value) = assign.split('=', 1)
        chain = get_key_chain(key)
        current = value_to_list(value)
        for key in chain:
            current = {key: current}
        return current

    return list(map(assign_to_dict, args.assigns))


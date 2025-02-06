
'''
Frontend for using logging module as terminal feedback facility.

$Id$
'''

import logging
import sys

logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.WARNING)
logging.addLevelName(logging.DEBUG, 'Debug')
logging.addLevelName(logging.INFO, 'Note')
logging.addLevelName(logging.WARNING, 'Hey')
logging.addLevelName(logging.ERROR, 'Oops')
logging.addLevelName(logging.CRITICAL, 'Sorry')

DEBUG = logging.DEBUG
INFO = logging.INFO
WARNING = logging.WARNING
ERROR = logging.ERROR
CRITICAL = logging.CRITICAL

debug = logging.debug
info = logging.info
warn = warning = logging.warning
error = logging.error
critical = logging.critical

def die(message, *args, **kwargs):
    status = kwargs.pop('status', 1)
    if kwargs:
        (key, dont_care) = kwargs.popitem()
        raise TypeError("die() got an unexpected keyword argument '%s'"%key)
    error(message, *args)
    sys.exit(status)

def setLevel(level):
    logging.getLogger().setLevel(level)

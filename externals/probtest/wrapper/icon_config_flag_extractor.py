import re
import sys

file_path = sys.argv[1]
builder = sys.argv[2]

regex = '.*{0}.*configureflags="(?P<flags>.*?)"'.format(builder)

with open(file_path) as f:
    data = f.read()
    try:
        m = re.search(regex, data, re.IGNORECASE).groupdict()
        print(m["flags"])
    except AttributeError:
        print("ERROR: did not find configure flags for {}".format(builder))

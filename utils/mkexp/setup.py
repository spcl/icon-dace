from setuptools import setup
import package_info

setup(
    name = package_info.name,
    version = package_info.version,
    description = 'Run-script generation for earth system models',
    long_description = open('README.txt').read(),
    author = 'Karl-Hermann Wieners',
    author_email = 'karl-hermann.wieners@mpimet.mpg.de',
    url = 'http://code.mpimet.mpg.de/projects/esmenv',
    py_modules = ['_version', 'configobj', 'validate', 'feedback', 'expargparse', 'expconfig', 'files', 'package_info', 'update'],
    scripts = ['mkexp', 'getexp', 'rmexp', 'diffexp', 'diffpath', 'cpexp', 'cppath', 'duexp', 'getconfig', 'editexp', 'upexp', 'setconfig', 'selconfig', 'namelist2config', 'files2config', 'importexp', 'unmergeconfig', 'compconfig', 'diffconfig'],
    data_files = [('share/doc/'+package_info.name, ['doc/mkexp.pdf', 'mkexp.bash'])],
    platforms = ['Posix'],
    license = 'LICENSE.txt',
    requires = ['Jinja2(>= 2.6)']
)

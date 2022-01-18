#!/usr/bin/env python3

import os
from distutils.core import setup
from distutils.command.install import INSTALL_SCHEMES

def fullsplit(path, result=None):
    """
    Split a pathname into components (the opposite of os.path.join) in a
    platform-neutral way.
    """
    if result is None:
        result = []
    head, tail = os.path.split(path)
    if head == '':
        return [tail] + result
    if head == path:
        return result
    return fullsplit(head, [tail] + result)

# Tell distutils not to put the data_files in platform-specific installation
# locations. See here for an explanation:
# http://groups.google.com/group/comp.lang.python/browse_thread/thread/35ec7b2fed36eaec/2105ee4d9e8042cb
for scheme in list(INSTALL_SCHEMES.values()):
    scheme['data'] = scheme['purelib']

packages = []
data_files = []
root_dir = os.path.dirname(__file__)
if root_dir != '':
    os.chdir(root_dir)

for dirpath, dirnames, filenames in os.walk('Cattery'):
    # Ignore dirnames that start with '.'
    dirnames[:] = [d for d in dirnames if not d.startswith('.') and d != '__pycache__']
    if '__init__.py' in filenames:
        packages.append('.'.join(fullsplit(dirpath)))
    elif filenames:
        data_files.append([dirpath, [os.path.join(dirpath, f) for f in filenames]])

install_requires = [
    'numpy>=1.16',
    'purr',
    'astropy>=3.0.0',
    'python_casacore',
    'scipy',
    'astro_kittens',
    'astro_pyxis',
    'six'
    # 'Timba' is not on pypi
]


setup(name='meqtrees_cattery',
      version='1.7.4',
      python_requires='>=3.0.0',
      description='MeqTrees-based frameworks for simulation and calibration of radio interferometers ',
      author='Oleg Smirnov',
      author_email='osmirnov@gmail.com',
      url='https://github.com/ska-sa/meqtrees-cattery',
      packages=packages,
      data_files=data_files,
      install_requires=install_requires
     )

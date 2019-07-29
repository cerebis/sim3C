import re
from setuptools import setup, find_packages

with open('README.md', 'r') as fh:
    long_description = fh.read()

version_str = None
VERSION_FILE = "sim3C/_version.py"
with open(VERSION_FILE, "rt") as vh:
    for _line in vh:
        mo = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", _line, re.M)
        if mo:
            version_str = mo.group(1)
            break

if version_str is None:
    raise RuntimeError("Unable to find version string in {}".format(VERSION_FILE))

setup(
    name='sim3C',
    description='Hi-C read-pair simulator',
    long_description=long_description,
    version=version_str,
    author='Matthew Z DeMaere',
    author_email='matt.demaere@gmail.com',
    platforms='Linux-86_x64',
    packages=find_packages(),
    url='https://github.com/cerebis/sim3C',
    license='GNU General Public License v3',
    include_package_data=True,
    zip_safe=False,

    install_requires=[
        'biopython',
        'intervaltree',
        'numba',
        'numpy',
        'scipy',
        'pyyaml',
        'tqdm'
    ],

    classifiers=[
        'Programming Language :: Python :: 3.6',
        'License :: OSI Approved :: GNU General Public License v3',
        'Operating System :: POSIX :: Linux',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Development Status :: 5 - Alpha'
    ],

    entry_points={
        'console_scripts': ['sim3C=sim3C.command_line:main'],
    }
)

import re
from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize
import numpy

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


extension = Extension("sim3C.faster", ["sim3C/faster.pyx"],
                      include_dirs=['venv/include', numpy.get_include()],
                      libraries=['pcg_random'],
                      library_dirs=['/Users/cerebis/git/sim3C/py3/venv/lib'],
                      extra_compile_args=['-std=c99'],
                      )

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

    ext_modules=cythonize(
        [extension],
        compiler_directives={'language_level': "3"},
        annotate=True,
        #force=True
    ),

    package_data = { 'sim3C': ['faster.pyx']},

    install_requires=[
        'biopython~=1.81',
        'intervaltree',
        'numba~=0.58',
        'numpy~=1.26',
        'scipy~=1.11',
        'dnaio',
        'pyyaml',
        'tqdm'
    ],

    classifiers=[
        'Programming Language :: Python :: 3.11',
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

from Cython.Build import cythonize
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext as build_ext_orig
import numpy
import os
import re
import subprocess


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


class TarballExtension(Extension, object):
    """
    Simple class for use by build_tarball
    """
    def __init__(self, name, url, parent_pkg):
        """
        :param name: a name for the extension. This is used for naming the tarball and extracted parent path
        :param url: the remote location of the tarball (GitHub)
        :param parent_pkg: parent containing package name
        """
        super(TarballExtension, self).__init__(name, sources=[])
        self.url = url
        self.parent_pkg = parent_pkg

        # attempt to use GNU tar and not Mac OSX bsd tar or the like
        self.tar_cmd = 'tar'
        if re.search(r'GNU', subprocess.check_output([self.tar_cmd, '--version']).decode()) is None:
            self.tar_cmd = 'gtar'
        if re.search(r'GNU', subprocess.check_output([self.tar_cmd, '--version']).decode()) is None:
            raise IOError('GNU tar was not found and installation requires special features')

    @property
    def tarball(self):
        """
        :return: a name for the tarball based on the extension name
        """
        return '{}_tarball.tar.gz'.format(self.name)


def curl_exists():
    try:
        subprocess.check_call(['curl', '--version'])
        return True
    except subprocess.CalledProcessError:
        return False


def wget_exists():
    try:
        subprocess.check_call(['wget', '--version'])
        return True
    except subprocess.CalledProcessError:
        return False


class build_tarball_ext(build_ext_orig, object):
    """
    Build a C/C++ Make projects from remote tarballs and place the binaries in proxigenomics_toolkit/external

    Build at install time allows easier support of runtime architectures which vary widely in age, making
    supplying a universal static binaries for external helpers difficult.
    """
    def run(self):
        for ext in self.extensions:
            if isinstance(ext, TarballExtension):
                self.build_tarball(ext)
            else:
                super(build_tarball_ext, self).run()

    def build_tarball(self, ext):
        # fetch the relevant commit from GitHub
        if curl_exists():
            self.spawn(['curl', '-L', ext.url, '-o', ext.tarball])
        elif wget_exists():
            self.spawn(['wget', '-O', ext.tarball, ext.url])
        else:
            raise IOError('Building {} requires either curl or wget be installed'.format(ext.parent_pkg))

        build_dir = os.path.join(self.build_lib, ext.parent_pkg)
        # rename parent folder to something simple and consistent
        self.spawn([ext.tar_cmd, '--transform=s,[^/]*,{},'.format(ext.name), '-xzvf', ext.tarball])
        # build
        self.spawn(['mkdir', '-p', 'sim3C/external/lib', 'sim3C/external/include'])
        self.spawn(['make', '-j4', '-C', ext.name])
        print(f'BUILD DIR = {build_dir} ABS: {os.path.abspath(build_dir)}')
        self.spawn(['make', '-C', ext.name, 'install', 'PREFIX=../sim3C/external'])


cython_extension = Extension("sim3C.faster", ["sim3C/faster.pyx"],
                             include_dirs=['sim3C/external/include', numpy.get_include()],
                             libraries=['pcg_random'],
                             library_dirs=['sim3C/external/lib'],
                             extra_compile_args=['-std=c99']
                             )

clib_extension = TarballExtension('pcg',
                                  'https://github.com/imneme/pcg-c/tarball/master',
                                  'sim3C'
                                  )

setup(
    name='sim3C',
    description='Hi-C read-pair simulator',
    long_description=long_description,
    version=version_str,
    author='Matthew Z DeMaere',
    author_email='matt.demaere@gmail.com',
    platforms='Linux-86_x64',
    packages=find_packages() + ['sim3C.Illumina_profiles'],
    url='https://github.com/cerebis/sim3C',
    license='GNU General Public License v3',
    include_package_data=True,
    zip_safe=False,

    cmdclass={
        'build_ext': build_tarball_ext
    },

    ext_modules=[clib_extension] + cythonize([cython_extension],
                                             compiler_directives={'language_level': "3"},
                                             annotate=True),

    package_data={'sim3C': ['faster.pyx']},

    install_requires=[
        'biopython~=1.81',
        'intervaltree',
        'numba~=0.58',
        'numpy~=1.26',
        'scipy~=1.11',
        'dnaio',
        'pyyaml',
        'tqdm',
        'toml'
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

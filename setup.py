#!/usr/bin/env python
# Copyright (C) 2012 Alex Nitz, Duncan Brown, Andrew Miller, Josh Willis
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

"""
setup.py file for PyCBC package
"""

import sys
import os, subprocess, shutil
import platform

from setuptools import Extension, setup, Command
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools import find_packages


requires = []
setup_requires = ['numpy>=1.16.0']
install_requires = setup_requires + [
    'cython>=0.29',
    'numpy>=1.16.0,!=1.19.0,<2.0.0',
    'scipy>=0.16.0',
    'astropy>=2.0.3,!=4.2.1,!=4.0.5',
    'matplotlib>=1.5.1',
    'mpld3>=0.3',
    'pillow',
    'h5py>=3.0.0,!=3.7.0',
    'jinja2',
    'Mako>=1.0.1',
    'beautifulsoup4>=4.6.0',
    'tqdm',
    'setuptools',
    'gwdatafind',
    'pegasus-wms.api >= 5.0.6',
    'python-ligo-lw >= 1.7.0',
    'ligo-segments',
    'lalsuite!=7.2',
    'lscsoft-glue>=1.59.3',
    'pykerr',
]

def find_files(dirname, relpath=None):
    def find_paths(dirname):
        items = []
        for fname in os.listdir(dirname):
            path = os.path.join(dirname, fname)
            if os.path.isdir(path):
                items += find_paths(path)
            elif not path.endswith(".py") and not path.endswith(".pyc"):
                items.append(path)
        return items
    items = find_paths(dirname)
    if relpath is None:
        relpath = dirname
    return [os.path.relpath(path, relpath) for path in items]

class cbuild_ext(_build_ext):
    def run(self):
        # At this point we can be sure pip has already installed numpy
        import numpy
        numpy_incl = numpy.get_include()

        for ext in self.extensions:
            if (hasattr(ext, 'include_dirs') and
                    numpy_incl not in ext.include_dirs):
                ext.include_dirs.append(numpy_incl)

        _build_ext.run(self)

def get_version_info():
    """Get VCS info and write version info to version.py.
    """
    from pycbc import _version_helper

    class vdummy(object):
        def __getattr__(self, attr):
            return ''

    # If this is a pycbc git repo always populate version information using GIT
    try:
        vinfo = _version_helper.generate_git_version_info()
    except:
        vinfo = vdummy()
        vinfo.version = '2.4.dev3'
        vinfo.release = 'False'

    version_script = f"""# coding: utf-8
# Generated by setup.py for PyCBC on {vinfo.build_date}.

# general info
version = '{vinfo.version}'
date = '{vinfo.date}'
release = '{vinfo.release}'
last_release = '{vinfo.last_release}'

# git info
git_hash = '{vinfo.hash}'
git_branch = '{vinfo.branch}'
git_tag = '{vinfo.tag}'
git_author = '{vinfo.author}'
git_committer = '{vinfo.committer}'
git_status = '{vinfo.status}'
git_builder = '{vinfo.builder}'
git_build_date = '{vinfo.build_date}'
git_verbose_msg = \"\"\"Version: {vinfo.version}
Branch: {vinfo.branch}
Tag: {vinfo.tag}
Id: {vinfo.hash}
Builder: {vinfo.builder}
Build date: {vinfo.build_date}
Repository status is {vinfo.status}\"\"\"

from pycbc._version import *
"""
    with open('pycbc/version.py', 'wb') as f:
        f.write(version_script.encode('utf-8'))

    from pycbc import version
    return version.version

class build_docs(Command):
    user_options = []
    description = "Build the documentation pages"
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        cmd = (
            "cd docs; "
            "cp Makefile.std Makefile; "
            "sphinx-apidoc -o ./ -f -A 'PyCBC dev team' -V '0.1' ../pycbc "
            "&& make html"
        )
        subprocess.check_call(cmd, stderr=subprocess.STDOUT, shell=True)

class build_gh_pages(Command):
    user_options = []
    description = "Build the documentation pages for GitHub"
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        cmd = (
            "mkdir -p _gh-pages/latest "
            "&& touch _gh-pages/.nojekyll "
            "&& cd docs; "
            "cp Makefile.gh_pages Makefile; "
            "sphinx-apidoc -o ./ -f -A 'PyCBC dev team' -V '0.1' ../pycbc "
            "&& make html"
        )
        subprocess.check_call(cmd, stderr=subprocess.STDOUT, shell=True)

cmdclass = {
    'build_docs': build_docs,
    'build_gh_pages': build_gh_pages,
    'build_ext': cbuild_ext
}

extras_require = {
    'cuda': [
        'pycuda>=2015.1',
        'scikit-cuda',
    ],
    'igwn': [
        'ciecplib>=0.7.0',
    ],
}

# do the actual work of building the package
VERSION = get_version_info()

cythonext = ['waveform.spa_tmplt',
             'waveform.utils',
             'types.array',
             'filter.matchedfilter',
             'vetoes.chisq']
ext = []
cython_compile_args = ['-O3', '-w', '-ffast-math',
                       '-ffinite-math-only']

if platform.machine() == 'x86_64':
    cython_compile_args.append('-msse4.2')
cython_link_args = []

# Mac's clang compiler doesn't have openMP support by default. Therefore
# disable openmp builds on MacOSX. Optimization should never really be a
# concern on that OS, and this line can be commented out if needed anyway.
# Mac's also alias gcc and can run into troubles getting libc correctly
if not sys.platform == 'darwin':
    cython_compile_args += ['-fopenmp']
    cython_link_args += ['-fopenmp']
else:
    cython_compile_args += ["-stdlib=libc++"]
    cython_link_args += ["-stdlib=libc++"]

for name in cythonext:
    e = Extension("pycbc.%s_cpu" % name,
                  ["pycbc/%s_cpu.pyx" % name.replace('.', '/')],
                  extra_compile_args=cython_compile_args,
                  extra_link_args=cython_link_args,
                  compiler_directives={'embedsignature': True})
    ext.append(e)

# Not all modules work like this:
e = Extension("pycbc.fft.fftw_pruned_cython",
              ["pycbc/fft/fftw_pruned_cython.pyx"],
              extra_compile_args=cython_compile_args,
              extra_link_args=cython_link_args,
              compiler_directives={'embedsignature': True})
ext.append(e)
e = Extension("pycbc.events.eventmgr_cython",
              ["pycbc/events/eventmgr_cython.pyx"],
              extra_compile_args=cython_compile_args,
              extra_link_args=cython_link_args,
              compiler_directives={'embedsignature': True})
ext.append(e)
e = Extension("pycbc.events.simd_threshold_cython",
              ["pycbc/events/simd_threshold_cython.pyx"],
              language='c++',
              extra_compile_args=cython_compile_args,
              extra_link_args=cython_link_args,
              compiler_directives={'embedsignature': True})
ext.append(e)
e = Extension("pycbc.filter.simd_correlate_cython",
              ["pycbc/filter/simd_correlate_cython.pyx"],
              language='c++',
              extra_compile_args=cython_compile_args,
              extra_link_args=cython_link_args,
              compiler_directives={'embedsignature': True})
ext.append(e)
e = Extension("pycbc.waveform.decompress_cpu_cython",
              ["pycbc/waveform/decompress_cpu_cython.pyx"],
              language='c++',
              extra_compile_args=cython_compile_args,
              extra_link_args=cython_link_args,
              compiler_directives={'embedsignature': True})
ext.append(e)
e = Extension("pycbc.inference.models.relbin_cpu",
              ["pycbc/inference/models/relbin_cpu.pyx"],
              language='c++',
              extra_compile_args=cython_compile_args,
              extra_link_args=cython_link_args,
              compiler_directives={'embedsignature': True})
ext.append(e)


setup(
    name = 'PyCBC',
    version = VERSION,
    description = 'Core library to analyze gravitational-wave data, find signals, and study their parameters.',
    long_description = open('README.md').read(),
    long_description_content_type='text/markdown',
    author = 'The PyCBC team',
    author_email = 'alex.nitz@gmail.org',
    url = 'http://www.pycbc.org/',
    download_url = f'https://github.com/gwastro/pycbc/tarball/v{VERSION}',
    keywords = [
        'ligo',
        'physics',
        'gravity',
        'signal processing',
        'gravitational waves'
    ],
    cmdclass = cmdclass,
    setup_requires = setup_requires,
    extras_require = extras_require,
    install_requires = install_requires,
    scripts  = find_files('bin', relpath='./'),
    packages = find_packages(),
    package_data = {
        'pycbc.workflow': find_files('pycbc/workflow'),
        'pycbc.results': find_files('pycbc/results'),
        'pycbc.neutron_stars': find_files('pycbc/neutron_stars')
    },
    ext_modules = ext,
    python_requires='>=3.10',
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ],
)

import os
import glob
import multiprocessing
from setuptools import setup, find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name='gff3toembl',
    version='1.0.4',
    description='Convert a GFF3 file to EMBL format for submission',
    long_description=read('README.md'),
    packages = find_packages(),
    author='Andrew J. Page',
    author_email='ap13@sanger.ac.uk',
    url='https://github.com/sanger-pathogens/gff3toembl',
    scripts=glob.glob('scripts/*'),
    test_suite='nose.collector',
    tests_require=['nose >= 1.3', 'mock'],
    license='GPLv3',
    classifiers=[
        "License :: OSI Approved :: GNU General Public License (GPLv3)",
        "Programming Language :: Python",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        ],
)

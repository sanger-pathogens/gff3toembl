# GFF3toEMBL
Converts GFF3 files from Prokka into a format suitable for submission to EMBL.

[![Build Status](https://travis-ci.org/sanger-pathogens/gff3toembl.svg?branch=master)](https://travis-ci.org/sanger-pathogens/gff3toembl)   
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/sanger-pathogens/gff3toembl/blob/master/LICENSE)   
[![status](http://joss.theoj.org/papers/9253390f38f4ce6b71674f433fa72afe/status.svg)](http://joss.theoj.org/papers/9253390f38f4ce6b71674f433fa72afe)  
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/recipes/gff3toembl/README.html)  
[![Container ready](https://img.shields.io/badge/container-ready-brightgreen.svg)](https://quay.io/repository/biocontainers/gff3toembl)  
[![Docker Build Status](https://img.shields.io/docker/build/sangerpathogens/gff3toembl.svg)](https://hub.docker.com/r/sangerpathogens/gff3toembl)  
[![Docker Pulls](https://img.shields.io/docker/pulls/sangerpathogens/gff3toembl.svg)](https://hub.docker.com/r/sangerpathogens/gff3toembl)  
[![codecov](https://codecov.io/gh/sanger-pathogens/gff3toembl/branch/master/graph/badge.svg)](https://codecov.io/gh/sanger-pathogens/gff3toembl)

## Contents
  * [Introduction](#introduction)
  * [Installation](#installation)
    * [Required dependencies](#required-dependencies)
    * [Docker](#docker)
    * [From source](#from-source)
    * [Running the tests](#running-the-tests)
  * [Usage](#usage)
    * [Example data](#example-data)
  * [License](#license)
  * [Feedback/Issues](#feedbackissues)
  * [Citation](#citation)
  * [Known Issues](#known-issues)

## Introduction
Submitting annoated genomes to EMBL is a very difficult and time consuming process. This software converts GFF3 files from the most commonly use prokaryote annotation tool Prokka into a format that is suitable for submission to EMBL. It has been used to prepare more than 30% of all annotated genomes in EMBL/GenBank.

__N.B.__ This implements some EMBL specific conventions and is not a generic conversion tool. It is also not a validator, so you need to pass in parameters which are acceptable to EMBL.

## Installation
GFF3toEMBL has the following dependencies:

### Required dependencies
* [Genometools](https://github.com/genometools/genometools/)

There are a number of ways to install GFF3toEMBL and details are provided below. If you encounter an issue when installing GFF3toEMBL please contact your local system administrator. If you encounter a bug please log it [here](https://github.com/sanger-pathogens/gff3toembl/issues) or email us at path-help@sanger.ac.uk.

### Docker
A docker container is provided with all of the dependancies setup and installed. To install the container:

`docker pull sangerpathogens/gff3toembl`

To run the script from within the container on test data (substituting /home/ubuntu/data for your own directory):

`docker run --rm -it  -v /home/ubuntu/data:/data sangerpathogens/gff3toembl  gff3_to_embl --output_filename /data/output_file.embl ABC 123 PRJ1234 ABC  /opt/gff3toembl-1.1.0/gff3toembl/tests/data/single_feature.gff`

### From source
This is for advanced users. The [homebrew recipe](https://raw.githubusercontent.com/andrewjpage/homebrew-science/gff3toembl/gff3toembl.rb), [Dockerfile](Dockerfile) and the [TravisCI install dependancies script](install_dependencies.sh) all contain steps to setup depenancies and install the software so might be worth looking at for hints.

- Install genometools including python bindings
- git clone git@github.com:sanger-pathogens/gff3toembl.git
- python setup.py install

### Running the tests
Run `python setup.py test`

## Usage
```
usage: gff3_to_embl [-h] [--authors AUTHORS] [--title TITLE]
                    [--publication PUBLICATION] [--genome_type GENOME_TYPE]
                    [--classification CLASSIFICATION]
                    [--output_filename OUTPUT_FILENAME]
                    [--locus_tag LOCUS_TAG]
                    [--translation_table TRANSLATION_TABLE]
                    [--chromosome_list CHROMOSOME_LIST] [--version]
                    organism taxonid project_accession description file

Converts prokaryote GFF3 annotations to EMBL for ENA submission. Cite
http://dx.doi.org/10.21105/joss.00080

positional arguments:
  organism              Organism
  taxonid               Taxon id
  project_accession     Accession number for the project
  description           Genus species subspecies strain of organism
  file                  GFF3 filename

optional arguments:
  -h, --help            show this help message and exit
  --authors AUTHORS, -i AUTHORS
                        Authors (in the EMBL RA line style)
  --title TITLE, -m TITLE
                        Title of paper (in the EMBL RT line style)
  --publication PUBLICATION, -p PUBLICATION
                        Publication or journal name (in the EMBL RL line
                        style)
  --genome_type GENOME_TYPE, -g GENOME_TYPE
                        Genome type (linear/circular)
  --classification CLASSIFICATION, -c CLASSIFICATION
                        Classification (PROK/UNC/..)
  --output_filename OUTPUT_FILENAME, -f OUTPUT_FILENAME
                        Output filename
  --locus_tag LOCUS_TAG, -l LOCUS_TAG
                        Overwrite the locus tag in the annotation file
  --translation_table TRANSLATION_TABLE, -n TRANSLATION_TABLE
                        Translation table
  --chromosome_list CHROMOSOME_LIST, -d CHROMOSOME_LIST
                        Create a chromosome list file, and use the supplied
                        name
  --version             show program's version number and exit
```

An example:
```
gff3_to_embl --authors 'John' --title 'Some title' --publication 'Some journal' \
             --genome_type 'circular' --classification 'PROK' \
             --output_filename /tmp/single_feature.embl --translation_table 11 \
             Organism 1234 'My project' 'My description' gff3toembl/tests/data/single_feature.gff
```

### Example data
The directory 'example_data' contains an input GFF file and the output file along with the command.

## License
GFF3toEMBL is free software, licensed under [GPLv3](https://github.com/sanger-pathogens/gff3toembl/blob/master/LICENSE).

## Feedback/Issues
Please report any issues to the [issues page](https://github.com/sanger-pathogens/gff3toembl/issues) or email path-help@sanger.ac.uk.

## Citation
If you use this software please cite:

__GFF3toEMBL: Preparing annotated assemblies for submission to EMBL__   
Andrew J. Page, Sascha Steinbiss, Ben Taylor, Torsten Seemann, Jacqueline A. Keane   
The Journal of Open Source Software, 1 (6) 2016. [doi: 10.21105/joss.00080](http://dx.doi.org/10.21105/joss.00080)

## Known Issues
This doesn't work with some versions of Genometools on Mac OS X; it appears to work with Genometools 1.5.4


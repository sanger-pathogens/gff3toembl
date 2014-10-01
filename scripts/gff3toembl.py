#!/usr/bin/env python
# encoding: utf-8

import sys
sys.path.append(".")
import argparse
import gff3toembl

parser = argparse.ArgumentParser(description='Iteratively detect recombinations')
parser.add_argument('input_file',       help='Input file in GFF3 format')
parser.add_argument('--outputfile',         '-o', help='Name of output file')
parser.add_argument('--verbose',          '-v', action='count', help='Turn on debugging')

runner  = Gff3toEmbl(parser.parse_args())
runner.parse_and_run()

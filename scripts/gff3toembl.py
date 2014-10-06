#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sysgff3toembl/EMBLWriter.py
import argparse              
from gff3toembl.EMBLWriter import EMBLWriter

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Converts prokaryote GFF3 annotations to ' + \
                                                 'EMBL for ENA submission.')
    parser.add_argument('file', metavar='file', type=str, nargs=1, help='GFF3 filename')
    
    parser.add_argument('--organism',           '-o', help='Organism')
    parser.add_argument('--taxon',              '-t', help='Taxon id', type=int)
    parser.add_argument('--project_accession',  '-a', help='Accession number for the project')
    parser.add_argument('--description',        '-d', help='Genus species subspecies strain of organism')
    parser.add_argument('--authors',            '-i', help='Authors', default = 'Pathogen Genomics')
    parser.add_argument('--title',              '-m', help='Title of paper',default = 'Annotated with Prokka')
    parser.add_argument('--publication',        '-p', help='Publication', default = 'Wellcome Trust Sanger Institute')
    parser.add_argument('--genome_type',        '-g', help='Genome type (linear/circular)', default = 'circular')
    parser.add_argument('--classification',     '-c', help='Classification (PROK/UNC/..)',  default = 'PROK')
    parser.add_argument('--submitter_name',     '-s', help='Submitter name',      default = 'Pathogen Informatics')
    parser.add_argument('--submitter_title',    '-b', help='Submitter title',     default = 'Direct submission')
    parser.add_argument('--submitter_location', '-l', help='Submitter location',  default = 'Wellcome Trust Sanger Institute')
    
    args = parser.parse_args()
    emblwriter = EMBLWriter.EMBLWriter()
    emblwriter.parse_and_run(args)

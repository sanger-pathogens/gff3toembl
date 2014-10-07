#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse              
from gff3toembl import EMBLWriter

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
    parser.add_argument('--submitter_location', '-e', help='Submitter location',  default = 'Wellcome Trust Sanger Institute')
    parser.add_argument('--output_filename',    '-f', help='Output filename',     default = 'output.embl')
    parser.add_argument('--locus_tag',          '-l', help='Overwrite the locus tag in the annotation file')
    
    args = parser.parse_args()
    emblwriter = EMBLWriter.EMBLWriter(args.file[0], args.organism, args.taxonid, args.project, args.description, args.authors, args.title,  args.publication, args.genome_type, args.classification, args.submitter_name, args.submitter_title,  args.submitter_location, args.output_filename, args.locus_tag )
    emblwriter.parse_and_run()
    
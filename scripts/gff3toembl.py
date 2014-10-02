#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from collections import defaultdict
from gt import FeatureNode, CommentNode, SequenceNode, RegionNode, \
               CustomVisitor, CustomStream, MetaNode, EOFNode, \
               FeatureNodeIteratorDepthFirst, GTError, GFF3InStream
from gff3toembl import convert

class EMBLConverter(CustomVisitor):

    def __init__(self):
        CustomVisitor.__init__(self)
        self.seqs = {}
        self.feats = defaultdict(lambda: [], {})
        self.regions = []
        self.converter =  convert.Convert()

    def visit_feature_node(self, fn):
        feature_string = self.converter.construct_feature(feature_type = fn.get_type(), start = fn.get_start(), end = fn.get_end(), strand = fn.get_strand(), feature_attributes = fn.attribs)
        if feature_string != '':
          self.feats[fn.get_seqid()].append(feature_string)

    def visit_region_node(self, rn):
        self.regions.append

    def visit_comment_node(self, cn):
        pass  # for now

    def visit_sequence_node(self, sn):
        self.seqs[sn.get_description()] = sn.get_sequence()

class VisitorStream(CustomStream):

    def __init__(self, instream, visitor):
        CustomStream.__init__(self)
        self.instream = instream
        self.visitor = visitor

    def next(self):
        node = self.instream.next_tree()
        if node:
            node.accept(self.visitor)
        return node

def output_seq(seq):
    converter = convert.Convert()
    sequence_string = converter.construct_sequence(seq)
    return sequence_string

def output_source(sequence_length, organism, taxonid):
    converter = convert.Convert()
    converter.source_template(sequence_length,organism, taxonid)
    return converter
    
def create_output(sequences, organism, taxonid, project, description, authors, title, publication, genome_type, classification, submitter_name, submitter_title, submitter_location):
    converter = convert.Convert()
    i = 1
    for seqid in sorted(sequences):
        target = sys.stdout
        target.write(converter.populated_header(len(conv.seqs[seqid]),  project, description, i, authors, title, publication, genome_type, classification, submitter_name, submitter_title, submitter_location ) )
        target.write(output_source(len(conv.seqs[seqid]), organism, taxonid))
        for feat in conv.feats[seqid]:
            target.write(feat)
        target.write(output_seq(conv.seqs[seqid]))
        target.write("//\n")
        i +=1
    
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

    ins = GFF3InStream(args.file[0])
    conv = EMBLConverter()
    vs = VisitorStream(ins, conv)
    try:
        while (vs.next_tree()):
            pass
    except Exception, e:
        print e
        exit(1)
    create_output_file(conv.seqs.keys(), args.organism, args.taxonid, args.project, args.description, args.authors, args.title, args.publication, args.genome_type, args.classification, args.submitter_name, args.submitter_title, args.submitter_location)


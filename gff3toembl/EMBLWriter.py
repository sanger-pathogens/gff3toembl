#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from collections import defaultdict
from gt import GFF3InStream
               
from gff3toembl import convert
from gff3toembl.EMBLConverter import EMBLConverter
from gff3toembl.VisitorStream import VisitorStream
 
class EMBLWriter():

    def __init__(self):
        self.converter =  convert.Convert()

    def output_seq(seq):
        sequence_string = self.converter.construct_sequence(seq)
        return sequence_string

    def output_source(sequence_length, organism, taxonid):
        source_string = self.converter.source_template(sequence_length,organism, taxonid)
        return source_string
    
    def create_output(sequences, organism, taxonid, project, description, authors, title, publication, genome_type, classification, submitter_name, submitter_title, submitter_location):
        i = 1
        for seqid in sorted(sequences):
            target = sys.stdout
            target.write(self.converter.populated_header(len(conv.seqs[seqid]),  project, description, i, authors, title, publication, genome_type, classification, submitter_name, submitter_title, submitter_location ) )
            target.write(self.output_source(len(conv.seqs[seqid]), organism, taxonid))
            for feat in conv.feats[seqid]:
                target.write(feat)
            target.write(self.output_seq(conv.seqs[seqid]))
            target.write("//\n")
            i +=1

    def parse_and_run(args):
        ins = GFF3InStream(args.file[0])
        conv = EMBLConverter()
        vs = VisitorStream(ins, conv)
        try:
            while (vs.next_tree()):
                pass
        except Exception, e:
            print e
            exit(1)
        self.create_output_file(conv.seqs.keys(), args.organism, args.taxonid, args.project, args.description, args.authors, args.title, args.publication, args.genome_type, args.classification, args.submitter_name, args.submitter_title, args.submitter_location)


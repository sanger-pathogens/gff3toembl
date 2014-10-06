#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from collections import defaultdict
from gt import GFF3InStream
               
from gff3toembl import convert
from gff3toembl.EMBLConverter import EMBLConverter
from gff3toembl.VisitorStream import VisitorStream
 
class EMBLWriter():

    def __init__(self, gff3_file, organism, taxonid, project, description, authors, title,  publication, genome_type, classification, submitter_name, submitter_title,  submitter_location, output_filename, locus_tag = None):
        self.locus_tag          = locus_tag
        self.converter          = convert.Convert(locus_tag)
        self.conv               = EMBLConverter()
        self.gff3_file          = gff3_file
        self.organism           = organism          
        self.taxonid            = taxonid           
        self.project            = project           
        self.description        = description       
        self.authors            = authors           
        self.title              = title             
        self.publication        = publication       
        self.genome_type        = genome_type       
        self.classification     = classification    
        self.submitter_name     = submitter_name    
        self.submitter_title    = submitter_title   
        self.submitter_location = submitter_location
        self.output_filename    = output_filename
 
 
    def output_seq(self, seq):
        sequence_string = self.converter.construct_sequence(seq)
        return sequence_string

    def output_source(self, sequence_length, organism, taxonid):
        source_string = self.converter.source_template(sequence_length,organism, taxonid)
        return source_string
    
    def create_output_file(self, sequences, organism, taxonid, project, description, authors, title, publication, genome_type, classification, submitter_name, submitter_title, submitter_location):
        i = 1
        target = open(self.output_filename, 'w')
        for seqid in sorted(sequences):
            target.write(self.converter.populated_header(len(self.conv.seqs[seqid]),  project, description, i, authors, title, publication, genome_type, classification, submitter_name, submitter_title, submitter_location ) )
            target.write(self.output_source(len(self.conv.seqs[seqid]), organism, taxonid))
            for feat in self.conv.feats[seqid]:
                target.write(feat)
            target.write(self.output_seq(self.conv.seqs[seqid]))
            target.write("//\n")
            i +=1
        target.close()

    def parse_and_run(self):
        ins = GFF3InStream(self.gff3_file)  
        vs = VisitorStream(ins, self.conv)
        try:
            while (vs.next_tree()):
                pass
        except Exception, e:
            print e
            exit(1)
        self.create_output_file(self.conv.seqs.keys(), self.organism, self.taxonid, self.project, self.description, self.authors, self.title, self.publication, self.genome_type, self.classification, self.submitter_name, self.submitter_title, self.submitter_location)


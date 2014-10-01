#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from collections import defaultdict
from gt import FeatureNode, CommentNode, SequenceNode, RegionNode, \
               CustomVisitor, CustomStream, MetaNode, EOFNode, \
               FeatureNodeIteratorDepthFirst, GTError, GFF3InStream

header = """\
ID   XXX; XXX; linear; genomic DNA; XXX; PRO; %d BP.
XX
AC   [ACCESSION];
XX
PR   Project:[PROJECT];
XX
DE   [DESCRIPTION]
XX
RN   [1]
RA   [AUTHOR];
RT   [TITLE];
RL   [LOCATION];
XX
CC   Data release policy http://www.sanger.ac.uk/legal/#t_2
XX
FH   Key             Location/Qualifiers
FH\
"""

source_template = """\
FT   source          1..%d
FT                   /organism="[ORGANISM]"
FT                   /mol_type="genomic DNA"
FT                   /db_xref="taxon:XXX"\
"""

class EMBLConverter(CustomVisitor):

    def __init__(self):
        CustomVisitor.__init__(self)
        self.seqs = {}
        self.feats = defaultdict(lambda: [], {})
        self.regions = []

    def visit_feature_node(self, fn):
        string = ""
        cmp1 = ''
        cmp2 = ''
        if fn.get_strand() == '-':
            cmp1 = 'complement('
            cmp2 = ')'
        string += "FT   %s%s%s%d..%d%s\n" % (fn.get_type(), ' ' * (16-len(fn.get_type())), cmp1, fn.get_start(), fn.get_end(), cmp2)
        for attr in fn.attribs.keys():
            string += "FT%s/%s=\"%s\"\n" % (' ' * 19, attr, fn.attribs[attr])
        self.feats[fn.get_seqid()].append(string)

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

def output_seq(header, seq, target):
    seq = seq.lower()
    a = seq.count('a')
    c = seq.count('c')
    g = seq.count('g')
    t = seq.count('t')
    o = len(seq) - a - c - g - t;
    print "SQ   Sequence %d BP; %d A; %d C; %d G; %d T; %d other;" % \
      (len(seq), a, c, g, t, o)
    target.write("     ")
    i = 1
    for j in range(len(seq)):
        target.write(seq[j])
        if (i) % 10 == 0:
            target.write(" ")
        if (i) % 60 == 0:
            target.write("%9s\n     " % (i))
        i += 1
    target.write(' '*(80-i%60-(i%60)/10-13) + "%9d\n" % (i - 1))

def output_source(length, target):
    target.write(source_template % length)
    target.write("\n")
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Converts prokaryote GFF3 annotations to ' + \
                                                 'EMBL for ENA submission.')
    parser.add_argument('file', metavar='file', type=str, nargs=1,
                         help='GFF3 filename')
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

    for seqid in sorted(conv.seqs.keys()):
        target = sys.stdout
        target.write(header % len(conv.seqs[seqid]))
        target.write("\n")
        output_source(len(conv.seqs[seqid]), target)
        for feat in conv.feats[seqid]:
            target.write(feat)
        output_seq(seqid, conv.seqs[seqid], target)
        target.write("//\n\n")

import sys
from gt import CustomVisitor
from collections import defaultdict
from gff3toembl import convert

class EMBLConverter(CustomVisitor):

    def __init__(self, converter):
        CustomVisitor.__init__(self)
        self.seqs = {}
        self.feats = defaultdict(lambda: [], {})
        self.regions = []
        self.converter =  converter

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

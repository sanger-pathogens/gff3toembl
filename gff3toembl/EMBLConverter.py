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
        self.features_seen = {}

    def visit_feature_node(self, feature_node):
        feature_string = self.converter.construct_feature(feature_type = feature_node.get_type(), start = feature_node.get_start(), end = feature_node.get_end(), strand = feature_node.get_strand(), feature_attributes = feature_node.attribs)
        if feature_string != '':
          feature_seq_coords = str(feature_node.get_seqid()) + "_"+ str(feature_node.get_start()) + "_" +str(feature_node.get_end())
          if feature_seq_coords not     in self.features_seen :          
            self.features_seen[feature_seq_coords] = 1
            self.feats[feature_node.get_seqid()].append(feature_string)

    def visit_region_node(self, region_node):
        pass  # for now

    def visit_comment_node(self, comment_node):
        pass  # for now

    def visit_sequence_node(self, sequence_node):
        self.seqs[sequence_node.get_description()] = sequence_node.get_sequence()

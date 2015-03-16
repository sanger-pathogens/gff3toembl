import sys
from gt import CustomVisitor
from collections import defaultdict
from gff3toembl import convert
from gff3toembl.EMBLContig import EMBLFeature

class EMBLConverter(CustomVisitor):

    def __init__(self, locus_tag = None, translation_table = 11):
        CustomVisitor.__init__(self)
        self.seqs = {}
        self.feats = defaultdict(lambda: [], {})
        self.regions = []
        self.locus_tag = locus_tag
        self.translation_table = translation_table
        self.features_seen = []

    def update_features_seen(self, sequence_id, feature_type, feature_start, feature_end):
      unique_feature_reference = "%s_%s_%s_%s" % (sequence_id, feature_type, feature_start, feature_end)
      if unique_feature_reference in self.features_seen:
        return True
      else:
        self.features_seen.append(unique_feature_reference)
        return False

    def visit_feature_node(self, feature_node):
      feature = EMBLFeature(feature_type = feature_node.get_type(), start = feature_node.get_start(),
                            end = feature_node.get_end(), strand = feature_node.get_strand(),
                            feature_attributes = feature_node.attribs,
                            locus_tag = self.locus_tag, translation_table = self.translation_table)
      feature_string = feature.format()
      if feature_string != None:
        if not self.update_features_seen(feature_node.get_seqid(), feature_node.get_type(),
                                         feature_node.get_start(), feature_node.get_end()):
            self.feats[feature_node.get_seqid()].append(feature_string)

    def visit_region_node(self, region_node):
        pass  # for now

    def visit_comment_node(self, comment_node):
        pass  # for now

    def visit_sequence_node(self, sequence_node):
        self.seqs[sequence_node.get_description()] = sequence_node.get_sequence()

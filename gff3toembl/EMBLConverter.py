import sys
from gt import CustomVisitor
from collections import defaultdict
import gff3toembl
from gff3toembl import convert
from gff3toembl.EMBLContig import EMBLContig

class EMBLConverter(CustomVisitor):

    def __init__(self,
                 authors="Pathogen Genomics",
                 classification="UNC",
                 genome_type="circular",
                 locus_tag=None,
                 organism=None,
                 project="",
                 publication="Unpublished",
                 translation_table=11
               ):
        CustomVisitor.__init__(self)
        self.contigs = {}
        self.locus_tag = locus_tag
        self.translation_table = translation_table

    def visit_feature_node(self, feature_node):
      sequence_id = feature_node.get_seqid()
      contig = self.contigs.get(sequence_id)
      if contig: # contig already exists, just try and update it
        contig.add_feature(sequence_id = sequence_id, feature_type = feature_node.get_type(), start = feature_node.get_start(),
                           end = feature_node.get_end(), strand = feature_node.get_strand(),
                           feature_attributes = feature_node.attribs,
                           locus_tag = self.locus_tag, translation_table = self.translation_table)
      else:
        contig = EMBLContig()
        successfully_added_feature = contig.add_feature(sequence_id = sequence_id, feature_type = feature_node.get_type(), start = feature_node.get_start(),
                           end = feature_node.get_end(), strand = feature_node.get_strand(),
                           feature_attributes = feature_node.attribs,
                           locus_tag = self.locus_tag, translation_table = self.translation_table)
        if successfully_added_feature:
          self.contigs[sequence_id] = contig
        else:
          pass # discard the contig because we didn't add a feature so it is empty

    def visit_region_node(self, region_node):
        pass  # for now

    def visit_comment_node(self, comment_node):
        pass  # for now

    def visit_sequence_node(self, sequence_node):
      sequence_id = sequence_node.get_description()
      contig = self.contigs.setdefault(sequence_id, EMBLContig())
      contig.add_sequence(sequence_node.get_sequence())

import os
import string
import re
import textwrap

class Convert(object):
    features_to_ignore = {'ncRNA': 1}
    feature_attributes_to_ignore = {'ID': 1, 'protein_id': 1}
    feature_attributes_translations = {'eC_number': 'EC_number'}
    feature_attributes_to_split_on_multiple_lines = {'inference': 1, 'EC_number': 1}
    feature_attributes_regex_to_change_feature_name = {}
    
    feature_attributes_inference_to_dbxref = {'similar to AA sequence:UniProtKB': 'UniProtKB/Swiss-Prot', 'protein motif:Pfam': 'PFAM', 'protein motif:CLUSTERS': "CDD", 'protein motif:Cdd': "CDD", 'protein motif:TIGRFAMs': "TIGRFAM"}
    
   
    def __init__(self, locus_tag = None, translation_table = 11):
        self.locus_tag = locus_tag
        self.translation_table = translation_table

    def blank_header(self):
      header = """\
ID   XXX; XXX; %s; genomic DNA; STD; %s; %d BP.
XX
AC   XXX;
XX
AC * _%s
XX
PR   Project:%s;
XX
DE   XXX;
XX
RN   [1]
RA   %s;
RT   "%s";
RL   %s.
XX
FH   Key             Location/Qualifiers
FH
"""
      return header
      
    def populated_header(self,
        num_bp=1,
        project="", 
        description="",
        contig_number=1, 
        authors="Pathogen Genomics", 
        title="Draft assembly annotated with Prokka",
        publication="Unpublished",
        genome_type="circular",
        classification="UNC",
        sequence_identifier=""
        ):

        header = self.blank_header()
        sequence_identifier_filtered  = re.sub(r'\W+', '', sequence_identifier)
        header_with_values = header % (genome_type, classification, num_bp,sequence_identifier_filtered, project,authors,title,publication)
        return header_with_values
        
    def source_template(self, sequence_length = None, organism = None, taxon_id = None, sequence_name = None):
        source_template = """\
FT   source          1..%d
FT                   /organism="%s"
FT                   /mol_type="genomic DNA"
FT                   /db_xref="taxon:%d"
FT                   /note="%s"
"""   % (sequence_length, organism,taxon_id,sequence_name)
        return source_template
        
    def construct_sequence(self,sequence):
      sequence_string = ''
      sequence_string += self.sequence_header(sequence)
      sequence_string += self.sequence_body(sequence)
      return sequence_string
    
    def sequence_header(self, sequence):
      sequence = sequence.lower()
      a = sequence.count('a')
      c = sequence.count('c')
      g = sequence.count('g')
      t = sequence.count('t')
      o = len(sequence) - a - c - g - t;
      return "SQ   Sequence %d BP; %d A; %d C; %d G; %d T; %d other;\n" % \
        (len(sequence), a, c, g, t, o)
      
    def sequence_body(self, sequence):
      sequence = sequence.lower()
      output = "     "
      i = 1
      for j in range(len(sequence)):
          output +=sequence[j]
          if (i) % 10 == 0:
              output += " "
          if (i) % 60 == 0 and i < len(sequence) :
              output += "%9s\n     " % (i)
          elif (i) % 60 == 0  and i == len(sequence):
             output += "%9s\n" % (i)
             return output
          i += 1

      if((i)%60 ==0):
        output += ' '*(66 -(((i-1)%60)/10) -((i-1)%60))  + "%9d\n" % (i - 1)
        return output
      else:
        output +=' '*(80-i%60-(i%60)/10-13) + "%9d\n" % (i - 1)
        return output

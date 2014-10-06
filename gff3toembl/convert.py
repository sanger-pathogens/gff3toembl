import os
import string

class Convert:
    features_to_ignore = {'ncRNA': 1}
    feature_attributes_to_ignore = {'ID': 1, 'protein_id': 1}
    feature_attributes_translations = {'eC_number': 'EC_number'}
    feature_attributes_to_split_on_multiple_lines = {'inference': 1, 'EC_number': 1}

    def __init__(self, locus_tag = None):
        self.locus_tag = locus_tag

    def blank_header(self):
      header = """\
ID   XXX; XXX; %s; genomic DNA; STD; %s; %d BP.
XX
AC   %s
XX
PR   Project:%s
XX
DE   %s contig %d
XX
RN   [1]
RA   %s
RT   "%s"
RL   %s
XX
RN   [2]
RA   %s
RT   "%s"
RL   %s
XX
CC   Data release policy http://www.sanger.ac.uk/legal/#t_2
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
        title="Draft assembly with annotation from Prokka",
        publication="Unpublished",
        genome_type="circular",
        classification="UNC",
        submitter_name="Pathogen Informatics",
        submitter_title="Direct submission",
        submitter_location="Sanger"):

        header = self.blank_header()
        header_with_values = header % (genome_type, classification, num_bp,project+str(num_bp)+str(contig_number), project, description, contig_number,authors,title,publication,submitter_name,submitter_title,submitter_location )
        return header_with_values
        
    def source_template(self, sequence_length = None, organism = None, taxon_id = None):
        source_template = """\
FT   source          1..%d
FT                   /organism="%s"
FT                   /mol_type="genomic DNA"
FT                   /db_xref="taxon:%d"
"""   % (sequence_length, organism,taxon_id)
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
        
    def feature_header(self, feature_type = None, start = None, end = None, strand = None):
      string = ""
      cmp1 = ''
      cmp2 = ''
      if strand == '-':
          cmp1 = 'complement('
          cmp2 = ')'
      string += "FT   %s%s%s%d..%d%s\n" % (feature_type, ' ' * (16-len(feature_type)), cmp1, start, end, cmp2)
      return string
      
    def construct_feature(self, feature_type = None, start = None, end = None, strand = None, feature_attributes = {}):
      feature = ''
      if feature_type in self.features_to_ignore:
        return feature
        
      feature += self.feature_header( feature_type ,start, end, strand )
      for attribute_key in feature_attributes.keys():
        feature += self.construct_feature_attribute( attribute_key = attribute_key, attribute_value = feature_attributes[attribute_key])
        
      return feature
      
    def update_locus_tag(self,attribute_value):
      if self.locus_tag == None:
        return attribute_value
      locus_tag_parts = attribute_value.split('_')
      new_attribute = self.locus_tag + '_' +str(locus_tag_parts[-1])
      return new_attribute
    
    def construct_feature_attribute(self,attribute_key = None, attribute_value = None):
      feature_string = ''
      if attribute_key in self.feature_attributes_to_ignore:      
        return feature_string
      if attribute_key in self.feature_attributes_translations:
        attribute_key = self.feature_attributes_translations[attribute_key]
      
      if attribute_key == 'locus_tag':
        attribute_value = self.update_locus_tag(attribute_value)
        
      split_attribute_values = attribute_value.split( ',')
      if attribute_key not in self.feature_attributes_to_split_on_multiple_lines:
        feature_string += self.create_multi_line_feature_attribute_string(attribute_key, split_attribute_values[0])
      else:
        for split_attribute_value in split_attribute_values:
          feature_string += self.create_multi_line_feature_attribute_string(attribute_key, split_attribute_value)
      return feature_string
      
      
    def create_multi_line_feature_attribute_string(self,attribute_key = None, attribute_value = None):
      feature_string = ''
      attribute_value = '"' + attribute_value + '"'
      
      # First line < first_line_size
      first_line_size = 55 - ( len(attribute_key))
      feature_string += "FT%s/%s=%s\n" % (' ' * 19, attribute_key, attribute_value[:first_line_size])
      if attribute_value[first_line_size:] == None:
        return feature_string
      attribute_value = attribute_value[first_line_size:]
      
      while(len(attribute_value) > 0):
        feature_string += "FT%s%s\n" % (' ' * 19, attribute_value[:57])
        if attribute_value[57:] == None:
          return feature_string
        attribute_value = attribute_value[57:]

      return feature_string    
      
      
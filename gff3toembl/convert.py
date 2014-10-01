import os

class Convert:
    def __init__(self, input_gff_file=None, output_embl_file=None):
      self.input_gff_file = input_gff_file
      self.output_embl_file = output_embl_file

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
FH\
"""
      return header
      
    def populated_header(self,
        num_bp=1, 
        accession="", 
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
        header_with_values = header % (genome_type,classification, num_bp,accession, project, description, contig_number,authors,title,publication,submitter_name,submitter_title,submitter_location )
        return header_with_values
      
    def create_output_file(self):
      pass
    
    def sequence_header(self, sequence):
      sequence = sequence.lower()
      a = sequence.count('a')
      c = sequence.count('c')
      g = sequence.count('g')
      t = sequence.count('t')
      o = len(sequence) - a - c - g - t;
      return "SQ   Sequence %d BP; %d A; %d C; %d G; %d T; %d other;" % \
        (len(sequence), a, c, g, t, o)
      
        
        
        
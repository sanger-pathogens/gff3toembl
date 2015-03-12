import re

class EMBLHeader(object):
  def __init__(self,
               authors="Pathogen Genomics",
               classification="UNC",
               genome_type="circular",
               organism=None,
               project="",
               publication="Unpublished",
               sequence_identifier="",
               sequence_length="",
               sequence_name=None,
               taxon_id=None,
               title="Draft assembly annotated with Prokka",
              ):
    self.authors=authors
    self.classification=classification
    self.genome_type=genome_type
    self.organism=organism
    self.project=project
    self.publication=publication
    self.sequence_identifier=self.remove_non_word_characters(sequence_identifier)
    self.sequence_length=sequence_length
    self.sequence_name=sequence_name
    self.taxon_id=taxon_id
    self.title=title
    self.header_template = """\
ID   XXX; XXX; {genome_type}; genomic DNA; STD; {classification}; {sequence_length} BP.
XX
AC   XXX;
XX
AC * _{sequence_identifier}
XX
PR   Project:{project};
XX
DE   XXX;
XX
RN   [1]
RA   {authors};
RT   "{title}";
RL   {publication}.
XX
FH   Key             Location/Qualifiers
FH
FT   source          1..{sequence_length}
FT                   /organism="{organism}"
FT                   /mol_type="genomic DNA"
FT                   /db_xref="taxon:{taxon_id}"
FT                   /note="{sequence_name}"
"""

  def remove_non_word_characters(self, sequence_identifier):
    return re.sub(r'\W+', '', sequence_identifier)

  def format(self):
    return self.header_template.format(**self.__dict__)

class EMBLHeader(object):
  def __init__(self):
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

  def format(self):
    return self.header_template.format(**self.__dict__)

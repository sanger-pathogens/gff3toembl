import re
from textwrap import TextWrapper

class EMBLContig(object):
  def __init__(self):
    self.header = None
    self.features = {}
    self.sequence = None

  def format(self):
    try:
      header = self.header.format()
    except AttributeError:
      raise ValueError("Could not format contig, no header data found")
    feature_strings = [feature.format() for feature in self.sorted_features()]
    features = "".join(feature_strings)
    try:
      sequence = self.sequence.format()
    except AttributeError:
      raise ValueError("Could not format contig, no sequence data found")
    return header + features + sequence

  def add_header(self, **kwargs):
    if self.header != None:
      raise ValueError("Contig already has header data")
    header = EMBLHeader(**kwargs)
    self.header = header

  def add_feature(self, sequence_id, **kwargs):
    feature = EMBLFeature(**kwargs)
    unique_feature_reference = "{}_{}_{}_{}".format(sequence_id, feature.feature_type, feature.start, feature.end)
    if unique_feature_reference in self.features:
      # we're already seen a feature in this region so don't add another
      return False
    elif feature.format() == None:
      # some feature types should be ignored; format() returns None in these cases
      return False
    else:
      self.features[unique_feature_reference] = feature
      return True

  def add_sequence(self, sequence_string):
    if self.sequence != None:
      raise ValueError("Contig already has sequence data")
    sequence = EMBLSequence(sequence_string)
    self.sequence = sequence

  def sorted_features(self):
    # Features should be sorted by start and then by end irrespective of strand
    def compare_features(feature_1, feature_2):
      if feature_1.start < feature_2.start:
        return -1
      elif feature_1.start > feature_2.start:
        return 1
      else:
        return feature_2.end - feature_1.end
    return sorted(self.features.values(), cmp=compare_features)

class EMBLFeature(object):
  inference_to_db_xref_map = {
          'similar to AA sequence:UniProtKB': 'UniProtKB/Swiss-Prot',
          'protein motif:Pfam': 'PFAM',
          'protein motif:CLUSTERS': "CDD",
          'protein motif:Cdd': "CDD",
          'protein motif:TIGRFAMs': "TIGRFAM"
  }

  def __init__(self, feature_type, start, end, strand, feature_attributes,
               locus_tag=None, translation_table=11):
    # Picks a feature builder and builds the feature
    # Most features are built with a default but some are either a little different or
    # should just be ignored
    feature_builder = self.pick_feature_builder(feature_type)
    feature_builder(feature_type=feature_type, start=start, end=end, strand=strand,
                    feature_attributes=feature_attributes, locus_tag=locus_tag,
                    translation_table=translation_table)

  def pick_feature_builder(self, feature_type):
    feature_builders = {
      'CDS': self.create_CDS_feature,
      'source': self.create_source_feature,
      'ncRNA': self.create_empty_feature
    }
    return feature_builders.get(feature_type, self.create_default_feature)

  def create_default_feature(self, feature_type, start, end, strand, feature_attributes, locus_tag, translation_table):
    self.feature_type = feature_type
    self.start = start
    self.end = end
    self.strand = strand
    self.locus_tag = locus_tag
    self.translation_table  = translation_table
    self.attributes = []
    for attribute_key, attribute_value in feature_attributes.items():
      attribute_creator = self.lookup_attribute_creator(attribute_key)
      new_attributes = attribute_creator(attribute_key, attribute_value)
      self.attributes += new_attributes

  def create_CDS_feature(self, **kwargs):
    self.create_default_feature(**kwargs)
    self.attributes += self.create_translation_table_attributes('transl_table', self.translation_table)

  def create_source_feature(self, feature_type, start, end, strand, feature_attributes, locus_tag, translation_table):
    self.feature_type = feature_type
    self.start = start
    self.end = end
    self.strand = strand
    self.locus_tag = locus_tag
    self.translation_table  = translation_table

    organism = feature_attributes['organism']
    db_xref = feature_attributes['db_xref']
    note = feature_attributes['note']

    # We hard code the order and composition of attributes for source features
    # Source features are only created as part of the header
    self.attributes = [("organism", organism), ("mol_type", "genomic DNA"), ("db_xref", db_xref), ("note", note)]

  def create_empty_feature(self, **kwargs):
    # Some features should be ignored.  This is how this is done
    self.format = lambda: None

  def format(self):
    coordinates = self.format_coordinates(self.start, self.end, self.strand)
    header_string = "FT   {feature_type: <16}{coordinates}".format( feature_type=self.feature_type,
                                                                     coordinates=coordinates)
    attribute_strings = [header_string]
    for attribute_key,attribute_value in self.attributes:
      attribute_strings.append(self.format_attribute(attribute_key, attribute_value))

    return '\n'.join(attribute_strings) + '\n'

  def format_attribute(self, key, value):
    # Looks up a formatter for an attribute and formats the attribute
    # Some attributes are formatted a little differently
    formatter = self.lookup_attribute_formatter(key)
    return formatter(key, value)

  def lookup_attribute_formatter(self, attribute_type):
    formatters = {
      'transl_table': self.translation_table_attribute_formatter
    }
    return formatters.get(attribute_type, self.default_attribute_formatter)

  def translation_table_attribute_formatter(self, key, value):
    # transl_table attributes do not have their values in quotes
    wrapper = TextWrapper()
    wrapper.initial_indent='FT                   '
    wrapper.subsequent_indent='FT                   '
    wrapper.width=79
    attribute_text_template='/{attribute_key}={attribute_value}'
    attribute_text=attribute_text_template.format(attribute_key=key, attribute_value=value)
    return wrapper.fill(attribute_text)

  def default_attribute_formatter(self, key, value):
    wrapper = TextWrapper()
    wrapper.initial_indent='FT                   '
    wrapper.subsequent_indent='FT                   '
    wrapper.width=79
    attribute_text_template='/{attribute_key}="{attribute_value}"'
    attribute_text=attribute_text_template.format(attribute_key=key, attribute_value=value)
    return wrapper.fill(attribute_text)

  def format_coordinates(self, start, end, strand):
    if strand == '-':
      return "complement({start}..{end})".format(start=start, end=end)
    else:
      return "{start}..{end}".format(start=start, end=end)

  def lookup_attribute_creator(self, attribute_key):
    # These functions take attributes and reformat them into a list
    # of (key, values) which are later formatted into strings by other
    # methods.  There is quite a lot of variation between these such as
    # whether to keep more than one value for a given attribute type.
    attribute_creator_table = {
      'product': self.create_product_attributes,
      'locus_tag': self.create_locus_tag_attributes,
      'eC_number': self.create_EC_number_attributes,
      'inference': self.create_inference_attributes,
      'protein_id': self.ignore_attributes,
      'ID': self.ignore_attributes
    }
    return attribute_creator_table.get(attribute_key, self.create_default_attributes)

  def create_default_attributes(self, attribute_key, attribute_value):
    all_attribute_values = attribute_value.split(',')
    first_attribute_value = all_attribute_values[0]
    return [(attribute_key, first_attribute_value)]

  def create_product_attributes(self, attribute_key, attribute_value):
    def remove_hypotheticals(value):
      return value != 'hypothetical protein'
    def replace_unknown_with_uncharacterised(value):
      return value.replace("nknown","ncharacterised")
    # attribute_value may be a comma deliminated list of values
    # only some of which might be valid
    attribute_values = attribute_value.split(',')
    attribute_values = filter(remove_hypotheticals, attribute_values)
    attribute_values = map(replace_unknown_with_uncharacterised, attribute_values)
    chosen_value = attribute_values[0] if len(attribute_values) > 0 else 'Uncharacterised protein'
    return [('product', chosen_value)]

  def create_locus_tag_attributes(self, attribute_key, attribute_value):
    if self.locus_tag == None:
      return [('locus_tag', attribute_value)]
    else:
      attribute_value_suffix = attribute_value.split('_')[-1]
      return [('locus_tag', "{}_{}".format(self.locus_tag, attribute_value_suffix))]

  def create_EC_number_attributes(self, attribute_key, attribute_value):
    attribute_values = attribute_value.split(',')
    def deduplicate_values(values):
      return list(set(values))
    attribute_values = deduplicate_values(attribute_values)
    return [('EC_number', value) for value in attribute_values]

  def create_inference_attributes(self, attribute_key, attribute_value):
    attribute_values = attribute_value.split(',')
    attributes = []
    for value in attribute_values:
      if self.should_convert_to_db_xref(value):
        attributes.append(('db_xref', self.convert_to_db_xref(value)))
      else:
        attributes.append(('inference', value))
    return attributes

  def ignore_attributes(self, attribute_key, attribute_value):
    return []

  def should_convert_to_db_xref(self, attribute_value):
    for search_text in self.inference_to_db_xref_map:
      if search_text in attribute_value:
        return True
    return False

  def convert_to_db_xref(self, attribute_value):
    for search_text, replacement_text in self.inference_to_db_xref_map.items():
      if search_text in attribute_value:
        return attribute_value.replace(search_text, replacement_text)
    raise ValueError("Failed to convert inference attribute '%s' to db_xref" % attribute_value)

  def create_translation_table_attributes(self, attribute_key, attribute_value):
    return [('transl_table', attribute_value)]

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

    source_attributes = self.build_source_attributes(organism, taxon_id, sequence_name)
    self.source_feature = EMBLFeature(feature_type='source', start=1, end=sequence_length,
                                      strand='+', feature_attributes=source_attributes)

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
"""

  def remove_non_word_characters(self, sequence_identifier):
    return re.sub(r'\W+', '', sequence_identifier)

  def format(self):
    return self.header_template.format(**self.__dict__) + self.source_feature.format()

  def build_source_attributes(self, organism, taxon_id, sequence_name):
    def empty_string_if_none(value):
      return value if value else ''
    organism = empty_string_if_none(organism)
    taxon_id = empty_string_if_none(taxon_id)
    sequence_name = empty_string_if_none(sequence_name)
    return {"organism": organism, "db_xref": "taxon:{}".format(taxon_id), "note": sequence_name}

class EMBLSequence(object):

  def __init__(self, sequence_string):
    nucleotide_counts = self.calculate_nucleotide_counts(sequence_string)
    self.header = self.format_header(nucleotide_counts)
    self.body = self.format_sequence_body(sequence_string)
    self.length = len(sequence_string)

  def format(self):
    return self.header + '\n' + self.body

  def calculate_nucleotide_counts(self, sequence):
    sequence = sequence.lower()
    counts = {}
    counts['a'] = sequence.count('a')
    counts['c'] = sequence.count('c')
    counts['g'] = sequence.count('g')
    counts['t'] = sequence.count('t')
    count_of_acgt = sum(counts.values())
    counts['other'] = len(sequence) - count_of_acgt
    return counts

  def format_header(self, nucleotide_counts):
    template = "SQ   Sequence {total} BP; {a} A; {c} C; {g} G; {t} T; {other} other;"
    total_counts = sum(nucleotide_counts.values())
    nucleotide_counts['total'] = total_counts
    return template.format(**nucleotide_counts)

  def format_sequence_body(self, sequence_string):
    sequence_string = sequence_string.lower()
    lines = self.split_sequence(sequence_string)
    def format_a_line(line):
      # a line looks like:
      # (["1234567890", "12345", '', '', '', ''], 15)
      # and should look like
      # "     1234567890 12345                                                         15"
      blocks_of_sequence, end_of_line = line
      format_arguments = blocks_of_sequence + [end_of_line]
      return "     {:<10} {:<10} {:<10} {:<10} {:<10} {:<10} {:>9}".format(*format_arguments)
    formatted_lines = map(format_a_line, lines)
    return '\n'.join(formatted_lines) + '\n'

  def split_line_of_sequence(self, line_of_sequence):
    # Turns "123456789012345" into ["1234567890", "12345", '', '', '', '']
    splits = []
    line_breaks = range(0, 60, 10)
    for line_break in line_breaks:
      split = line_of_sequence[line_break:line_break+10]
      splits.append(split)
    return splits

  def split_sequence(self, sequence_string):
    splits = []
    sequence_length = len(sequence_string)
    for start_of_line in range(0, sequence_length, 60):
      # might not actually be the end of the line if the line isn't long enough
      end_of_line = start_of_line + 60
      line_of_sequence = sequence_string[start_of_line:end_of_line]
      length_of_line = len(line_of_sequence)
      end_of_line = start_of_line + length_of_line # actually end of the line
      splits.append((self.split_line_of_sequence(line_of_sequence), end_of_line))
    return splits

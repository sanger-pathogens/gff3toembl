import unittest
import pdb
from gff3toembl.EMBLContig import EMBLFeature, EMBLHeader

class TestEMBLHeader(unittest.TestCase):

  def test_initialize_header_object(self):
    header = EMBLHeader(
      authors="John Doe",
      classification="UNC",
      genome_type="circular",
      organism="My organism",
      project="PRJ1234",
      publication="Unpublished",
      sequence_identifier="**contig123",
      sequence_length=8,
      sequence_name="chromX",
      taxon_id=5678,
      title="My title"
    )

    self.assertEqual(header.authors, "John Doe")
    self.assertEqual(header.classification, "UNC")
    self.assertEqual(header.genome_type, "circular")
    self.assertEqual(header.organism, "My organism")
    self.assertEqual(header.project, "PRJ1234")
    self.assertEqual(header.publication, "Unpublished")
    self.assertEqual(header.sequence_identifier, "contig123") # Removed non-word characters
    self.assertEqual(header.sequence_length, 8)
    self.assertEqual(header.sequence_name, "chromX")
    self.assertEqual(header.taxon_id, 5678)
    self.assertEqual(header.title, "My title")

  def test_remove_non_word_characters(self):
    header = EMBLHeader()
    test_cases = [
      ('foo', 'foo'),
      ('#foo', 'foo'),
      ('#fo!o', 'foo'),
      ('fo#o', 'foo'),
      ('foo##', 'foo'),
      ("*!#foo", 'foo')
    ]
    for test_string, expected_result in test_cases:
      self.assertEqual(header.remove_non_word_characters(test_string), expected_result)

  def test_format(self):
    header = EMBLHeader()

    header.authors="John Doe"
    header.classification="UNC"
    header.genome_type="circular"
    header.organism="My organism"
    header.project="PRJ1234"
    header.publication="Unpublished"
    header.sequence_identifier="contig123"
    header.sequence_length=1234
    header.sequence_name="chromX"
    header.taxon_id=5678
    header.title="My title"

    expected_header = """\
ID   XXX; XXX; circular; genomic DNA; STD; UNC; 1234 BP.
XX
AC   XXX;
XX
AC * _contig123
XX
PR   Project:PRJ1234;
XX
DE   XXX;
XX
RN   [1]
RA   John Doe;
RT   "My title";
RL   Unpublished.
XX
FH   Key             Location/Qualifiers
FH
FT   source          1..1234
FT                   /organism="My organism"
FT                   /mol_type="genomic DNA"
FT                   /db_xref="taxon:5678"
FT                   /note="chromX"
"""
    calculated_header = header.format()
    for calculated_line,expected_line in zip(calculated_header.split('\n'), expected_header.split('\n')):
      self.assertEqual(calculated_line, expected_line)
    self.assertEqual(len(calculated_header), len(expected_header))

class TestEMBLFeature(unittest.TestCase):

  def test_format(self):
    feature = EMBLFeature()
    feature.feature_type = "feature_type"
    feature.start = 1
    feature.end = 10
    feature.strand = ''
    feature.attributes = [
      ("attributeA", "foo"),
      ("attributeB", 'bar'),
      ("attributeB", 'baz')
    ]
    calculated_string = feature.format()
    expected_string = """\
FT   feature_type    1..10
FT                   /attributeA="foo"
FT                   /attributeB="bar"
FT                   /attributeB="baz"
"""

    for calculated_line,expected_line in zip(calculated_string.split('\n'), expected_string.split('\n')):
      self.assertEqual(calculated_line, expected_line)
    self.assertEqual(len(calculated_string), len(expected_string))

  def test_should_ignore_feature_type(self):
    feature = EMBLFeature()
    self.assertTrue(feature.should_ignore_feature('ID'))
    self.assertTrue(feature.should_ignore_feature('protein_id'))
    self.assertFalse(feature.should_ignore_feature('other'))

  def test_reformat_feature_type(self):
    feature = EMBLFeature()
    self.assertEqual(feature.reformat_feature_type('eC_number'), 'EC_number')
    self.assertEqual(feature.reformat_feature_type('anything_else'), 'anything_else')

  def test_format_attribute(self):
    feature = EMBLFeature()
    calculated_string = feature.format_attribute('attributeA', 'foo')
    expected_string = 'FT                   /attributeA="foo"'
    self.assertEqual(calculated_string, expected_string)

  def test_lookup_attribute_creator(self):
    feature = EMBLFeature()
    self.assertEqual(feature.lookup_attribute_creator('some_key'),
                     feature.create_default_attributes)
    self.assertEqual(feature.lookup_attribute_creator('product'),
                     feature.create_product_attributes)
    self.assertEqual(feature.lookup_attribute_creator('locus_tag'),
                     feature.create_locus_tag_attributes)
    self.assertEqual(feature.lookup_attribute_creator('eC_number'),
                     feature.create_EC_number_attributes)

  def test_create_product_attributes(self):
    feature = EMBLFeature()
    test_cases = [
      ('abc,efg,hij', [('product', "abc")]),
      ('hypothetical protein,efg,hij', [('product', "efg")]),
      ('efg,hypothetical protein,hij', [('product', "efg")]),
      ('hypothetical protein,hypothetical protein,hij', [('product', "hij")]),
      ('hypothetical protein', [('product', "Uncharacterised protein")]),
      ('hypothetical protein,Unknown protein abc', [('product', "Uncharacterised protein abc")]),
      ('hypothetical protein,unknown protein abc', [('product', "uncharacterised protein abc")])
    ]
    for test_case, expected_result in test_cases:
      calculated_result = feature.create_product_attributes('product', test_case)
      self.assertEqual(calculated_result, expected_result)

  def test_create_locus_tag_attributes(self):
    feature = EMBLFeature()
    calculated_attributes = feature.create_locus_tag_attributes('locus_tag', 'ABC123')
    expected_attributes = [('locus_tag', 'ABC123')]
    self.assertEqual(calculated_attributes, expected_attributes)

    calculated_attributes = feature.create_locus_tag_attributes('locus_tag', 'ABC123_XYZ')
    expected_attributes = [('locus_tag', 'ABC123_XYZ')]
    self.assertEqual(calculated_attributes, expected_attributes)

    calculated_attributes = feature.create_locus_tag_attributes('locus_tag', 'ABC_123_XYZ')
    expected_attributes = [('locus_tag', 'ABC_123_XYZ')]
    self.assertEqual(calculated_attributes, expected_attributes)

    feature.locus_tag = "some_other_tag"

    calculated_attributes = feature.create_locus_tag_attributes('locus_tag', 'ABC123')
    expected_attributes = [('locus_tag', 'some_other_tag_ABC123')]
    self.assertEqual(calculated_attributes, expected_attributes)

    calculated_attributes = feature.create_locus_tag_attributes('locus_tag', 'ABC123_XYZ')
    expected_attributes = [('locus_tag', 'some_other_tag_XYZ')]
    self.assertEqual(calculated_attributes, expected_attributes)

    calculated_attributes = feature.create_locus_tag_attributes('locus_tag', 'ABC_123_XYZ')
    expected_attributes = [('locus_tag', 'some_other_tag_XYZ')]
    self.assertEqual(calculated_attributes, expected_attributes)

    feature.locus_tag = None

    calculated_attributes = feature.create_locus_tag_attributes('locus_tag', 'ABC123')
    expected_attributes = [('locus_tag', 'ABC123')]
    self.assertEqual(calculated_attributes, expected_attributes)

    calculated_attributes = feature.create_locus_tag_attributes('locus_tag', 'ABC123_XYZ')
    expected_attributes = [('locus_tag', 'ABC123_XYZ')]
    self.assertEqual(calculated_attributes, expected_attributes)

    calculated_attributes = feature.create_locus_tag_attributes('locus_tag', 'ABC_123_XYZ')
    expected_attributes = [('locus_tag', 'ABC_123_XYZ')]
    self.assertEqual(calculated_attributes, expected_attributes)

  def test_create_EC_number_attributes(self):
    feature = EMBLFeature()
    calculated_attributes = feature.create_EC_number_attributes('eC_number', '123')
    expected_attributes = [('EC_number', '123')]
    self.assertEqual(calculated_attributes, expected_attributes)

    calculated_attributes = feature.create_EC_number_attributes('eC_number', '123,ABC')
    expected_attributes = [('EC_number', '123'), ('EC_number', 'ABC')]
    self.assertEqual(calculated_attributes, expected_attributes)

    calculated_attributes = feature.create_EC_number_attributes('eC_number', '123,123')
    expected_attributes = [('EC_number', '123')]
    self.assertEqual(calculated_attributes, expected_attributes)

    calculated_attributes = feature.create_EC_number_attributes('eC_number', '123,ABC,123')
    expected_attributes = [('EC_number', '123'), ('EC_number', 'ABC')]
    self.assertEqual(calculated_attributes, expected_attributes)

    calculated_attributes = feature.create_EC_number_attributes('eC_number', 'B,A,A')
    expected_attributes = [('EC_number', 'B'), ('EC_number', 'A')]
    self.assertEqual(calculated_attributes, expected_attributes)

    calculated_attributes = feature.create_EC_number_attributes('eC_number', 'B,A,B')
    expected_attributes = [('EC_number', 'B'), ('EC_number', 'A')]
    self.assertEqual(calculated_attributes, expected_attributes)

  def test_format_coordinates(self):
    feature = EMBLFeature()
    calculated_coordinates = feature.format_coordinates(1, 10, '')
    expected_coordinates = '1..10'
    self.assertEqual(calculated_coordinates, expected_coordinates)

    calculated_coordinates = feature.format_coordinates(1, 10, '-')
    expected_coordinates = 'compliment(1..10)'
    self.assertEqual(calculated_coordinates, expected_coordinates)

    calculated_coordinates = feature.format_coordinates(1, 10, '***NONSENCE***')
    expected_coordinates = '1..10'
    self.assertEqual(calculated_coordinates, expected_coordinates)

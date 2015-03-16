import unittest
import mock
from gff3toembl.EMBLConverter import EMBLConverter
from gff3toembl.EMBLContig import EMBLContig

class TestEMBLConverter(unittest.TestCase):

  def mock_feature_node(self, seqid, feature_type, start, end, strand, attributes):
    mock_methods = {
      'get_seqid.return_value': seqid,
      'get_type.return_value': feature_type,
      'get_start.return_value': start,
      'get_end.return_value': end,
      'get_strand.return_value': strand
    }
    return mock.Mock(attribs={'attribute_key': attributes}, **mock_methods)

  @mock.patch('gff3toembl.EMBLContig.EMBLFeature')
  def test_visit_features_node_feature_ignored(self, embl_feature_mock):
    converter = EMBLConverter(None)

    feature_node = self.mock_feature_node(1, 'ignored_type', 1, 100, '', {'attr_k1': 'attr_v1'})
    embl_feature_mock.return_value.format.return_value = None
    converter.visit_feature_node(feature_node)
    self.assertEqual(converter.contigs, {})

  @mock.patch('gff3toembl.EMBLContig.EMBLFeature')
  def test_visit_features_node_feature_included(self, embl_feature_mock):
    converter = EMBLConverter(None)

    feature_node = self.mock_feature_node(1, 'type_1', 1, 100, '', {'attr_k1': 'attr_v1'})
    embl_feature_mock.return_value.format.return_value = 'Feature_string'
    converter.visit_feature_node(feature_node)
    self.assertEqual(converter.contigs.keys(), [1])
    self.assertIsInstance(converter.contigs[1], EMBLContig)

  @mock.patch('gff3toembl.EMBLContig.EMBLFeature')
  def test_visit_features_node_one_feature_included_one_ignored(self, embl_feature_mock):
    converter = EMBLConverter(None)

    feature_node = self.mock_feature_node(1, 'type_1', 1, 100, '', {'attr_k1': 'attr_v1'})
    embl_feature_mock.return_value.format.return_value = 'Feature_string'
    converter.visit_feature_node(feature_node)

    feature_node = self.mock_feature_node(2, 'ignored_type', 101, 200, '', {'attr_k1': 'attr_v1'})
    embl_feature_mock.return_value.format.return_value = None
    converter.visit_feature_node(feature_node)
    self.assertEqual(converter.contigs.keys(), [1])
    self.assertIsInstance(converter.contigs[1], EMBLContig)

  @mock.patch('gff3toembl.EMBLContig.EMBLFeature')
  def test_visit_features_node_two_different_features(self, embl_feature_mock):
    converter = EMBLConverter(None)

    feature_node = self.mock_feature_node(1, 'type_1', 1, 100, '', {'attr_k1': 'attr_v1'})
    embl_feature_mock.return_value.format.return_value = 'Feature_string'
    converter.visit_feature_node(feature_node)

    feature_node = self.mock_feature_node(2, 'type_2', 101, 200, '', {'attr_k1': 'attr_v1'})
    embl_feature_mock.return_value.format.return_value = "Another_feature_string"
    converter.visit_feature_node(feature_node)
    self.assertEqual(converter.contigs.keys(), [1, 2])
    self.assertIsInstance(converter.contigs[1], EMBLContig)
    self.assertIsInstance(converter.contigs[2], EMBLContig)

  @mock.patch('gff3toembl.EMBLContig.EMBLFeature')
  def test_visit_features_node_repeated_features(self, embl_feature_mock):
    converter = EMBLConverter(None)

    feature_node = self.mock_feature_node(1, 'type_1', 1, 100, '', {'attr_k1': 'attr_v1'})
    embl_feature_mock.return_value.format.return_value = 'Feature_string'
    converter.visit_feature_node(feature_node)

    feature_node = self.mock_feature_node(1, 'type_1', 1, 100, '', {'attr_k1': 'attr_v1'})
    embl_feature_mock.return_value.format.return_value = 'Feature_string'
    converter.visit_feature_node(feature_node)
    self.assertEqual(converter.contigs.keys(), [1])
    self.assertIsInstance(converter.contigs[1], EMBLContig)

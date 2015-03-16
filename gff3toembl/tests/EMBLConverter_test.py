import unittest
import mock
from gff3toembl.EMBLConverter import EMBLConverter

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

  def test_update_features_seen(self):
    converter = EMBLConverter(None)
    assert converter.update_features_seen("SEQ_ID", "FEATURE_TYPE", 1, 100) == False
    assert converter.update_features_seen("SEQ_ID", "FEATURE_TYPE", 1, 100) == True
    assert converter.update_features_seen("ANOTHER_SEQ_ID", "FEATURE_TYPE", 1, 100) == False
    assert converter.update_features_seen("SEQ_ID", "FEATURE_TYPE", 10, 110) == False
    assert converter.update_features_seen("SEQ_ID", "ANOTHER_FEATURE_TYPE", 1, 100) == False
    assert converter.update_features_seen("ANOTHER_SEQ_ID", "ANOTHER_FEATURE_TYPE", 1001, 1101) == False
    assert converter.update_features_seen("SEQ_ID", "FEATURE_TYPE", 1, 100) == True

  @mock.patch('gff3toembl.EMBLConverter.EMBLFeature')
  def test_visit_features_node_feature_ignored(self, embl_feature_mock):
    converter = EMBLConverter(None)

    self.assertItemsEqual(converter.features_seen, [])
    self.assertItemsEqual(converter.feats, {})
    feature_node = self.mock_feature_node(1, 'ignored_type', 1, 100, '', {'attr_k1': 'attr_v1'})
    embl_feature_mock.return_value.format.return_value = None
    converter.visit_feature_node(feature_node)
    self.assertItemsEqual(converter.features_seen, [])
    self.assertItemsEqual(converter.feats, {})

  @mock.patch('gff3toembl.EMBLConverter.EMBLFeature')
  def test_visit_features_node_feature_included(self, embl_feature_mock):
    converter = EMBLConverter(None)

    self.assertItemsEqual(converter.features_seen, [])
    self.assertItemsEqual(converter.feats, {})
    feature_node = self.mock_feature_node(1, 'type_1', 1, 100, '', {'attr_k1': 'attr_v1'})
    embl_feature_mock.return_value.format.return_value = 'Feature_string'
    converter.visit_feature_node(feature_node)
    self.assertItemsEqual(converter.features_seen, ['1_type_1_1_100'])
    self.assertItemsEqual(converter.feats, {1: 'Feature_string'})

  @mock.patch('gff3toembl.EMBLConverter.EMBLFeature')
  def test_visit_features_node_one_feature_included_one_ignored(self, embl_feature_mock):
    converter = EMBLConverter(None)

    self.assertItemsEqual(converter.features_seen, [])
    self.assertItemsEqual(converter.feats, {})
    feature_node = self.mock_feature_node(1, 'type_1', 1, 100, '', {'attr_k1': 'attr_v1'})
    embl_feature_mock.return_value.format.return_value = 'Feature_string'
    converter.visit_feature_node(feature_node)
    self.assertItemsEqual(converter.features_seen, ['1_type_1_1_100'])
    self.assertItemsEqual(converter.feats, {1: 'Feature_string'})

    feature_node = self.mock_feature_node(2, 'ignored_type', 101, 200, '', {'attr_k1': 'attr_v1'})
    embl_feature_mock.return_value.format.return_value = None
    converter.visit_feature_node(feature_node)
    self.assertItemsEqual(converter.features_seen, ['1_type_1_1_100'])
    self.assertItemsEqual(converter.feats, {1: 'Feature_string'})

  @mock.patch('gff3toembl.EMBLConverter.EMBLFeature')
  def test_visit_features_node_two_different_features(self, embl_feature_mock):
    converter = EMBLConverter(None)

    self.assertItemsEqual(converter.features_seen, [])
    self.assertItemsEqual(converter.feats, {})
    feature_node = self.mock_feature_node(1, 'type_1', 1, 100, '', {'attr_k1': 'attr_v1'})
    embl_feature_mock.return_value.format.return_value = 'Feature_string'
    converter.visit_feature_node(feature_node)
    self.assertItemsEqual(converter.features_seen, ['1_type_1_1_100'])
    self.assertItemsEqual(converter.feats, {1: 'Feature_string'})

    feature_node = self.mock_feature_node(2, 'type_2', 101, 200, '', {'attr_k1': 'attr_v1'})
    embl_feature_mock.return_value.format.return_value = "Another_feature_string"
    converter.visit_feature_node(feature_node)
    self.assertItemsEqual(converter.features_seen, ['1_type_1_1_100', '2_type_2_101_200'])
    self.assertItemsEqual(converter.feats, {1: 'Feature_string', 2: 'Another_feature_string'})

  @mock.patch('gff3toembl.EMBLConverter.EMBLFeature')
  def test_visit_features_node_repeated_features(self, embl_feature_mock):
    converter = EMBLConverter(None)

    self.assertItemsEqual(converter.features_seen, [])
    self.assertItemsEqual(converter.feats, {})
    feature_node = self.mock_feature_node(1, 'type_1', 1, 100, '', {'attr_k1': 'attr_v1'})
    embl_feature_mock.return_value.format.return_value = 'Feature_string'
    converter.visit_feature_node(feature_node)
    self.assertItemsEqual(converter.features_seen, ['1_type_1_1_100'])
    self.assertItemsEqual(converter.feats, {1: 'Feature_string'})

    feature_node = self.mock_feature_node(1, 'type_1', 1, 100, '', {'attr_k1': 'attr_v1'})
    embl_feature_mock.return_value.format.return_value = 'Feature_string'
    converter.visit_feature_node(feature_node)
    self.assertItemsEqual(converter.features_seen, ['1_type_1_1_100'])
    self.assertItemsEqual(converter.feats, {1: 'Feature_string'})

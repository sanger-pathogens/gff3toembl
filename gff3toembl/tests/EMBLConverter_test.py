import unittest
from gff3toembl.EMBLConverter import EMBLConverter

class TestEMBLConverter(unittest.TestCase):

  def test_update_features_seen(self):
    converter = EMBLConverter(None)
    assert converter.update_features_seen("SEQ_ID", "FEATURE_TYPE", 1, 100) == False
    assert converter.update_features_seen("SEQ_ID", "FEATURE_TYPE", 1, 100) == True
    assert converter.update_features_seen("ANOTHER_SEQ_ID", "FEATURE_TYPE", 1, 100) == False
    assert converter.update_features_seen("SEQ_ID", "FEATURE_TYPE", 10, 110) == False
    assert converter.update_features_seen("SEQ_ID", "ANOTHER_FEATURE_TYPE", 1, 100) == False 
    assert converter.update_features_seen("ANOTHER_SEQ_ID", "ANOTHER_FEATURE_TYPE", 1001, 1101) == False
    assert converter.update_features_seen("SEQ_ID", "FEATURE_TYPE", 1, 100) == True

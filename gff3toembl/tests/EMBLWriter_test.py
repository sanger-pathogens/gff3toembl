import unittest
import sys
import os
import filecmp
from gff3toembl.EMBLWriter import EMBLWriter
from gff3toembl import convert

modules_dir = os.path.dirname(os.path.abspath(convert.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestEMBLWriter(unittest.TestCase):

    def test_blank_header(self):
        '''test that the script will convert from GFF3 to EMBL'''
        
        emblwriter = EMBLWriter(os.path.join(data_dir,'single_feature.gff'), 
           'Organism', 
           1234, 
           'My project', 
           'My description', 
           'John', 
           'Some title',  
           'Some journal', 
           'circular', 
           'PROK', 
           'Jane',
           'My institute',  
           'UK' )
        print emblwriter.parse_and_run()
        assert  emblwriter.parse_and_run() == expected_output

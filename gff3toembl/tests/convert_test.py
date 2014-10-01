import unittest
import sys
import os
import filecmp
from gff3toembl import convert

modules_dir = os.path.dirname(os.path.abspath(convert.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestConvert(unittest.TestCase):
    def test_overall_large_conversion(self):
        '''test converting a large gff file into embl'''
        converter = convert.Convert(input_gff_file=os.path.join(data_dir, 'large_convert_test.gff'), output_embl_file='large_convert_test.embl')
        converter.create_output_file()
        self.assertTrue(filecmp.cmp('large_convert_test.embl', os.path.join(data_dir, 'expected_large_convert_test.embl')))
        os.unlink('large_convert_test.embl')
    
    def test_blank_header(self):
        '''test that you can get the correct template header out'''
        converter = convert.Convert(input_gff_file=os.path.join(data_dir, 'large_convert_test.gff'), output_embl_file='large_convert_test.embl')
        expected_header = """\
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
        assert converter.blank_header() == expected_header
    
    def test_populate_header(self):
        converter = convert.Convert(input_gff_file=os.path.join(data_dir, 'large_convert_test.gff'), output_embl_file='large_convert_test.embl')
        actual_populated_header = converter.populated_header(num_bp=1234, 
          accession="ABC123", 
          project="PRJ1234", 
          description="One line description",
          contig_number=1, 
          authors="John Doe", 
          title="My title",
          publication="Unpublished",
          genome_type="circular",
          classification="UNC",
          submitter_name="Jane Doe",
          submitter_title="Direct submission",
          submitter_location="Sanger"
          )
        expected_populated_header = """\
ID   XXX; XXX; circular; genomic DNA; STD; UNC; 1234 BP.
XX
AC   ABC123
XX
PR   Project:PRJ1234
XX
DE   One line description contig 1
XX
RN   [1]
RA   John Doe
RT   "My title"
RL   Unpublished
XX
RN   [2]
RA   Jane Doe
RT   "Direct submission"
RL   Sanger
XX
CC   Data release policy http://www.sanger.ac.uk/legal/#t_2
XX
FH   Key             Location/Qualifiers
FH\
"""
        assert actual_populated_header == expected_populated_header
        
        
import unittest
import sys
import os
import filecmp
from gff3toembl import convert

modules_dir = os.path.dirname(os.path.abspath(convert.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestConvert(unittest.TestCase):

    def test_blank_header(self):
        '''test that you can get the correct template header out'''
        converter = convert.Convert()
        expected_header = """\
ID   XXX; XXX; %s; genomic DNA; STD; %s; %d BP.
XX
AC   XXX;
XX
AC * _%s
XX
PR   Project:%s;
XX
DE   XXX;
XX
RN   [1]
RA   %s;
RT   "%s";
RL   %s.
XX
FH   Key             Location/Qualifiers
FH
"""
        assert converter.blank_header() == expected_header
    
    def test_populate_header(self):
        converter = convert.Convert()
        actual_populated_header = converter.populated_header(num_bp=1234, 
          project="PRJ1234", 
          description="One line description",
          contig_number=1, 
          authors="John Doe", 
          title="My title",
          publication="Unpublished",
          genome_type="circular",
          classification="UNC",
          sequence_identifier="contig123"
          )
        expected_populated_header = """\
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
"""
        assert actual_populated_header == expected_populated_header
        
    def test_construct_sequence(self):
        converter = convert.Convert()
        assert converter.construct_sequence("AAAACCCGGTNN") == """\
SQ   Sequence 12 BP; 4 A; 3 C; 2 G; 1 T; 2 other;
     aaaacccggt nn                                                            12
"""

    def test_source_template(self):
        converter = convert.Convert()
        assert converter.source_template(1234,"My organism", 5678,"chromX") == """\
FT   source          1..1234
FT                   /organism="My organism"
FT                   /mol_type="genomic DNA"
FT                   /db_xref="taxon:5678"
FT                   /note="chromX"
"""

    def test_sequence_header(self):
        converter = convert.Convert()
        assert converter.sequence_header("AAAACCCGGTNN") == "SQ   Sequence 12 BP; 4 A; 3 C; 2 G; 1 T; 2 other;\n"
        assert converter.sequence_header("AAAAaaaaAAAA") == "SQ   Sequence 12 BP; 12 A; 0 C; 0 G; 0 T; 0 other;\n"
        assert converter.sequence_header("------------") == "SQ   Sequence 12 BP; 0 A; 0 C; 0 G; 0 T; 12 other;\n"
        assert converter.sequence_header("acgtACGTtttT") == "SQ   Sequence 12 BP; 2 A; 2 C; 2 G; 6 T; 0 other;\n"
        

    def test_sequence_body(self):
        converter = convert.Convert()
        assert converter.sequence_body("tctgacaatcgctttctt") == """\
     tctgacaatc gctttctt                                                      18
"""
        assert converter.sequence_body("TTTAAAACCCCGGGtttcccgggaaa") == """\
     tttaaaaccc cgggtttccc gggaaa                                             26
"""
        assert converter.sequence_body("tctgacaatcgctttctttaaaaagaaactattgtcgagaatttgcattagcaatatcactttgtcaaaaagatgtttgaatgttaaataaacattcaaaactgaatacaatatgtca") == """\
     tctgacaatc gctttcttta aaaagaaact attgtcgaga atttgcatta gcaatatcac        60
     tttgtcaaaa agatgtttga atgttaaata aacattcaaa actgaataca atatgtca         118
"""

        assert converter.sequence_body("tctgacaatcgctttctttaaaaagaaactattgtcgagaatttgcattagcaatatcactttgtcaaaaagatgtttgaatgttaaataaacattcaaaactgaatacaatatgtcac") == """\
     tctgacaatc gctttcttta aaaagaaact attgtcgaga atttgcatta gcaatatcac        60
     tttgtcaaaa agatgtttga atgttaaata aacattcaaa actgaataca atatgtcac        119
"""     

        print converter.sequence_body("tctgacaatcgctttctttaaaaagaaactattgtcgagaatttgcattagcaatatcactttgtcaaaaagatgtttgaatgttaaataaacattcaaaactgaatacaatatgtcaca")
        assert converter.sequence_body("tctgacaatcgctttctttaaaaagaaactattgtcgagaatttgcattagcaatatcactttgtcaaaaagatgtttgaatgttaaataaacattcaaaactgaatacaatatgtcaca") == """\
     tctgacaatc gctttcttta aaaagaaact attgtcgaga atttgcatta gcaatatcac        60
     tttgtcaaaa agatgtttga atgttaaata aacattcaaa actgaataca atatgtcaca       120
"""  

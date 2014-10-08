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
AC * _%s
XX
PR   Project:%s
XX
DE   %s contig %d
XX
RN   [1]
RA   %s
RT   "%s"
RL   %s.
XX
RN   [2]
RA   %s
RT   "%s"
RL   %s.
XX
RN   [3]
RA   Torsten Seemann;
RT   "Prokka: rapid prokaryotic genome annotation"
RL    Bioinformatics. 2014 Jul 15;30(14):2068-9.;
XX
CC   Data release policy http://www.sanger.ac.uk/legal/#t_2
XX
FH   Key             Location/Qualifiers
FH
"""
        assert converter.blank_header() == expected_header
    
    def test_update_locus_tag(self):
      converter = convert.Convert(locus_tag = 'new_locus_tag')
      assert converter.update_locus_tag("1234_5#6_789") == "new_locus_tag_789"
    
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
          submitter_name="Jane Doe",
          submitter_title="Direct submission",
          submitter_location="Sanger"
          )
        expected_populated_header = """\
ID   XXX; XXX; circular; genomic DNA; STD; UNC; 1234 BP.
XX
AC * _PRJ123412341
XX
PR   Project:PRJ1234
XX
DE   One line description contig 1
XX
RN   [1]
RA   John Doe
RT   "My title"
RL   Unpublished.
XX
RN   [2]
RA   Jane Doe
RT   "Direct submission"
RL   Sanger.
XX
RN   [3]
RA   Torsten Seemann;
RT   "Prokka: rapid prokaryotic genome annotation"
RL    Bioinformatics. 2014 Jul 15;30(14):2068-9.;
XX
CC   Data release policy http://www.sanger.ac.uk/legal/#t_2
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
        assert converter.source_template(1234,"My organism", 5678) == """\
FT   source          1..1234
FT                   /organism="My organism"
FT                   /mol_type="genomic DNA"
FT                   /db_xref="taxon:5678"
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

    def test_construct_feature_header(self):
        converter = convert.Convert()
        assert converter.feature_header(feature_type = 'tRNA', start = 174883, end = 174959, strand = '-') == "FT   tRNA            complement(174883..174959)\n"
        assert converter.feature_header(feature_type = 'CDS', start = 163111, end = 163365, strand = '+')  == "FT   CDS             163111..163365\n"
        
    def test_construct_feature(self):
        converter = convert.Convert()
        assert converter.construct_feature(feature_type = 'tRNA', start = 174883, end = 174959, strand = '-',feature_attributes =  {'locus_tag': 'ABC123','eC_number': '12,34' }) == """\
FT   tRNA            complement(174883..174959)
FT                   /locus_tag="ABC123"
FT                   /EC_number="12"
FT                   /EC_number="34"
FT                   /transl_table=11
"""

    def test_construct_feature_locus_tag_update(self):
        converter = convert.Convert(locus_tag = 'new_locus_tag')
        assert converter.construct_feature(feature_type = 'tRNA', start = 174883, end = 174959, strand = '+',feature_attributes =  {'locus_tag': 'ABC_123' }) == """\
FT   tRNA            174883..174959
FT                   /locus_tag="new_locus_tag_123"
FT                   /transl_table=11
"""

    def test_create_db_xref_from_inference(self):
      converter = convert.Convert()
      assert converter.construct_feature_attribute(attribute_key = 'inference', attribute_value = 'similar to AA sequence:UniProtKB:Q2G282') == """\
FT                   /db_xref="UniProtKB/Swiss-Prot:Q2G282"
"""
      assert converter.construct_feature_attribute(attribute_key = 'inference', attribute_value = 'protein motif:Pfam:PF01475.13') == """\
FT                   /db_xref="PFAM:PF01475.13"
"""

      assert converter.construct_feature_attribute(attribute_key = 'inference', attribute_value = 'protein motif:CLUSTERS:PRK09462') == """\
FT                   /db_xref="CDD:PRK09462"
"""
      assert converter.construct_feature_attribute(attribute_key = 'inference', attribute_value = 'protein motif:TIGRFAMs:TIGR01327') == """\
FT                   /db_xref="TIGRFAM:TIGR01327"
"""

      assert converter.construct_feature_attribute(attribute_key = 'inference', attribute_value = 'protein motif:Cdd:COG1932') == """\
FT                   /db_xref="CDD:COG1932"
"""


    def test_construct_feature_attribute(self):
        converter = convert.Convert()
        
        # Key with a single value
        assert converter.construct_feature_attribute(attribute_key = 'locus_tag', attribute_value = 'ABC123') == "FT                   /locus_tag=\"ABC123\"\n"
        
        # Keys to ignore
        assert converter.construct_feature_attribute(attribute_key = 'ID', attribute_value = '123') == ''
        assert converter.construct_feature_attribute(attribute_key = 'protein_id', attribute_value = 'ABC123') == ''
        
        # Key to translate between GFF name and EMBL name
        assert converter.construct_feature_attribute(attribute_key = 'eC_number', attribute_value = '1.2.3.4') == "FT                   /EC_number=\"1.2.3.4\"\n"
        
        # Keys to split over multiple lines
        assert converter.construct_feature_attribute(attribute_key = 'eC_number', attribute_value = '2.4.2.-,2.4.-.-') == "FT                   /EC_number=\"2.4.2.-\"\nFT                   /EC_number=\"2.4.-.-\"\n"
        assert converter.construct_feature_attribute(attribute_key = 'inference', attribute_value = 'ab initio prediction:Prodigal:2.60,similar to AA sequence:RefSeq:YP_005742575.1') == "FT                   /inference=\"ab initio prediction:Prodigal:2.60\"\nFT                   /inference=\"similar to AA sequence:RefSeq:YP_005742575.1\"\n"

        # Only first value should be taken by default
        assert converter.construct_feature_attribute(attribute_key = 'product', attribute_value = 'product 1, product 2') == "FT                   /product=\"product 1\"\n"
        
        # Very long attributes should be split over multiple lines
        assert converter.construct_feature_attribute(attribute_key = 'product', attribute_value = 'abc efg hij klm nop qrs tuvw xyz abc efg hij klm nop qrs tuvw xyz') == """\
FT                   /product="abc efg hij klm nop qrs tuvw xyz abc efg hij klm
FT                   nop qrs tuvw xyz"
"""
        assert converter.construct_feature_attribute(attribute_key = 'product', attribute_value = 'abc efg hij klm nop qrs tuvw xyz abc efg hij klm nop qrs tuvw xyz abc efg hij klm nop qrs tuvw xyz abc_efg hij klm nop qrs tuvw xyz') == """\
FT                   /product=\"abc efg hij klm nop qrs tuvw xyz abc efg hij kl
FT                   m nop qrs tuvw xyz abc efg hij klm nop qrs tuvw xyz abc_e
FT                   fg hij klm nop qrs tuvw xyz\"
"""
        
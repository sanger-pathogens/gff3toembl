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
        assert  emblwriter.parse_and_run() == """\
ID   XXX; XXX; circular; genomic DNA; STD; PROK; 240 BP.
XX
AC   My project2401
XX
PR   Project:My project
XX
DE   My description contig 1
XX
RN   [1]
RA   John
RT   "Some title"
RL   Some journal
XX
RN   [2]
RA   Jane
RT   "My institute"
RL   UK
XX
CC   Data release policy http://www.sanger.ac.uk/legal/#t_2
XX
FH   Key             Location/Qualifiers
FH
FT   source          1..240
FT                   /organism="Organism"
FT                   /mol_type="genomic DNA"
FT                   /db_xref="taxon:1234"
FT   CDS             complement(1..210)
FT                   /product="Peroxide stress regulator PerR%2C FUR family"
FT                   /inference="ab initio prediction:Prodigal:2.60"
FT                   /inference="similar to AA sequence:RefSeq:YP_005742566.1"
FT                   /inference="similar to AA sequence:UniProtKB:Q2G282"
FT                   /inference="protein motif:CLUSTERS:PRK09462"
FT                   /inference="protein motif:Pfam:PF01475.13"
FT                   /locus_tag="8233_4#93_02128"
FT                   /gene="perR"
SQ   Sequence 240 BP; 76 A; 54 C; 36 G; 74 T; 0 other;
     tctgacaatc gctttcttta aaaagaaact attgtcgaga atttgcatta gcaatatcac        60
     tttgtcaaaa agatgtttga atgttaaata aacattcaaa actgaataca atatgtcacg       120
     ttattccgca tcttctgaag aagatgttcc gaatatatcc ttagaaagga ggtgatccag       180
     ccgcaccttc cgatacggct accttgttac gacttcaccc caatcatttg tcccaccttc       240
//
"""
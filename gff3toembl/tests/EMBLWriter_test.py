import unittest
import sys
import os
from gff3toembl.EMBLWriter import EMBLWriter

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data')

class TestEMBLWriter(unittest.TestCase):

    def compare_files(self, calculated_filename, expected_filename):
        with open(expected_filename, 'r') as expected_file:
          expected_string = expected_file.read()
        with open(calculated_filename, 'r') as calculated_file:
          calculated_string = calculated_file.read()
        for calculated, expected in zip(calculated_string.split('\n'), expected_string.split('\n')):
          self.assertEqual(calculated, expected)
        self.assertEqual(len(calculated_string), len(expected_string))

    def test_single_feature(self):
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
           'single_feature.embl', None,11,  None )
        emblwriter.parse_and_run()
        self.compare_files('single_feature.embl', os.path.join(data_dir, 'expected_single_feature.embl'))
        os.remove('single_feature.embl')

    def test_single_feature_new_locus_tag(self):
        '''test that the script will convert from GFF3 to EMBL and change the locus tag'''
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
           'single_feature.embl', 'new_locus_tag', 11, None )
        emblwriter.parse_and_run()
        self.compare_files('single_feature.embl', os.path.join(data_dir, 'expected_single_feature_new_locus_tag.embl'))
        os.remove('single_feature.embl')

    def test_single_feature_translation_table(self):
        '''test that the script will convert from GFF3 to EMBL and change the locus tag'''
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
           'single_feature.embl', None, 1, None )
        emblwriter.parse_and_run()
        self.compare_files('single_feature.embl', os.path.join(data_dir, 'expected_single_feature_translation_table.embl'))
        os.remove('single_feature.embl')

    def test_large_conversion(self):
        '''test a large gff3 file converts to EMBL'''
        emblwriter = EMBLWriter(os.path.join(data_dir,'large_annotation.gff'),
           'Organism',
           1234,
           'My project',
           'My description',
           'John',
           'Some title',
           'Some journal',
           'circular',
           'PROK',
           'large_annotation.embl', None, 11, None )
        emblwriter.parse_and_run()
        self.compare_files('large_annotation.embl', os.path.join(data_dir, 'expected_large_annotation.embl'))
        os.remove('large_annotation.embl')


    def test_chromosome_list_conversion(self):
       '''test chromosome list creation'''
       emblwriter = EMBLWriter(os.path.join(data_dir,'chromosome_list.gff'),
          'Organism',
          1234,
          'ABC',
          'My description',
          'John',
          'Some title',
          'Some journal',
          'circular',
          'PROK',
          'chromosome_list.embl', None, 11, 'chromosome_list.txt' )
       emblwriter.parse_and_run()
       self.compare_files('chromosome_list.txt', os.path.join(data_dir, 'expected_chromosome_list.txt'))
       os.remove('chromosome_list.embl')
       os.remove('chromosome_list.txt')


    def test_remove_duplicate_tags(self):
       '''test remove duplicate tags '''
       emblwriter = EMBLWriter(os.path.join(data_dir,'duplicate_coords.gff'),
          'Organism',
          1234,
          'ABC',
          'My description',
          'John',
          'Some title',
          'Some journal',
          'circular',
          'PROK',
          'duplicate_coords.embl', None, 11, None )
       emblwriter.parse_and_run()
       self.compare_files('duplicate_coords.embl', os.path.join(data_dir, 'expected_duplicate_coords.embl'))
       os.remove('duplicate_coords.embl')



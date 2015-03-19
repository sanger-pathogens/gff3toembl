import sys
import subprocess
from gt import GFF3InStream
import os
import re

from gff3toembl.EMBLConverter import EMBLConverter
from gff3toembl.VisitorStream import VisitorStream

class EMBLWriter(object):

    def __init__(self, gff3_file, organism, taxonid, project, description, authors, title,  publication, genome_type, classification,  output_filename, locus_tag = None, translation_table = 11, chromosome_list = None):
        self.locus_tag          = locus_tag
        self.translation_table  = translation_table
        self.conv               = EMBLConverter(locus_tag, translation_table)
        self.gff3_file          = gff3_file
        self.organism           = organism
        self.taxonid            = taxonid
        self.project            = project
        self.authors            = authors
        self.title              = title
        self.publication        = publication
        self.genome_type        = genome_type
        self.classification     = classification
        self.output_filename    = output_filename
        self.chromosome_list    = chromosome_list
        self.fixed_gff_file     = str(self.gff3_file)+"_fixed.gff"

    def create_output_file(self, organism, taxonid, project, authors, title, publication, genome_type, classification):
        target = open(self.output_filename, 'w')
        for sequence_identifier, contig in sorted(self.conv.contigs.items()):
            contig.add_header(
              authors = authors,
              classification = classification,
              genome_type = genome_type,
              organism = organism,
              project = project,
              publication = publication,
              sequence_identifier = sequence_identifier,
              sequence_length = contig.sequence.length,
              sequence_name = sequence_identifier,
              taxon_id = taxonid,
              title = title,
            )
            target.write(contig.format())
            target.write("//\n")
        target.close()

    def create_chromosome_list(self, chromosome_list_filename, embl_filename):
        if chromosome_list_filename == None:
          return
        if not os.path.exists(embl_filename):
          return

        embl_file = open(embl_filename, 'r')
        chromosome_list_file = open(chromosome_list_filename, 'w')
        object_accessions = []

        for embl_line in embl_file.readlines():
          m = re.match("AC \* _(\w+)", embl_line)
          if m != None and m.group(1):
            object_accessions.append(m.group(1))

        for index, object_accession in enumerate(object_accessions):
          chromosome_name = str((index+1))
          # TODO make it work for more than just Bacteria
          chromosome_type = "Chromosome"
          if index > 0:
            chromosome_type = "Plasmid"
          chromosome_list_file.write(object_accession + "\t" + chromosome_name + "\t" + chromosome_type + "\n")

    def sort_and_tidy_gff_file(self):
        try:
          subprocess.check_call("gt gff3 -force -sort -retainids -tidy -o "+str(self.fixed_gff_file)+" "+str(self.gff3_file), shell=True)
        except:
          sys.exit("Failed to sort and tidy gff file with GT")

    def parse_and_run(self):
        self.sort_and_tidy_gff_file()
        ins = GFF3InStream(self.fixed_gff_file)
        vs = VisitorStream(ins, self.conv)
        try:
            while (vs.next_tree()):
                pass
        except Exception, e:
            print e
            exit(1)
        self.create_output_file(self.organism, self.taxonid, self.project, self.authors, self.title, self.publication, self.genome_type, self.classification)
        self.create_chromosome_list(self.chromosome_list, self.output_filename)
        os.remove(self.fixed_gff_file)


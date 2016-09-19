---
title: 'GFF3toEMBL: Preparing annotated assemblies for submission to EMBL'
tags:
  - bioinformatics
  - annotation
  - EMBL
authors:
 - name: Andrew J. Page
   orcid: 0000-0001-6919-6062
   affiliation: Pathogen Informatics, Wellcome Trust Sanger Institute
 - name: Ben Taylor
   affiliation: Pathogen Informatics, Wellcome Trust Sanger Institute
 - name: Jacqueline A. Keane
   orcid: 0000-0002-2021-1863
   affiliation: Pathogen Informatics, Wellcome Trust Sanger Institute
  
date: 19 Sept 2016
bibliography: paper.bib
---

# Summary
An essential part of open reproducable research in genomics is the deposition of annoated de novo assembled genomes in public archives such as EMBL/GenBank [@BLAXTER2016]. The interfaces provided by the major archives do not allow for data to be easily submitted on a large scale without substantial prior knowledge on the part of the submitter. This has lead to a situation where less than 15% of all sequenced bacteria have corresponding public assemblies. We address this by providing GFF3toEMBL, which converts the output of the most commonly used annotation tool, PROKKA [@SEEMANN2014], and converts it to a format suitable for submisison to EMBL. It is robust, fast, memory efficient and well tested, and has been used to submit more than 30% of all annoated genomes in EMBL/GenBank [@PAGE2016]. It is a small, but essential missing step in making genomic research more open and reproducable.

# References

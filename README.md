# gff3toembl
Converts gff3 files to emble files for uploading to EBI.

[![Build Status](https://travis-ci.org/sanger-pathogens/gff3toembl.svg?branch=master)](https://travis-ci.org/sanger-pathogens/gff3toembl)

NB this implements some EBI specific conventions and is not a generic conversion tool.

## Installation
- Install HomeBrew/LinuxBrew
- brew tap homebrew/science
- brew install https://raw.githubusercontent.com/andrewjpage/homebrew-science/gff3toembl/gff3toembl.rb

## Example usage
Run the following to get usage:
`gff3_to_embl -h`

An example:
```
gff3_to_embl --authors 'John' --title 'Some title' --publication 'Some journal' \
             --genome_type 'circular' --classification 'PROK' \
             --output_filename /tmp/single_feature.embl --translation_table 11 \
             Organism 1234 'My project' 'My description' gff3toembl/tests/data/single_feature.gff
```

## Tests
Run `python setup.py test`

## Known Issues
This doesn't work with some versions of Genometools on Mac OS X; it appears to work with Genometools 1.5.4

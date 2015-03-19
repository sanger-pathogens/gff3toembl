# gff3toemble
Converts gff3 files to emble files for uploading to EBI.

NB this implements some EBI specific conventions and is not a generic conversion tool.

## Installation
- install Genometools, including python binding
- make sure the library is in your `LD_LIBRARY_PATH` and your python is in `PYTHONPATH`
- `git clone <this repo>`
- `python setup.py install`

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
Run `python setup.py nosetests`

## Known Issues
Genometools doesnt currently work on OSX.

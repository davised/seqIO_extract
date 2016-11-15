# seqIO\_extract.py


# Installation

No installation required. After cloning the git repository, move the script to a directory in your PATH or add the cloned directory to your PATH. Install BioPython and its dependencies and you will be good to go. 

## Dependencies

Required 

* [python v. 2.7.8+](https://www.python.org)
* [Biopython v. 1.65+](http://biopython.org/DIST/docs/install/Installation.html#sec12)

# Usage

The purpose of this script is to extract matching sequences from FASTA, FASTQ, or Genbank formatted files. Basic syntax is:

`seqIO_extract.py input.fasta list.txt`

You can pass individual IDs on the command line as well:

`seqIO_extract.py --raw input.fasta id0 id1 id2`

Searching genbank files is slightly more complex, as you can choose the search type and output type, using either locus tags or protein IDs. Only CDS entries are examined using this script. You can also use the 'all' argument to convert Genbank to FASTA with coding sequences (protein by default, nucleotide with the --ffn flag) or contig/chromosome output with the --fna flag.

For example-

Genbank to protein FASTA (locus tags for IDs):

`seqIO_extract.py input.gbk all > output.faa`

Genbank to protein FASTA (protein IDs for IDs):

`seqIO_extract.py --outname protein_id input.gbk all > output.faa`

Genbank to nucleotide CDS FASTA (locus tags for IDs):

`seqIO_extract.py --ffn input.gbk all > output.ffn`

Genbank to nucleotide FASTA (contig/chromosome names for IDs):

`seqIO_extract.py --fna input.gbk all > output.fna`

You can also extract the upstream and downstream sequences from matching genes, for example 1kb fragments:

`seqIO_extract.py --us 1000 --ds 1000 input.gbk list.txt > output.fasta`

Lastly, the script allows for inexact matching with the --matchtype inexact flag. You can use this to supply substrings and, in effect, extract a range of sequences. For example:

`seqIO_extract.py --matchtype inexact --raw input.fasta 100 > output.fasta`

This would match 1000, 1001, 1002, 1003, 1004 ... 1100 2100 ... 9100. etc.

The --gi flag was added to extract accession, locus tag, and GI numbers from a genbank file. This is mainly to facilitate switching between use of each of the three identifiers. Once GI numbers are phased out from NCBI, this feature will be deprecated.

# Example Dataset

To be included...

# History

v1.1.2 - 2016-11-15 - Added ability to pass from STDIN/pipe.

v1.1.1 - 2016-10-13 - Added inexact matching for genbank file records. record.id is actually VERSION line of the genbank file.

v1.1.0 - 2016-10-13 - Added ability to search for LOCUS headers in genbank file to output a subset of the genbank file or fasta nucleotide data.

v1.0.0 - 2016-09-21 - First revision released to GitHub.

# Credits

* [Edward Davis](mailto:davised.dev@gmail.com)
* [Dr. Jeff Chang](mailto:changj@science.oregonstate.edu)

# License

Please see the LICENSE file for more information.

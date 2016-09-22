#!/usr/bin/env python
import sys
import argparse
import csv
import os.path

from Bio import SeqIO

version = '1.0.0'
date = 'September 21, 2016'


def extant_file(x):
    """
    'Type' for argparse - checks that file exists but does not open.
    """
    if not os.path.exists(x):
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: x does not exist
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x

class PrintVersion(argparse.Action):
    def __init__(self, nargs=0, **kwargs):
        if nargs != 0:
            raise ValueError('nargs for PrintVersion must be 0; it is just a flag.')
        super(PrintVersion, self).__init__(nargs=nargs, **kwargs)
    def __call__(self, parser, values, namespace, option_string=None):
        sys.stderr.write('{} version {}\n'.format(__file__,version))
        sys.stderr.write('Submitted {} to github\n'.format(date))
        parser.exit()

parser = argparse.ArgumentParser(description='Extract sequence based on ID matching from a FASTA or Genbank file.')
parser.add_argument('infile',help='Sequence file (FASTA or Genbank) of interest.',type=extant_file)
parser.add_argument('--raw',help='List of locus tags, not file with list of locus tags',action='store_true',default=False)
parser.add_argument('listfiles',help='File or list (if --raw is invoked) of IDs of genes of interest to extract; \'all\' will print all records in file',type=str,nargs='+')
parser.add_argument('--filetype',help='File type provided (default = auto). FASTQ parsing not recommended with this tool.',default='auto',choices=['FASTA','FASTQ','GBK','auto'])
parser.add_argument('--ffn',help='Print nucleotide sequence for each CDS of a genbank file (gbk to ffn conversion)',action='store_true',default=False)
parser.add_argument('--fna',help='Print nucleotide sequence for the source sequence of a genbank file (gbk to fna conversion, use in conjunction with \'all\' positional argument)',action='store_true',default=False)
parser.add_argument('--outname',help='Type of identifier to print to output from genbank file (default = locus_tag)',choices=['locus_tag','protein_id'],default='locus_tag',type=str)
parser.add_argument('--searchname',help='Type of identifier to search in genbank file (default = locus_tag)',choices=['locus_tag','protein_id'],default='locus_tag',type=str)
parser.add_argument('--us',help='Extract this length of sequence from the upstream region of a gene',type=int)
parser.add_argument('--ds',help='Extract this length of sequence from the downstream region of a gene',type=int)
parser.add_argument('--matchtype',help='Include exact (default) or inexact matches. Inexact matching takes longer.',choices=['exact','inexact'],default='exact',type=str)
parser.add_argument('-v','--verbose',help='Print progress messages. Helpful for debugging or double-checking output',action='store_true')
parser.add_argument('-V','--version',help='Print version message and quit.',action=PrintVersion)
args = parser.parse_args()
tags = {}
matches = {}
count = 0

def parse_opts(x, y):
    x.all = False
    if x.listfiles[0] == 'all':
        x.all = True
    else:
        if x.raw == False:
            for listfile in x.listfiles:
                with open(listfile) as l:
                    for line in l:
                        line = line.rstrip('\n')
                        y[line] = 1
        else:
            for tag in args.listfiles:
                y[tag] = 1
    if x.filetype == 'auto':
        with open(x.infile, 'r') as f:
            line = f.readline().strip()
            if line[0] == '>':
                x.filetype = 'FASTA'
            elif line[0] == '@':
                x.filetype = 'FASTQ'
            elif 'LOCUS' in line:
                x.filetype = 'GBK'
            else:
                sys.exit('Unable to determine filetype of {}. Quitting...\n'.format(x.infile))


def parse_fasta(x, y):
    k = 0
    infile = SeqIO.parse(x, y)
    for record in infile:
        k += 1
        if args.all == True or args.matchtype == 'exact' and record.id in tags:
            matches[record.id] = 1
        elif args.matchtype == 'inexact':
            for tag in tags:
                if tag in record.id:
                    matches[record.id] = tag
                    break
        if record.id in matches:
            SeqIO.write(record, sys.stdout, y)
    return k

def parse_genbank(x):
    k = 0
    infile = SeqIO.parse(x, 'genbank')
    if not args.fna:
        for record in infile:
            for feature in record.features:
                if feature.type == 'CDS':
                    k += 1
                    locus_tag = feature.qualifiers['locus_tag'][0]
                    try:
                        protein_id = feature.qualifiers['protein_id'][0]
                    except KeyError:
                        protein_id = 'NULL'
                    outname = locus_tag
                    searchname = locus_tag
                    if args.searchname == 'protein_id':
                        searchname = protein_id
                    if args.outname == 'protein_id':
                        outname = protein_id
                    if args.all == True or args.matchtype == 'exact' and searchname in tags:
                        matches[outname] = searchname
                    elif args.matchtype == 'inexact':
                        for tag in tags:
                            if tag in searchname:
                                matches[outname] = searchname
                                break
                    if outname in matches:
                        if args.us or args.ds:
                            #Print nucleotide sequence upstream and downstream of sequence
                            #This does not take into account the strand
                            if args.us:
                                start = feature.location.start
                                print '>{}_upstream_{}'.format(locus_tag,args.us)
                                print record[start-args.us:start].seq
                            print '>{}'.format(locus_tag)
                            print feature.extract(record).seq
                            if args.ds:
                                end = feature.location.end
                                print '>{}_downstream_{}'.format(locus_tag,args.ds)
                                print record[end:end+args.ds].seq
                        elif args.ffn == False:
                            #Protein output
                            try:
                                translation = feature.qualifiers['translation'][0]
                            except KeyError:
                                sys.stderr.write('No translation available for {}. Attempting to translate.\n'.format(outname))
                                try:
                                    translation = feature.extract(record).translate(table=1).seq
                                except:
                                    sys.stderr.write('Unable to translate {}.\n'.format(outname))
                                else:
                                    print '>' + outname
                                    print translation
                            else:
                                print '>' + outname
                                print translation
                        else:
                            #Nucleotide output for coding sequences
                            print '>' + outname
                            print feature.extract(record).seq
    else:
        #All nucleotide output (contigs/chromosomes)
        for record in infile:
            k += 1
            if args.all == True:
                SeqIO.write(record, sys.stdout, 'fasta')
            else:
                sys.exit('No input set to print to stdout. Use the \'all\' positional argument to print all of the nucleotide fasta formatted data from a genbank file')
    return k

parse_opts(args,tags)

if args.verbose:
    if not args.all:
        sys.stderr.write('Looking for {} sequences matches from {}.\n'.format(len(tags), args.infile))
    else:
        sys.stderr.write('Printing all sequences from {}.\n'.format(args.infile))

if args.filetype == 'FASTA':
    count = parse_fasta(args.infile, 'fasta')
elif args.filetype == 'GBK':
    count = parse_genbank(args.infile)
elif args.filetype == 'FASTQ':
    count = parse_fasta(args.infile, 'fastq')

if args.verbose:
    if not args.all:
        sys.stderr.write('Found {}/{} {} matches out of {} in {}.\n'.format(len(matches),len(tags),args.matchtype,count,args.infile))
    else:
        sys.stderr.write('Printed {} sequences from {}.\n'.format(count,args.infile))

#!/usr/bin/env python
from __future__ import print_function
import sys
import argparse
import csv
import os.path
import logging
from signal import signal, SIGPIPE, SIGINT, SIG_DFL
signal(SIGPIPE, SIG_DFL)
signal(SIGINT, SIG_DFL)

from Bio import SeqIO

version = '1.6.0'
date = 'April 18, 2017'

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
        eprint('{} version {}'.format(__file__,version))
        eprint('Submitted {} to github'.format(date))
        parser.exit()

def init_logger(args):
    logger = logging.getLogger(__name__)
    ch = logging.StreamHandler()
    if args.debug:
        logger.setLevel(logging.DEBUG)
        ch.setLevel(logging.DEBUG)
    elif args.verbose:
        logger.setLevel(logging.INFO)
        ch.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)
        ch.setLevel(logging.WARNING)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    return(logger)

def run_argparse():
    parser = argparse.ArgumentParser(description='Extract sequence based on ID matching from a FASTA or Genbank file.')
    parser.add_argument('infile',help='Sequence file (FASTA or Genbank) of interest.',type=extant_file)
    parser.add_argument('--raw',help='List of locus tags, not file with list of locus tags',action='store_true',default=False)
    parser.add_argument('listfiles',help='File or list (if --raw is invoked) of IDs of genes of interest to extract; \'all\' will print all records in file; - will take either filenames or list from STDIN',type=str,nargs='+')
    parser.add_argument('--filetype',help='File type provided (default = auto). FASTQ parsing not recommended with this tool.',default='auto',choices=['FASTA','FASTQ','GBK','auto'])
    parser.add_argument('--ffn',help='Print nucleotide sequence for each CDS of a genbank file (gbk to ffn conversion)',action='store_true',default=False)
    parser.add_argument('--fna',help='Print nucleotide sequence for the source sequence of a genbank file (can select records using the VERSION ID as well) (gbk to fna conversion, use in conjunction with \'all\' positional argument)',action='store_true',default=False)
    parser.add_argument('--gbk',help='Print genbank output from genbank file. Cannot convert FASTA to genbank format.',action='store_true',default=False)
    parser.add_argument('--outname',help='Type of identifier to print to output from genbank file (default = locus_tag)',choices=['locus_tag','protein_id','gene'],default='locus_tag',type=str)
    parser.add_argument('--searchname',help='Type of identifier to search in genbank file (default = locus_tag)',choices=['locus_tag','protein_id','gene'],default='locus_tag',type=str)
    parser.add_argument('--us',help='Extract this length of sequence from the upstream region of a gene',type=int)
    parser.add_argument('--ds',help='Extract this length of sequence from the downstream region of a gene',type=int)
    parser.add_argument('--matchtype',help='Include exact (default) or inexact matches. Inexact matching takes longer.',choices=['exact','inexact'],default='exact',type=str)
    parser.add_argument('-v','--verbose',help='Print progress messages.',action='store_true')
    parser.add_argument('--debug',help='Print debugging messages.',action='store_true')
    parser.add_argument('--desc',help='Print description in output - product name from genbank (if available).',action='store_true')
    parser.add_argument('--nodesc',help='Print only sequence ID in output, not description (for FASTA files).',action='store_true')
    parser.add_argument('-V','--version',help='Print version message and quit.',action=PrintVersion)
    parser.add_argument('--gi',help='Print GI number for matches.',action='store_true')
    args = parser.parse_args()
    args.logger = init_logger(args)
    args.base = os.path.basename(args.infile)
    return args

def parse_opts(args, tags):
    args.all = False
    if args.listfiles[0] == 'all':
        args.all = True
    else:
        if args.listfiles[0] == '-':
            #Read from STDIN
            args.logger.debug('Getting search terms from STDIN, "-" given as input')
            args.listfiles = [x.strip() for x in sys.stdin.read().split()]
            args.listfiles = filter(None, args.listfiles)
        if args.raw == False:
            #Read from file(s)
            args.logger.debug('Reading search terms from files:')
            for listfile in args.listfiles:
                args.logger.debug(' + {}'.format(listfile))
                with open(listfile) as l:
                    for line in l:
                        line = line.rstrip('\n')
                        tags.add(line)
            args.logger.debug('Finished reading search terms.')
        else:
            #Take raw tags from cmd line
            args.logger.debug('Taking raw values from command line as search terms')
            for tag in args.listfiles:
                tags.add(line)
    if args.filetype == 'auto':
        #Detect input file type
        args.logger.debug('Automatically detecting file type.')
        with open(args.infile, 'r') as f:
            line = f.readline().strip()
            if line[0] == '>':
                args.logger.debug('FASTA format detected.')
                args.filetype = 'FASTA'
            elif line[0] == '@':
                args.logger.debug('FASTQ format detected.')
                args.filetype = 'FASTQ'
            elif 'LOCUS' in line:
                args.logger.debug('Genbank format detected.')
                args.filetype = 'GBK'
            else:
                args.logger.critical('Unable to determine filetype of {}. Exiting\n'.format(args.infile))
                sys.exit()
    return(args, tags)


def parse_fasta(args, matches, tags):
    k = 0
    args.logger.debug('Starting search of file {}.'.format(args.base))
    infile = SeqIO.parse(args.infile, args.filetype.lower())
    for record in infile:
        if args.all == True or args.matchtype == 'exact' and record.id in tags:
            matches[record.id] = 1
        elif args.matchtype == 'inexact':
            for tag in tags:
                if tag in record.id:
                    matches[record.id] = tag
                    break
        if record.id in matches:
            if args.nodesc:
                print('>' + record.id)
                print(record.seq)
            else:
                SeqIO.write(record, sys.stdout, args.filetype.lower())
        k += 1
    args.logger.debug('Finished searching file {}.'.format(args.base))
    return k

def parse_genbank(args, matches, tags):
    j = 0
    k = 0
    args.logger.debug('Starting search of file {}.'.format(args.base))
    infile = SeqIO.parse(args.infile, 'genbank')
    if not args.fna and not args.gbk:
        for record in infile:
            for feature in record.features:
                if feature.type == 'CDS':
                    j += 1
                    try:
                        locus_tag = feature.qualifiers['locus_tag'][0]
                    except KeyError:
                        locus_tag = 'NULL'
                        if args.searchname == 'locus_tag' or args.outname == 'locus_tag':
                            args.logger.warning('No locus tag found for gene {}'.format(j))
                            args.logger.warning('Change both searchname and outname to include this entry.')
                            args.logger.warning('Skipping entry for now.')
                            continue
                    #Skip pseudo genes
                    if 'pseudo' in feature.qualifiers:
                        args.logger.warning('{} is a pseudogene. Skipping.'.format(locus_tag))
                        continue
                    #Get possible output and search terms
                    try:
                        protein_id = feature.qualifiers['protein_id'][0]
                    except KeyError:
                        protein_id = 'NULL'
                    try:
                        product = feature.qualifiers['product'][0]
                    except KeyError:
                        product = 'NULL'
                    try:
                        gene = feature.qualifiers['gene'][0]
                    except KeyError:
                        gene = 'NULL'
                    #Set search and output names
                    outname = locus_tag
                    searchname = locus_tag
                    if args.searchname == 'protein_id':
                        searchname = protein_id
                    elif args.searchname == 'gene':
                        searchname = gene
                    if args.outname == 'protein_id':
                        outname = protein_id
                    elif args.outname == 'gene' and gene != 'NULL':
                        outname = gene
                    if args.desc and product != 'NULL':
                        outname = outname + ' {}'.format(product)
                    if args.all == True or args.matchtype == 'exact' and searchname in tags:
                        matches[outname] = searchname
                    elif args.matchtype == 'inexact':
                        for tag in tags:
                            if tag in searchname:
                                matches[outname] = searchname
                                break
                    if outname in matches:
                        if args.gi:
                            gi = ''
                            for item in feature.qualifiers['db_xref']:
                                if 'GI:' in item:
                                    gi = item.split(':')[1]
                                    gi = gi.replace('"','')
                            if gi:
                                print('\t'.join([protein_id,locus_tag,gi]))
                            continue
                        if args.us or args.ds:
                            #Print nucleotide sequence upstream and downstream of sequence
                            #This does not take into account the strand
                            if args.us:
                                start = feature.location.start
                                print('>{}_upstream_{}'.format(locus_tag,args.us))
                                print(record[start-args.us:start].seq)
                            print('>{}'.format(locus_tag))
                            print(feature.extract(record).seq)
                            if args.ds:
                                end = feature.location.end
                                print('>{}_downstream_{}'.format(locus_tag,args.ds))
                                print(record[end:end+args.ds].seq)
                        elif args.ffn == False:
                            #Protein output
                            try:
                                translation = feature.qualifiers['translation'][0]
                            except KeyError:
                                eprint('No translation available for {}. Attempting to translate.'.format(outname))
                                try:
                                    translation = feature.extract(record).seq.translate(table=1)
                                except:
                                    eprint('Unable to translate {}.'.format(outname))
                                    continue
                            print('>' + outname)
                            print(translation)
                        else:
                            #Nucleotide output for coding sequences
                            print('>' + outname)
                            print(feature.extract(record).seq)
                    k += 1
    else:
        #All nucleotide output (contigs/chromosomes)
        outfmt = 'fasta'
        if args.gbk:
            outfmt = 'genbank'

        for record in infile:
            if args.all == True or args.matchtype == 'exact' and record.id in tags:
                matches[record.id] = record.id
            elif args.matchtype == 'inexact':
                for tag in tags:
                    if tag in record.id:
                        matches[record.id] = record.id
                        break
            if record.id in matches:
                if args.nodesc and outfmt == 'fasta':
                    print('>' + record.id)
                    print(record.seq)
                else:
                    SeqIO.write(record, sys.stdout, outfmt)
            k += 1
    args.logger.debug('Finished searching file {}.'.format(args.base))
    return k



def main():
    #Parse the arguments
    args = run_argparse()
    tags = set()

    args, tags = parse_opts(args,tags)

    matches = {}
    count = 0
    if not args.all:
        args.logger.info('Looking for {} sequences matches from {}.'.format(len(tags), args.base))
    else:
        args.logger.info('Printing all sequences from {}.'.format(args.base))

    #Parse the files
    if args.filetype == 'FASTA' or args.filetype == 'FASTQ':
        count = parse_fasta(args, matches, tags)
    elif args.filetype == 'GBK':
        count = parse_genbank(args, matches, tags)
    else:
        args.logger.critical('No filetype given. Should never occur. Exiting...')
        sys.exit()

    #Print status messages
    if not args.all:
        args.logger.info('Found {}/{} {} matches out of {} in {}.'.format(len(matches),len(tags),args.matchtype,count,args.base))
        seqs = [x for x in tags if x not in matches.keys()]
        if len(seqs) > 0:
            args.logger.info('Missing these sequences:')
            for seq in seqs:
                args.logger.info('{}'.format(seq))
        else:
            args.logger.info('All target sequences were found!')
    else:
        args.logger.info('Printed {} sequences from {}.'.format(count,args.base))

if __name__ == '__main__':
    main()

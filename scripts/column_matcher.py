#!/usr/bin/python3
__author__ = 'Przemek Decewicz'
__date__ = '06.01.2023'
__version__= '0.1.3'

from argparse import ArgumentParser
from collections import OrderedDict
from sys import argv

def read_file(filepath, delim, col, verbose):
    """
    Reads a file into a dict
    :param filepath: path to file
    :param delim: a character for splitting data into columns
    :param col: a column number to use as main key
    :param verbose: whether to print what was considered as key
    :return dict:
    """

    d = OrderedDict()
    with open(filepath) as inf:
        header = inf.readline().strip().split(delim)
        if header[0].startswith("#"):
            header[0] = header[0][1:]
        d['HEADER'] = header
        for line in inf:
            line = line.split(delim)
            line[-1] = line[-1].strip()
            if line[col] not in d:
                d[line[col]] = [line]
            else:
                d[line[col]].append(line)
            if verbose: print(f'  - Added {line[col]}.')

    return d

def main():
    args = ArgumentParser(
        description = 'Simply matches data from two files based on selected columns. REQUIRES FURTHER CLEANING.',
        epilog=f"Version: {__version__} ({__date__}).")

    args.add_argument('-r', '--reference',
                      type = str,
                      help = 'Path to input reference file to which data from --subject file.',
                      required = True)

    args.add_argument('-s', '--subject',
                      nargs = '+',
                      type = str,
                      help = 'Path to input file(s) with extra metadata to add to --reference file.',
                      required = True)

    args.add_argument('-m', '--many_subjects',
                      choices = ('unique', 'different'),
                      default = 'different',
                      help = 'Choose how to use multiple subjects file provided. "unique" - collect the information from all --subject files and then update --reference or "different" - add data from each --subject file as new columns. [Default: %(default)s]')

    args.add_argument('-o', '--outfile',
                      type = str,
                      help = 'Path to output file.',
                      required = True)

    args.add_argument('-d', '--ref_del',
                      type = str,
                      default = '\t',
                      help = 'Character(s) to use as delimiter in --reference file. [Default: %(default)s]')

    args.add_argument('-D', '--sub_del',
                      type = str,
                      default = '\t',
                      help = 'Character(s) to use as delimiter in --subject file. [Default: %(default)s]')

    args.add_argument('-x', '--out_del',
                      type = str,
                      default = '\t',
                      help = 'Character(s) to use as delimiter in --output file. [Default: %(default)s]')

    args.add_argument('-c', '--ref_col',
                      type = int,
                      help = 'Number of column from --reference file used for matching with --subject file.',
                      required = True)

    args.add_argument('-C', '--sub_col',
                      type = int,
                      help = 'Number of column from --subject file used fot matching with --reference.')

    args.add_argument('--skip_missing',
                      action = 'store_true',
                      help = 'If set, when the --subject file does not contain an according entry, empty cells will be added rather than breaking the run.')

    args.add_argument('--add_us',
                      action = 'store_true',
                      help = 'If set, adds unmatched subjects when there is no reference entry that could be matched.')
    
    args.add_argument('--header_prefix',
                      type = str,
                      default = '',
                      help = 'If set, adds a prefix to the header of the output file. [Default: %(default)s]')

    args.add_argument('--verbose',
                      action = 'store_true',
                      help = 'If set, all matches will be printed.')

    if len(argv[1:]) == 0:
        args.print_usage()
        args.exit()

    try:
        args = args.parse_args()
    except:
        args.exit()

    # Fix column numbers
    ref_col = args.ref_col - 1
    sub_col = args.sub_col - 1

    # Read both input files
    sub = {}

    print(f"Reading reference file")
    ref = read_file(args.reference, args.ref_del, ref_col, args.verbose)
    print(f'Read {len(ref)} lines in reference file.')

    if args.many_subjects == "different":
        print(f"-- Using \"different\" mode --")
        ref_length = len(ref['HEADER'])
        new_header = ref.pop('HEADER')
        for infile in args.subject:
            print(f'Reading {infile} file.')
            sub = read_file(infile, args.sub_del, sub_col, args.verbose)
            sub_length = len(sub['HEADER'])
            print(f'Read {len(sub)} lines in subject file.')

            # Write output file
            ref_unmatched = []
            sub_matched = []
            new_header += sub.pop('HEADER')
            for reference, lines in ref.items():
                for i,line in enumerate(lines):
                    try:
                        ref[reference][i] += sub[reference][0]
                        sub_matched.append(reference)
                    except KeyError as e:
                        if args.skip_missing:
                            line += ['']*sub_length
                            ref_unmatched.append(reference)
                        else:
                            print(f"Tried to match for '{reference}'")
                            print(f"Line {i}: {line[:ref_col]}")
                            if reference not in sub: 
                                print(f"Reference key not in subject file")
                                print(f"These are first two keys in subject: " + ", ".join(list(sub.keys())[:2]))
                            raise e
                if args.verbose: print(f'  Matched {reference}.')
            if args.add_us:
                not_su = 0
                if args.verbose: print(f"Adding umatched subject records.")
                for k, v in sub.items():
                    if v not in sub_matched:
                        ref[k] = ['']*ref_length + v
                        not_su += 1
            print('Done!\n')
            print(f"Unmatched {len(ref_unmatched)} entries for reference file.")
            print(f"Matched {len(sub_matched)} times with {len(sub)} subject entries.")
    elif args.many_subjects == "unique":
        print("Working on it :)")

    print("Writing output file")
    with open(args.outfile, 'w') as outf:
        outf.write(args.header_prefix + args.out_del.join(new_header) + "\n")
        for lines in ref.values():
            for l in lines:
                outf.write(args.out_del.join(l) + "\n")
    # print(sub_matched[:10])

if __name__ == '__main__':
    main()

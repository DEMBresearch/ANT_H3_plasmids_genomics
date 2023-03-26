#!/usr/bin/python3
__author__ = 'Przemek Decewicz'
__date__ = '2022.01.31'
__version__ = '0.1.9'

from Bio import SeqIO
from Bio import BiopythonParserWarning
from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor
from glob import glob
from os import makedirs, path
from re import sub, DOTALL
from sys import argv

import binascii
import gzip
import warnings

def is_gzip_file(f):
    """
    This is an elegant solution to test whether a file is gzipped by reading the first two characters.
    I also use a version of this in fastq_pair if you want a C version :)
    See https://stackoverflow.com/questions/3703276/how-to-tell-if-a-file-is-gzip-compressed for inspiration
    :param f: the file to test
    :return: True if the file is gzip compressed else false
    """
    with open(f, 'rb') as i:
        return binascii.hexlify(i.read(2)) == b'1f8b'

def prepare_seq(header_components, delimiter, sequence):
    """
    Prepares a sequence to be written.
    :param header_components: a list of header elements
    :param delimiter: a string used for joining a list of header elemnts
    :param sequence: any sequence to be written
    :return: a string representation of a single sequence record
    """

    # TODO consider replacing illegal characters and shortening of the sequence header
    # illegal = r'[@\'/\\`]'
    # header = sub(illegal, '.', '>' + delimiter.join(header_components))
    # header = sub(r'(.{60})', '\\1\n', header, 0, DOTALL)
    header = '>' + delimiter.join(header_components)

    return '%s\n%s\n' % (
        header,
        sub(r'(.{60})', '\\1\n', str(sequence), 0, DOTALL).strip()
        )

def extract(infiles, outdir, get_file_name, delimiter, write_record, write_proteins, write_genes, write_prot_nucl, write_intergenic, add_cds_record_num, add_coords, underscore, qualifiers, add_description, add_annotation, override_genes, skip_done):

    if 'record.id' in qualifiers:
        get_record_id = True
        qualifiers.remove('record.id')
    else: 
        get_record_id = False


    if add_annotation:
        add_annotation = add_annotation.split(',') if ',' in add_annotation else [add_annotation]

    file_cnt = 0

    for gb_file in infiles:
        added_desc = False
        indir = path.split(gb_file)[0]
        gb_file = path.basename(gb_file)
        gb_file_base = gb_file.rsplit('.', 1)[0]
        file_cnt += 1
        
        if get_file_name: 
            gb_name = delimiter.join(gb_file.rsplit('.', 1)[0]) #.split('_', 1)


        if skip_done:
            ctrl = []

            if write_record:
                ctrl.append(1) if path.isfile(path.join(outdir, 'records', gb_file_base + '.fna')) else ctrl.append(0)
            if write_proteins:
                ctrl.append(1) if path.isfile(path.join(outdir, 'prot', gb_file_base + '.faa')) else ctrl.append(0)
            if write_prot_nucl:
                ctrl.append(1) if path.isfile(path.join(outdir, 'prot_nucl', gb_file_base + '.fna')) else ctrl.append(0)
            if write_intergenic:
                ctrl.append(1) if path.isfile(path.join(outdir, 'intergenic', gb_file_base + '.fna')) else ctrl.append(0)
            if write_genes:
                ctrl.append(1) if path.isfile(path.join(outdir, 'genes', gb_file_base + '.fna')) else ctrl.append(0)

            if sum(ctrl) == len(ctrl): continue


        if write_record:
            nucl_out = open(path.join(outdir, 'records', gb_file_base + '.fna'), 'w')
        if write_proteins:
            prot_out = open(path.join(outdir, 'prot', gb_file_base + '.faa'), 'w')
        if write_prot_nucl:
            prot_nucl_out = open(path.join(outdir, 'prot_nucl', gb_file_base + '.fna'), 'w')
        if write_intergenic:
            int_out = open(path.join(outdir, 'intergenic', gb_file_base + '.fna'), 'w')
            int_regions = []
            last = 0
        if write_genes:
            genes_out = open(path.join(outdir, 'genes', gb_file_base + '.fna'), 'w')

        file_format = 'genbank'
        if is_gzip_file(path.join(indir, gb_file)):
            handle = gzip.open(path.join(indir, gb_file), 'rt')
        else:
            handle = open(path.join(indir, gb_file), 'r')
        # try:
        #     line = inf.readline()
        #     if 'LOCUS' in line[:10]:
        #         file_format = 'genbank'
        #     elif '{' in line[0]:
        #         file_format = 'ASN.1'
        #         print('REPEAT: %s' % gb_file)
        #     elif line.startswith("ID"):
        #         file_format = 'embl'
        # except UnicodeDecodeError as e:
        #     print('\n\nerror\n', gb_file, '\n')
        #     raise e
        
        #SeqIO.convert(indir + gb_file, file_format, nucl_out, 'fasta')



        prot_cnt = 1
        desc_header_components = []
        try:
            for record in SeqIO.parse(handle, file_format):

                cds_record_num = 0

                if write_intergenic:
                    for f in record.features[::-1]:
                        if f.type == 'gene':
                            last_gene_end = f.location.end
                            break
                for f in record.features:
                    header_components = [] + desc_header_components

                    if add_annotation:
                        # if 'taxonomy' in add_annotation and 'taxonomy' in record.annotations:
                        #     desc_header_components.append(','.join(record.annotations['taxonomy']))

                        if f.type == 'source' and not added_desc:
                            for desc in add_annotation:
                                if desc in f.qualifiers:
                                    if type(f.qualifiers[desc]) == str:
                                        desc_header_components.append(f.qualifiers[desc])
                                    else:
                                        desc_header_components.append(','.join(f.qualifiers[desc]))
                                else:
                                    desc_header_components.append('missing %s' % desc)
                            added_desc = True
                            header_components += desc_header_components
                    

                    
                    if f.type == 'gene' and write_genes and not override_genes:
                        if get_file_name:
                            header_components.append(path.basename(gb_file).rsplit('.', 1)[0])

                        if get_record_id:
                            header_components.append(record.id)

                        if add_description:
                            header_components.append(record.description)

                        for q in ['locus_tag', 'note']:
                            if q in f.qualifiers:
                                header_components.append(f.qualifiers[q][0])
                                continue
                            else:
                                header_components.append('N/A')
                        if add_coords:
                            header_components.append('%s..%s(%i)' % (str(f.location.start), str(f.location.end), f.location.strand))

                        genes_out.write(prepare_seq(header_components, delimiter, f.extract(record).seq))

                    if f.type == 'CDS' and write_proteins:

                        cds_record_num += 1

                        if get_file_name:
                            header_components.append(path.basename(gb_file).rsplit('.', 1)[0])

                        if get_record_id:
                            header_components.append(record.id)

                        if add_description:
                            header_components.append(record.description)

                        if add_cds_record_num:
                            header_components.append(f"cds_{cds_record_num}")

                        for q in qualifiers:
                            if q in f.qualifiers:
                                header_components.append(f.qualifiers[q][0])
                                continue
                            elif q == 'protein_id':
                                header_components.append('prot_%05d' % prot_cnt)
                                prot_cnt += 1
                            #elif q == 'product':
                            #    header_components.append('N/A')
                            else:
                                header_components.append('N/A')

                        if write_genes and override_genes:
                            gene_components = [c for c in header_components]
                            for q in ['locus_tag', 'note']:
                                if q in f.qualifiers:
                                    gene_components.append(f.qualifiers[q][0])
                                    continue
                                else:
                                    gene_components.append('N/A')
                            if add_coords:
                                gene_components.append('%s..%s(%i)' % (str(f.location.start), str(f.location.end), f.location.strand))

                            genes_out.write(prepare_seq(gene_components, delimiter, f.extract(record).seq))

                        if add_coords:
                            header_components.append('%s..%s(%i)' % (str(f.location.start), str(f.location.end), f.location.strand))

                        try:
                            seq = f.extract(record).seq
                        except ValueError:
                            print(delimiter.join(header_components))
                            continue
                        if 'translation' not in f.qualifiers:
                            aa_seq = seq.translate(table=11, to_stop=True)
                        else:
                            aa_seq = f.qualifiers['translation'][0]
                                
                        if len(aa_seq) == 0: continue

                        if not underscore:
                            prot_out.write(prepare_seq(header_components, delimiter, aa_seq))
                            if write_prot_nucl:
                                prot_nucl_out.write(prepare_seq(header_components, delimiter, seq))
                        else:
                            #print(f.qualifiers['note'][0])
                            #print(gb_name)
                            uhc = [sub(' ', underscore, x) for x in header_components]
                            prot_out.write(prepare_seq(uhc, delimiter, aa_seq))
                            if write_prot_nucl:
                                prot_nucl_out.write(prepare_seq(uhc, delimiter, seq))
                    
                    if write_intergenic and (f.type == 'CDS' or f.type == 'tRNA' or f.type == 'rRNA' or f.type == 'tmRNA' or f.type == 'ncRNA'): # f.type == 'gene' or 
                        # if f.location.start > 0:
                        int_regions.append((int(last), int(f.location.start)))
                        last = f.location.end
                        if last == last_gene_end:
                            int_regions.append((int(last), int(len(record.seq))))



                if write_intergenic:
                    header_components = []
                    if get_file_name:
                        header_components.append(gb_file.rsplit('.', 1)[0])

                    if get_record_id:
                        header_components.append(record.id)

                    if add_description:
                        header_components.append(record.description)

                    ir_cnt = 0
                    for ir in int_regions:
                        if ir[1] - ir[0] > 0: ir_cnt += 1
                        # if ir[1] - ir[0] < 0:
                        #     print(record.id, ir, ir[1] - ir[0]) 
                        if ir[1] - ir[0] >= write_intergenic:
                            # print(ir)
                            ihc = header_components + ['ir_%05d' % ir_cnt] + ['%s..%s' % (ir[0], ir[1])] if add_coords else [] 
                            header = '|'.join()
                            ir_seq = sub(r'(.{60})', '\\1\n', str(record.seq[ir[0]:ir[1]]), 0, DOTALL)
                            if not underscore:
                                int_out.write(prepare_seq(ihc, delimiter, ir_seq))
                            else:
                                uihc = [sub(' ', underscore, x) for x in ihc]
                                int_out.write(prepare_seq(uihc, delimiter, ir_seq))
                    last = 0
                    int_regions = []

                if write_record:
                    if get_file_name:
                        record_header = [gb_file.rsplit('.', 1)[0], record.id]
                        if underscore: 
                            record_header = [sub(' ', underscore, x) for x in record_header]
                        nucl_out.write(prepare_seq(record_header, delimiter, record.seq))
                    else:
                        nucl_out.write(prepare_seq([record.id], delimiter, record.seq)    )
        except ValueError as e:
            print('ValueError fpr %s' % gb_file)
            raise e
            
        if write_record:
            nucl_out.close()
        if write_proteins:
            prot_out.close()
        if write_prot_nucl:
            prot_nucl_out.close()
        if write_genes:
            genes_out.close()
        if write_intergenic:
            int_out.close()
        
        # print('Processed %3d/%3d files (%s).' % (file_cnt, len(infiles), gb_file))

    return True
    
    
def multi_extract(args):

    infiles, outdir, get_file_name, delimiter, record, proteins, genes, prot_nucl, intergenic, add_cds_record_num, add_coords, underscore, qualifiers, add_description, add_annotation, override_genes, skip_done = args
    return extract(infiles, outdir, get_file_name, delimiter, record, proteins, genes, prot_nucl, intergenic, add_cds_record_num, add_coords, underscore, qualifiers, add_description, add_annotation, override_genes, skip_done)

def main():
    args = ArgumentParser(description=
                          """This program extracts nucleotide sequences of the records and aminoacid sequences of CDSs from GenBank  files located in selected directory.""",
                          epilog=f"Version: {__version__} ({__date__})."
    )

    args.add_argument('-f', '--infiles',
                      type = str,
                      nargs = '*',
                      help = 'Path to infile(s) used for sequence extraction.')

    args.add_argument('-i', '--indir',
                      type = str,
                      help = 'Path to input directory where GenBank files are located. Overrides --infiles.')

    args.add_argument('-o', '--outdir',
                      type = str,
                      help = 'Path to output directory.',
                      required = True)
                      
    args.add_argument('-u', '--underscore',
                      type = str,
                      help = 'A character used for changing whitespaces in sequence header, eg. in product. Generally, _ or ~ would work well.')    

    args.add_argument('-x', '--delimiter',
                      type = str,
                      default = '|',
                      help = 'A character or set of characters used to separate parts of sequence header. [Default: %(default)s]')
                      
    args.add_argument('-n', '--add_file_name',
                      action = 'store_true',
                      help = 'If set while file name splited on first _ will be added to sequence header with the order: gb_file name, note, protein id, product')    
                      
    args.add_argument('-pn', '--prot_nucl',
                      action = 'store_true',
                      help = 'If set nucleotide suquences of CDS features will be extracted. [Default: do not extract]')    
                      
    args.add_argument('-ig', '--intergenic',
                      type = int,
                      help = 'Set the minimum size for intergenic regions to be extracted. [Default: do not extract]')

    args.add_argument('-C', '--add_cds_record_num',
                      action = 'store_true',
                      help = 'If set number of CDS per replicon will be added to sequence header. [Default: do not add]')    

    args.add_argument('-c', '--add_coords',
                      action = 'store_true',
                      help = 'If set coordinates of each sequence will be added to sequence header. [Default: do not add]')        

    args.add_argument('-g', '--genes',
                      action = 'store_true',
                      help = 'If set genes sequences will be extracted. [Default: off]')    

    args.add_argument('-og', '--override_genes',
                      action = 'store_true',
                      help = 'If set genes features will be overriden and extracted based on CDS features. Requires --genes flag. [Default: off]')    

    args.add_argument('-p', '--proteins',
                      action = 'store_true',
                      help = 'If set protein (CDS) sequences will be extracted. [Default: off]')    

    args.add_argument('-r', '--record',
                      action = 'store_true',
                      help = 'If set record sequences will be extracted. [Default: off]')    

    args.add_argument('-a', '--add_annotation',
                      type = str,
                      help = 'If set, adds source field annotation qulifiers values to sequences headers. [Default: off]')

    args.add_argument('-d', '--add_description',
                      action = 'store_true',
                      help = 'If set description field will be added to sequences headers.')          

    args.add_argument('-q', '--qualifiers',
                      type = str,
                      help = 'Comma-delimited list of qualifiers to add to sequence headers. [Default: record.id,protein_id,product]',
                      default = 'record.id,protein_id,product')

    args.add_argument('-t', '--threads',
                      type = int,
                      help = 'If more than 1 provided, multiple files will be processed simultaneously. [Default: 1]',
                      default = 1)

    args.add_argument('-s', '--skip_done',
                     action = 'store_true',
                     help = 'If set, records already processed will be ommited.')

    args.add_argument('--show_warnings',
                     action = 'store_true',
                     help = 'By default Biopytohn\'s warning are turned off. Use the flag to show them.')

    if len(argv[1:]) == 0:
        args.print_usage()
        args.exit()

    try:
        args = args.parse_args()
    except:
        args.exit()

    if not path.isdir(args.outdir):makedirs(args.outdir)

    if args.record:
        if not path.isdir(path.join(args.outdir, 'records')): makedirs(path.join(args.outdir, 'records'))
    if args.proteins:
        if not path.isdir(path.join(args.outdir, 'prot')): makedirs(path.join(args.outdir, 'prot'))
    if args.prot_nucl: 
        if not path.isdir(path.join(args.outdir, 'prot_nucl')): makedirs(path.join(args.outdir, 'prot_nucl'))
    if args.intergenic:
        if not path.isdir(path.join(args.outdir, 'intergenic')): makedirs(path.join(args.outdir, 'intergenic'))
    if args.genes:
        if not path.isdir(path.join(args.outdir, 'genes')): makedirs(path.join(args.outdir, 'genes'))


    if not args.infiles and not args.indir:
        print('Provide either --infiles or --indir ...')
        exit(1)

    if not args.show_warnings:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BiopythonParserWarning)

    if args.indir:
        infiles = sorted(glob(path.join(args.indir, '*.gb*')))
    else:
        infiles = sorted(args.infiles)
    
    qualifiers = args.qualifiers.split(',') if ',' in args.qualifiers else [args.qualifiers]

    extract_args = [[[infile], args.outdir, args.add_file_name, args.delimiter, args.record, args.proteins, args.genes, args.prot_nucl, args.intergenic, args.add_cds_record_num, args.add_coords, args.underscore, qualifiers, args.add_description, args.add_annotation, args.override_genes, args.skip_done] for infile in infiles]

    print(f"Processing {len(infiles)} files with {args.threads} parallel processes.")
    with ProcessPoolExecutor(max_workers = args.threads) as executor:
        cnt = 0
        for status in executor.map(multi_extract, extract_args):
            if status: cnt += 1
            print(f"  Processed {cnt}/{len(extract_args)} files.", end = '\r', flush = True)
    print(f"  Processed {cnt}/{len(extract_args)} files.\nDone!\n")
    

if __name__ == '__main__':
    main()

#!/usr/bin/python3
__author__ = "Przemyslaw Decewicz"
__date__ = "2022.09.13"
__version__ = "0.1"

from argparse import ArgumentParser
from Bio import SeqIO

import os
import re
import subprocess
import sys

# def get_seq_header(header, header_parts):

#     return sub(r'[_|]', '.', '.'.join(header.split('|')[:header_parts]))[:30]

# def read_matching_karyotypes(filepath):

#     matching_karyotypes = {}
#     with open(filepath) as inf:
#         for line in inf:
#             k, n = line.strip().split('\t')
#             matching_karyotypes[get_seq_header(k, 1)] = sub(r'[_ |]', '.', n)

#     return matching_karyotypes
def info_colors(colors):
    """
    Prints ranges for each color
    :param colors: a dict of colors
    """
    
    info_str = ""
    min_val = 0
    for val, color in colors.items():
        info_str += f" [{min_val} - {val}) - {color}\n"
        min_val = val
        
    return info_str

def fix_label(string, n = '_'):
    """
    Converts whitespaces into desired character
    :param string: a string
    :param n: a character to change whitespace into
    :return string: modified string
    """
    
    return re.sub(r'\s+', n, string)

def call_blast(program, query, outfile, evalue, minqcov, minid):
    """
    Calls BLAST to get pairwise sequence similarities.
    :param query: query file
    :param outfile: output file
    :param evalue: E-value threshold
    :param minqcov: minimum query coverage
    :param minid: minimum identity
    :return: path to BLAST output
    """

    blast_cmd = f"{program} -query {query} -subject {query} -out {outfile} -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp\" -evalue {evalue} "
    if program == 'blastn':
        blast_cmd += f" -perc_identity {minid}"
    if program == 'blastp':
        blast_cmd += f" -qcov_hsp_perc {minqcov}"
    subprocess.call(blast_cmd, shell=True)

    return 

def get_coords(header, program, start, end):
    """
    Gets coordinates of a sequence from its header.
    :param header: sequence header
    :param program: BLAST program
    :param start: start position
    :param end: end position
    :return: sequence coordinates
    """

    if program == 'blastp':
        query_coords = ' '.join(header.split('|')[-1].split('(')[0].split('..'))
    else:
        query_coords = '%s %s' % (start, end)
    
    return re.sub(r'<|>', '', query_coords)

def parse_blast(blast_out, program, minid, minqcov, keep_self_nucl, colors, skip_str):
    """
    Parses BLAST output.
    :param blast_out: BLAST output file
    :param program: BLAST program
    :param minid: minimum identity
    :param minqcov: minimum query coverage
    :param keep_self_nucl: keep self hits for nucleotide sequences
    :param colors: color scheme
    :return records: list of records
    :return proteins: list of proteins
    :return hits: list of BLAST hits
    """

    res = {}
    print(f' - Reading BLAST output file: {blast_out}.')
    cnt = 0
    cnt_ident = 0
    cnt_below = 0
    cnt_rec = 0
    gcnt = 0
    records = [] # for updating karyotypes and links
    min_pident = 101
    max_pident = 0
    with open(blast_out) as inf:
        for line in inf:
            gcnt += 1
            line = line.split()
            query, subject, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qcovhsp = line
            pident = float(pident)
            qcovhsp = float(qcovhsp)

            # get query and subject names
            query_rec = query.split('|')[0]
            subject_rec = subject.split('|')[0]
            
            # in case extra filtering must be applied
            query_skip = query.split(skip_str)[0]
            subject_skip = subject.split(skip_str)[0]

            # filter
            if not keep_self_nucl and ((query == subject) or (query_rec == subject_rec) or (query_skip == subject_skip)): 
                cnt_ident += 1
                # record header names
                # if query_rec not in records: records.append(query_rec)
                if query not in records: records.append(query)
                continue
            else:
                # record header names
                # if query_rec not in records: records.append(query_rec)
                # if subject_rec not in records: records.append(subject_rec)
                if query not in records: records.append(query)
                if subject not in records: records.append(subject)
            if qcovhsp < minqcov:
                cnt_below += 1
                continue
            if pident < minid:
                cnt_below += 1
                continue            

            # get query coords
            query_coords = get_coords(query, program, qstart, qend)
            subject_coords = get_coords(subject, program, sstart, send)

            # pick the color for ribbon
            for min_col_value, color in colors.items():
                if pident <= min_col_value:
                    ribbon_color = color
                    break

            if pident < min_pident: min_pident = pident
            if pident > max_pident: max_pident = pident

            pq = pident if program == 'blastn' else pident * qcovhsp

            # if any pair identified so far, add new
            # pair = (query_rec, subject_rec)
            # pair_rev = (subject_rec, query_rec)
            pair = (query, subject)
            pair_rev = (subject, query)
            if pair not in res and pair_rev not in res:
                res[pair] = [[pq, query_coords, subject_coords, ribbon_color]]
                cnt += 1
            # if query and subject already in
            elif pair in res:
                res[pair].append([pq, query_coords, subject_coords, ribbon_color])
                cnt += 1
            else:
                cnt_rec += 1


    print(f' - Considering {cnt}/{gcnt} BLAST hits between {len(records)} sequences.')
    print(f' - % identity ranged from {min_pident} to {max_pident}')
    print(f' - Skipped:\n\t{cnt_ident} - identical hits.\n\t{cnt_below} - below {minqcov} mincov threshols\n\t{cnt_rec} - reciprocated')

    print(f' - Sorting BLAST hits.')
    filtered_blast_hits = []
    for pair, hits in res.items():
        pair = list(pair)
        filtered_blast_hits.extend([hit + pair for hit in hits])
    
    return res, records, filtered_blast_hits

def conf_housekeeping():
    
    return "<<include etc/housekeeping.conf>>\n<<include etc/colors_fonts_patterns.conf>>\n\n"

def conf_colors(my_colors):
    
    # write colors
    conf = ""
    conf += f"<colors>\n"
    for c, rgb in sorted(my_colors.items()):
        conf += f"{c} = {rgb}\n"
    conf += f"</colors>\n\n"
    
    return conf

def conf_karyotypes(ordered_records, nucl, nucl_labels_file, karyo_file):

    print(f" - Writing karyotypes file for {len(ordered_records)} sequences.")
    fasta = list(SeqIO.parse(nucl, 'fasta'))
    matching_names = {}
    for rec in fasta:
        nucl_id = rec.id #.split('|')[0]
        matching_names[nucl_id] = [rec.description, len(rec.seq)]
    if nucl_labels_file:
        nucl_labels = {}
        with open(nucl_labels_file) as inf:
            for line in inf:
                line = line.strip().split('\t')
                nucl_labels[line[0]] = fix_label(line[1])
    
    with open(karyo_file, 'w') as outf:
        for nucl_id in ordered_records:
            if nucl_labels_file:
                label = nucl_labels[nucl_id]
            else:
                label = nucl_id
            outf.write(f'chr - {nucl_id} {label} 0 {matching_names[nucl_id][1]} vvlgrey\n')
            outf.write(f'band {nucl_id} beg beg 0 {round(matching_names[nucl_id][1] * 0.05)} lgreen\n')
            outf.write(f'band {nucl_id} end end {matching_names[nucl_id][1] - round(matching_names[nucl_id][1] * 0.05)} {matching_names[nucl_id][1]} lred\n')

    return f"karyotype           = {karyo_file}\n\n"

def conf_links(seq_type, links_file, hits):
    
    print(f" - Writing {len(hits)} {seq_type} links in {links_file}.")
    z_value = 5
    
    with open(links_file, 'w') as outf:
        z_step = 95/len(hits)
        for x in sorted(hits):
            pq, query_coords, subject_coords, ribbon_color, query_name, subject_name = x
            z_value += z_step
            outf.write(f"{query_name} {query_coords} {subject_name} {subject_coords} color={ribbon_color},z={z_value},radius=0.98r\n")
    conf = ""
    # write links
    conf += f"<links>\n"
    conf += f"radius              = 0.66r\n"
    conf += f"bezier_radius       = 0.10r\n"
    conf += f"ribbon              = yes\n"
    conf += f"color               = black\n"
    conf += f"stroke_thickness    = 1\n"
#         conf += f"show_link_label     = yes\n"
#         conf += f"link_label_font     = condensed\n"
#         conf += f"link_label_size     = 10p\n"
#         conf += f"link_label_color    = black\n"
#         conf += f"link_label_parallel = yes\n"
#         conf += f"link_label_case     = upper\n"
#         conf += f"link_label_snuggle  = yes\n"
#         conf += f"link_label_radius   = dims(ideogram,radius_outer) + 0.075r\n"
    conf += f"link_thickness      = 1\n"
    conf += f"link_color          = black\n"
    conf += f"link_dims           = 0.075r,0.075r,0.075r,0.075r\n"
    conf += f"link_dims_default   = 0.075r,0.075r,0.075r,0.075r\n"
    conf += f"link_cap            = yes\n"
    conf += f"link_cap_style      = 0\n"
    conf += f"link_cap_color      = black\n"
    conf += f"link_cap_thickness  = 1\n"
    conf += f"link_cap_dims       = 0.075r,0.075r,0.075r,0.075r\n"
    conf += f"link_cap_dims_default = 0.075r,0.075r,0.075r,0.075r\n"
    conf += f"link_inverted       = no\n"
    conf += f"link_inverted_thickness = 1\n"
    conf += f"link_inverted_color = black\n"
    conf += f"link_inverted_dims  = 0.075r,0.075r,0.075r,0.075r\n"
    conf += f"link_inverted_dims_default = 0.075r,0.075r,0.075r,0.075r\n"
    conf += f"link_inverted_cap   = yes\n"
    conf += f"link_inverted_cap_style = 0\n"
    conf += f"link_inverted_cap_color = black\n"
    conf += f"link_inverted_cap_thickness = 1\n"
    conf += f"link_inverted_cap_dims = 0.075r,0.075r,0.075r,0.075r\n"
    conf += f"link_inverted_cap_dims_default = 0.075r,0.075r,0.075r,0.075r\n"
    conf += f"<link>\n"
    conf += f"show                = yes\n"
    conf += f"file                = {links_file}\n"
    conf += f"record_limit        = 100000\n"
    conf += f"record_limit_type   = total\n"
    conf += f"record_limit_action = hide\n"
    conf += f"</link>\n"
    conf += f"</links>\n"
    conf += f"\n"
        
    
    return conf
    
def conf_ideograms(i_label_size):
    
    conf = ""
    # write ideogram
    conf += "<ideogram>\n"
    conf += f"<spacing>\n"
    conf += f"default             = 0.006r\n"
    conf += f"break               = 10u\n"
    conf += f"axis_break_at_edge  = no\n"
    conf += f"axis_break          = no\n"
    conf += f"axis_break_style    = 1\n"
    conf += f"</spacing>\n"
    conf += f"thickness           = 12p\n"
    conf += f"stroke_thickness    = 1\n"
    conf += f"stroke_color        = black\n"
    conf += f"fill                = yes\n"
    conf += f"fill_color          = lgrey\n"
    conf += f"radius              = 0.66r\n"
    conf += f"show_label          = yes\n"
    conf += f"label_with_tag      = yes\n"
    conf += f"label_font          = condensed\n"
    conf += f"label_color         = grey\n"
    conf += f"label_radius        = dims(ideogram,radius) + 0.104r\n"
    conf += f"label_size          = {i_label_size}p\n"
    conf += f"band_stroke_thickness = 0\n"
    conf += f"show_bands          = yes\n"
    conf += f"fill_bands          = yes\n"
    conf += f"</ideogram>\n"
    conf += f"\n"
    
    return conf

def conf_ticks():
    
    conf = ""
    # write ticks
    conf += f"show_ticks          = yes\n"
    conf += f"show_tick_labels    = no\n"
    conf += f"<ticks>\n"
    conf += f"radius              = dims(ideogram,radius_outer)\n"
    conf += f"multiplier          = 1e-6\n"
    conf += f"<tick>\n"
    conf += f"spacing             = 100u\n"
    conf += f"rspacing            = 0.1\n"
    conf += f"spacing_type        = relative\n"
    conf += f"size                = 5p\n"
    conf += f"thickness           = 1p\n"
    conf += f"color               = vvvvdgrey\n"
    conf += f"show_label          = no\n"
    conf += f"label_size          = 10p\n"
    conf += f"label_offset        = 5p\n"
    conf += f"format              = %d\n"
    conf += f"</tick>\n"
    conf += f"</ticks>\n"
    conf += f"\n"
    
    return conf
             
def conf_chromosomes(ordered_records):
    
    conf = ""
    # write chromosomes
    chromosomes = ','.join(ordered_records) 
    conf += f"<chromosomes>\n"
    conf += f"chromosomes_units   = 1\n"
    conf += f"chromosomes_display_default = yes\n"
    conf += f"chromosomes_order   = {chromosomes}\n"
    conf += f"radius              = 0.66r\n"
    conf += f"stroke_thickness    = 1\n"
    conf += f"stroke_color        = black\n"
    conf += f"fill                = yes\n"
    conf += f"fill_color          = lgrey\n"
    conf += f"show_label          = yes\n"
    conf += f"label_font          = condensed\n"
    conf += f"label_color         = grey\n"
    conf += f"label_radius        = dims(ideogram,radius) + 0.104r\n"
    conf += f"band_stroke_thickness = 0\n"
    conf += f"show_bands          = yes\n"
    conf += f"fill_bands          = yes\n"
    conf += f"</chromosomes>\n"
    conf += f"\n"
    
    return conf

def conf_image(output, image_file, image_size):
    
    conf = ""
    # write image
    conf += f"<image>\n"
    conf += f"dir                 = {output}\n"
    conf += f"file                = {image_file}\n"
    conf += f"radius              = {image_size}p\n"
    conf += f"angle_offset        = -90\n"
    conf += f"svg                 = yes\n"
    conf += f"png                 = yes\n"
    conf += f"png_use_dynamic_colormap = yes\n"
    conf += f"24bit               = yes\n"
    conf += f"height              = {image_size + round(image_size * 0.05)}p\n"
    conf += f"width               = {image_size + round(image_size * 0.05)}p\n"
    conf += f"auto_alpha_colors   = yes\n"
    conf += f"auto_alpha_steps    = 5\n"
    conf += f"background_color    = white\n"
    conf += f"stroke_color        = black\n"
    conf += f"stroke_thickness    = 1\n"
    conf += f"raster              = yes\n"
    conf += f"raster_chromosomes  = yes\n"
    conf += f"raster_black_and_white = yes\n"
    conf += f"raster_black_and_white_threshold = 0.5\n"
    conf += f"raster_black_and_white_transparent = yes\n"
    conf += f"raster_black_and_white_transparent_threshold = 0.5\n"
    conf += f"raster_resolution   = 300\n"
    conf += f"raster_font_antialias = yes\n"
    conf += f"raster_shape_antialias = yes\n"
    conf += f"raster_image_filter = lanczos\n"
    conf += f"raster_chromosome_filter = lanczos\n"
    conf += f"raster_shape_filter = lanczos\n"
    conf += f"raster_font_filter = lanczos\n"
    conf += f"raster_font_size    = 10\n"
    conf += f"raster_font_dpi     = 300\n"
    conf += f"</image>\n"
    conf += f"\n"
             
    return conf
                                                        
def conf_highlights(prot, highlights_file):
    
    print(f" - Writing highlights file to {highlights_file}")
    with open(highlights_file, 'w') as outf:
        for prot in SeqIO.parse(prot, 'fasta'):
            nucl_id, prot_id, coords = prot.id.split('|')
            start, stop = get_coords(prot.id, 'blastp', 0, 0).split()
            strand = prot.id[:-1].split('(')[1]
            if strand == '1':
                c = "colour-rainbow-2_a2" # if not hyp else "grey_a2"
                r = f"r0=1.05r,r1=1.10r,fill_color={c}"
            else:
                c = "colour-rainbow-7_a2" # if not hyp else "grey_a2"
                r = f"r0=1.00r,r1=1.05r,fill_color={c}"
    #         if strand == '1':
    #             rl = "r0=1.10r,r1=1.15r"
    #         else:
    #             rl = "r0=0.90r,r1=0.95r"
    #         rl = "r0=1.10r,r1=1.15r"

            outf.write(f"{nucl_id}\t{start}\t{stop}\t{r}\n")
    #         if not hyp:
    #             product = line[3]
    #             outftxt.write(f"{nucl_id}\t{start}\t{stop}\t{rl}\n")
    
    conf = ""
    # write highlights
    conf += f"<highlights>\n"
    conf += f"<highlight>\n"
    conf += f"file                = {highlights_file}\n"
    conf += f"#fill_color=blue_a2\n"
    conf += f"ideogram            = no\n"
    conf += f"stroke_thickness    = 0.5\n"
    conf += f"stroke_color        = vvvgrey\n"
    conf += f"</highlight>\n"
    conf += f"</highlights>\n"
    conf += f"\n"

    return conf

def call_circos(conf_file, log_file):
    
    print(f" - Calling circos")
    cmd = f"circos -noparanoid -conf {conf_file} > {log_file}"
    print(f" - Command: {cmd}")
    subprocess.call(cmd, shell=True)

def generate_labels(table_file, output, name, seq_type):
    """
    Generates bands file with descriptions.
    :param table_file: a tab-delimited two column file with ideogram id and label for it
    :param output: output directory
    :param name: output name
    :param seq_type: sequence type
    """
    
    from re import sub
    infile = 'prots.fasta'
    highlights_file = os.path.join(output, f"{name}.{seq_type}.highlights.txt")
    textbands_file = os.path.join(output, f"{name}.{seq_type}.textbands.txt")

    with open(outfile, 'w') as outf, open(textfile, 'w') as outftxt:
        with open(infile) as inf:
            for line in inf:
                if line.startswith('>'):
                    line = line.strip().split('|')
                    genome = line[0][1:]
                    genome = sub('_', '.', genome)
    #                 product = sub('~', ' ', line[3])
                    coords, strand = line[-1][:-1].split('(')
                    coords = sub('[<>]', '', coords)
                    start, stop = coords.split('..')
                    hyp = True if line[3] == "hypothetical~protein" else False
                    if strand == '1':
                        c = "colour-rainbow-2_a2" if not hyp else "grey_a2"
                        r = f"r0=1.05r,r1=1.10r,fill_color={c}"
                    else:
                        c = "colour-rainbow-7_a2" if not hyp else "grey_a2"
                        r = f"r0=1.00r,r1=1.05r,fill_color={c}"
                    if strand == '1':
                        rl = "r0=1.10r,r1=1.15r"
                    else:
                        rl = "r0=0.90r,r1=0.95r"
                    rl = "r0=1.10r,r1=1.15r"

                    outf.write(f"{genome}\t{start}\t{stop}\t{r}\n")
                    if not hyp:
                        product = line[3]
                        outftxt.write(f"{genome}\t{start}\t{stop}\t{product}\t{rl}\n")
    

def main():
    args = ArgumentParser(
        description = "Takes sequence files and performs all against all BLAST searches to generate CIRCOS plots.",
        usage = f"%(prog)s --nucl records.fasta --name  TEST --output TESTING_CIRCOS_WRAPPER --prot proteins.fasta -pi 60 --order plasmids.nucl.labels.tsv --skip_searches --add_nucl_label plasmids.nucl.labels.tsv --add_prot_highlights",
        epilog = f"v{__version__} ({__date__}) written by {__author__}")

    # DEFAULT SETTINGS
    args.add_argument('-n', '--nucl', type = str, dest = 'nucl', required = True, help = 'Input nucleotide FASTA file. >REC_ID')
    args.add_argument('-p', '--prot', type = str, dest = 'prot', help = 'Input protein FASTA file. Requires a matching nucleotide FASTA file. >REC_ID|PROT_ID|START..STOP(1/-1)')
    # args.add_argument('-g', '--genbank', type = str, dest = 'genbank', required = True, help = 'Input GenBank file(s).', nargs='+') # not supported yet
    args.add_argument('-o', '--output', type = str, dest = 'output', help = 'Output directory. Default: CIRCOS_WRAPPER_OUTDIR', default = 'CIRCOS_WRAPPER_OUTDIR')
    args.add_argument('-N', '--name', type = str, dest = 'name', help = 'Output name. Default: %(default)s', default = 'CIRCOS')
    # args.add_argument('-m', '--mode', type = str, dest = 'mode', choices = ('nucl', 'prot'), help = 'What data are searched. Prots require headers of type >RECORD_ID|PROT_ID|1..1000(1). Default: %(default)s', default = 'nucl') # 'all' is not supported yet
    # SEARCH SETTINGS
    args_search = args.add_argument_group('Search settings')
    args_search.add_argument('-ne', '--nucl_evalue', type = float, help = 'E-value cutoff. Default: %(default)s', default = 1e-100)
    args_search.add_argument('-nq', '--nucl_minqcov', type = int, help = 'Minimum query coverage. Default: %(default)s', default = 0.0)
    args_search.add_argument('-ni', '--nucl_minid', type = int, help = 'Minimum identity. Default: %(default)s', default = 30.0)
    args_search.add_argument('-pe', '--prot_evalue', type = float, help = 'E-value cutoff. Default: %(default)s', default = 1e-10)
    args_search.add_argument('-pq', '--prot_minqcov', type = int, help = 'Minimum query coverage. Default: %(default)s', default = 75.0)
    args_search.add_argument('-pi', '--prot_minid', type = int, help = 'Minimum identity. Default: %(default)s', default = 30.0)
    args_search.add_argument('-k', '--keep_self_nucl', action = 'store_true', help = 'Keep BLAST hits to the same nucleotide record. Applies for both NUCL and PROT')
    args_search.add_argument('-s', '--skip_str_nucl', type = str, default = '|', help = 'Additionally check if hits are from the same group, i.e. based on sequence header split. Applieis for both NUCL and PROT')
    args_search.add_argument('-b', '--skip_searches', action = 'store_true', help = 'Skip BLAST searches if files are present.')
    # VISUAL SETTINGS
    args_visual = args.add_argument_group('Visual settings')
    args_visual.add_argument('--colors', type = str, help = 'Colors to use and ranges for such regarding percent identity of sequences. Default: %(default)s', default = '65:quick-silver_a2,75:blue-ncs_a2,85:lime-green_a2,95:yellow-orange_a2,100:amaranth_a2')
    args_visual.add_argument('--format', type = str, choices = ('png', 'svg'), help = 'Choose format for output file: svg or png. Default %(default)s', default = 'svg')
    args_visual.add_argument('--i_label_size', type = int, default = '12', help = 'Change ideograms label size. Default: %(default)s')
    args_visual.add_argument('--image_size', type = int, default = '800', help = 'Change image size (in pixels). Default: %(default)s')
    args_visual.add_argument('--add_prot_highlights', action = 'store_true', help = 'Adds bands for proteins.')  
    args_visual.add_argument('--add_prot_label', type = str, help = 'Tab-delimited table with protein_id and label for corresponding protein bands.')
    args_visual.add_argument('--add_nucl_label', type = str, help = 'Tab-delimited table with nucl_id and label for it.')  
    args_visual.add_argument('--order', type = str, help = 'A file with nucl_ids to order the ideograms. If not provided it will be based on identity.')
    

    if len(sys.argv[1:]) == 0:
        args.print_help()
        sys.exit()

    try:
        args = args.parse_args()
    except:
        sys.exit()

    # my colors https://coolors.co/a09e98-0a8cc7-73ca44-eacfff-fba945-dd454c
    my_colors = {
        "red": "247, 42, 66",
        "green": "51, 204, 94",
        "blue": "54, 116, 217",
        "orange": "255, 136, 0",
        "lime": "186, 255, 0",
        "colour-rainbow-1": "255, 0, 0",
        "colour-rainbow-2": "255, 127, 0",
        "colour-rainbow-3": "255, 255, 0",
        "colour-rainbow-4": "0, 255, 0",
        "colour-rainbow-5": "0, 0, 255",
        "colour-rainbow-6": "75, 0, 130",
        "colour-rainbow-7": "143, 0, 255",
        "grey-rainbow-1": "230, 230, 230",
        "grey-rainbow-2": "200, 200, 200",
        "grey-rainbow-3": "170, 170, 170",
        "grey-rainbow-4": "140, 140, 140",
        "grey-rainbow-5": "110, 110, 110",
        "grey-rainbow-6": "80, 80, 80",
        "grey-rainbow-7": "50, 50, 50",
        "quick-silver": "160, 158, 152",
        "blue-ncs": "10, 140, 199",
        "lime-green": "115, 202, 68",
        "thistle": "234, 207, 255",
        "yellow-orange": "251, 169, 69",
        "english-vermillion": "221, 69, 76",
        "amaranth": "228, 34, 76",
        "dark-red": "141, 0, 0",
    }

    # check if input files exist
    if not os.path.isfile(args.nucl):
        print('ERROR: Nucleotide FASTA file does not exist.')
        sys.exit()
    if args.prot and not os.path.isfile(args.prot):
        print('ERROR: Protein FASTA file does not exist.')
        sys.exit()

    # check if input files are not empty
    if os.stat(args.nucl).st_size == 0:
        print('ERROR: Nucleotide FASTA file is empty.')
        sys.exit()
    if args.prot and os.stat(args.prot).st_size == 0:
        print('ERROR: Protein FASTA file is empty.')
        sys.exit()

    # check if output directory exists
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    # check if colors present in our preset
    colors = {float(x.split(':')[0]): x.split(':')[1] for x in args.colors.split(',')}
    for color in colors.values():
        if color.rsplit('_', 1)[0] not in my_colors:
            print('ERROR: Color {} not found in our preset.'.format(color))
            sys.exit()
    print("Ribbon Colors correspond to the following % identity ranges:")
    print(info_colors(colors))
    
    # write run_info 
    run_info_file = os.path.join(args.output, 'run_info.txt')
    with open(run_info_file, 'w') as outf:
        outf.write("Ribbon Colors correspond to the following % identity ranges:\n")
        outf.write(info_colors(colors))
        outf.write("Params used for search are as follows:\n")
        outf.write(f"Nucleotide: e-value <= {args.nucl_evalue}, % identity >= {args.nucl_minid}, query HSP coverage >= {args.nucl_minqcov}.\n")
        if args.prot:
            outf.write(f"Proteins: e-value <= {args.prot_evalue}, % identity >= {args.prot_minid}, query HSP coverage >= {args.prot_minqcov}.\n")

    # call BLAST
    print('BLASTing nucleotide sequences...')
    bout_nucl_file = os.path.join(args.output, f"{args.name}.blastn.e{args.nucl_evalue}_p{args.nucl_minid}_q{args.nucl_minqcov}.tsv")
    if (not args.skip_searches and os.path.isfile(bout_nucl_file)) or not os.path.isfile(bout_nucl_file):
        call_blast('blastn', args.nucl, bout_nucl_file, args.nucl_evalue, args.nucl_minqcov, args.nucl_minid)
    if args.prot:
        print('BLASTing protein sequences...')
        bout_prot_file = os.path.join(args.output, f"{args.name}.blastp.e{args.prot_evalue}_p{args.prot_minid}_q{args.prot_minqcov}.tsv")
        if (not args.skip_searches and os.path.isfile(bout_prot_file)) or not os.path.isfile(bout_prot_file):
            call_blast('blastp', args.prot, bout_prot_file, args.prot_evalue, args.prot_minqcov, args.prot_minid)

    # parse BLAST output
    print('Parsing BLAST output...')
    nucl_res, nucl_records, nucl_hits = parse_blast(bout_nucl_file, 'blastn', args.nucl_minid, args.nucl_minqcov, args.keep_self_nucl, colors, args.skip_str_nucl)
    if args.prot:
        prot_res, prot_records, prot_hits = parse_blast(bout_prot_file, 'blastp', args.prot_minid, args.nucl_minqcov, args.keep_self_nucl, colors, args.skip_str_nucl)

    # create config file
    print('Creating CIRCOS config file for nucleotide sequences...')
    if args.order:
        print(f" - Reading order file")
        with open(args.order) as inf:
            order = []
            for line in inf:
                order.append(line.split()[0])
        if len(set(order)) == len(nucl_records):
            ordered_records = order
        else:
            print(f"ERROR: {len(set(order))} records ids were provided instead of {len(nucl_records)} required")
            sys.exit()
    else:
        ordered_records = nucl_records
        
    # PREPARE CONF FILES
    seq_type = 'nucl'
    karyo_file = os.path.join(args.output, f"{args.name}.{seq_type}_karyotype")
    links_file = os.path.join(args.output, f"{args.name}.{seq_type}_links")
    conf_file = os.path.join(args.output, f"{args.name}.{seq_type}_conf")
    log_file = os.path.join(args.output, f"{args.name}.{seq_type}.log")
    image_file = f"{args.name}.{seq_type}"
    conf = conf_housekeeping()
    conf += conf_colors(my_colors)
    conf += conf_karyotypes(ordered_records, args.nucl, args.add_nucl_label, karyo_file)
    conf += conf_links(seq_type, links_file, nucl_hits)
    conf += conf_ideograms(args.i_label_size)
    conf += conf_chromosomes(ordered_records)
    conf += conf_ticks()
    conf += conf_image(args.output, image_file, args.image_size)
    # write conf file
    with open(conf_file, 'w') as outf:
        outf.write(conf)
    # call circos
    print('Calling CIRCOS...')
    call_circos(conf_file, log_file)

    if args.prot:
        print('Creating CIRCOS config file for protein sequences...')
        print(f" - proteins were shared between {len(prot_records)}/{len(nucl_records)} records")
        seq_type = 'prot'
        karyo_file = os.path.join(args.output, f"{args.name}.{seq_type}_karyotype")
        links_file = os.path.join(args.output, f"{args.name}.{seq_type}_links")
        conf_file = os.path.join(args.output, f"{args.name}.{seq_type}_conf")
        log_file = os.path.join(args.output, f"{args.name}.{seq_type}.log")
        highlights_file = os.path.join(args.output, f"{args.name}.{seq_type}_highlights")
        image_file = f"{args.name}.{seq_type}"
        conf = conf_housekeeping()
        conf += conf_colors(my_colors)
        conf += conf_karyotypes(ordered_records, args.nucl, args.add_nucl_label, karyo_file)
        conf += conf_links(seq_type, links_file, prot_hits)
        conf += conf_ideograms(args.i_label_size)
        conf += conf_chromosomes(ordered_records)
        conf += conf_ticks()
        conf += conf_image(args.output, image_file, args.image_size)
        if args.add_prot_highlights:
            conf += conf_highlights(args.prot, highlights_file)
            
        # write conf file
        with open(conf_file, 'w') as outf:
            outf.write(conf)
        # call circos
        print('Calling CIRCOS...')
        call_circos(conf_file, log_file)

    print(f'Done!')


if __name__ == '__main__':
    main()


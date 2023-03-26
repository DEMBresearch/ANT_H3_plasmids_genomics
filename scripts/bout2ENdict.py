#!/usr/bin/python3
__author__ = 'Przemek Decewicz'
__date__ = '17.06.2021'

from argparse import ArgumentParser
from os import makedirs, path
import re
from sys import argv

def reformat(infile, fields, grouping, evalue, qcovhsp, pid, identical, full_cure, get_full_header, add_pq, consider_once):
    
    nodes       = set() #TODO rewrite to dictionary: nodes should not 
    edges       = {} # all edges based on parts of sequence headers, eventually full sequence headers
    fcure       = {} # all edges based on full sequence headers
    node_groups = {}
    node_once   = {}

    print('\nReading BLAST output file...')

    # get bout res number 
    file_length = 0
    with open(infile) as inf:
        for line in inf:
            file_length += 1

    print('The file %s contains %i BLAST hits.' % (infile, file_length))

    cnt = 0
    with open(infile) as inf:
        line = inf.readline().strip().split('\t')
        while line != ['']:
            cnt += 1
            if cnt % 10000 == 0: 
                print('  Read  %0.2f%% (%i / %i) BLAST hits.' % (cnt/file_length * 100, cnt, file_length), end='\r', flush=True)

            if not grouping: # consider whole sequence headers
                node0 = line[0]
                node1 = line[1]
                nodes.add(node0)
                nodes.add(node1)
            else: # split header and pick the column
                node0_header = line[0].split('|')
                node0 = node0_header[grouping - 1]
                node0_ext = [node0] + (node0_header[grouping:] if get_full_header else []) # deleted + node0_header[:grouping - 1]  after first part
                node1_header = line[1].split('|')
                node1 = node1_header[grouping - 1]
                node1_ext = [node1] + (node1_header[grouping:] if get_full_header else []) # deleted + node1_header[:grouping - 1] 
                nodes.add(tuple(node0_ext))
                nodes.add(tuple(node1_ext))

                # consider once
                if consider_once:
                    edge_once = (node0, node1, line[0])
                    try:
                        node_once[edge_once]
                        move_on = True
                    except KeyError:
                        node_once[edge_once] = None
                        move_on = False
                    if move_on: 
                        line = inf.readline().strip().split('\t')
                        continue


                try:
                    node_groups[node0][line[0]] += 1
                    node_groups[(node0, node1)]['total'] += 1
                    node_groups[node0]['unique'].add(line[0])
                except KeyError:
                    try:
                        node_groups[node0][line[0]] = 1
                        node_groups[node0]['unique'].add(line[0])
                    except KeyError:
                        node_groups[node0] = {line[0]: 1} 
                        node_groups[node0]['unique'] = set()
                        node_groups[node0]['unique'].add(line[0])

                    try:
                        node_groups[(node0, node1)]['total'] += 1
                    except KeyError:
                        node_groups[(node0, node1)] = {'total': 1}

                try:
                    node_groups[node1][line[1]] += 1
                    node_groups[(node1, node0)]['total'] += 1
                    node_groups[node1]['unique'].add(line[1])
                except KeyError:
                    try:
                        node_groups[node1][line[1]] = 1
                        node_groups[node1]['unique'].add(line[1])
                    except KeyError:
                        node_groups[node1] = {line[1]: 1}
                        node_groups[node1]['unique'] = set()
                        node_groups[node1]['unique'].add(line[1])

                    try:
                        node_groups[(node1, node0)]['total'] += 1
                    except KeyError:
                        node_groups[(node1, node0)] = {'total': 1}

            if identical:
                if node0 == node1:
                    line = inf.readline().strip().split('\t')
                    continue

            if 'Evalue=float' in fields:
                e = float(line[fields.index('Evalue=float')])
            else:
                e = 0.0

            if 'QueryCoverage=integer' in fields:
                q = float(line[fields.index('QueryCoverage=integer')])
            else:
                q = 0

            if 'PercIdentity=float' in fields:
                p = float(line[fields.index('PercIdentity=float')])
            else:
                p = 0.0

            if e > evalue or q < qcovhsp or p < pid:
                line = inf.readline().strip().split('\t')
                continue

            # multiple qcovhsp * pident
            if add_pq:
                pq = p * q / 100
                line += [str(pq)]

            if full_cure:
                fcure[(line[0], line[1])] = ''
            # if grouping == 0:
            #   edges[(node0, node1)] = '\t'.join([node0_header[0], node1_header[0]] + line[2:]) + '\n'
            # else:

            edge = (node0, node1)
            try:
                # edges[(node0, node1)].append('\t'.join([node0, node1] + line[2:]) + '\n')
                edges[edge].append('\t'.join(line))
            except KeyError:
                # edges[(node0, node1)] = ['\t'.join([node0, node1] + line[2:]) + '\n'] # jesli n0 i n1 te same to zawsze bedzie brany pod uwage tylko ostatnie !!!!
                edges[edge] = ['\t'.join(line)]
            line = inf.readline().strip().split('\t')

    print('  Read  100%% (%i / %i) BLAST hits.  ' % (cnt, file_length))
    print(f"  Kept {len(fcure)}({len(fcure)/file_length * 100:0.2f}%) BLAST hits based on provided thresholds.")

    node_once = None

    return nodes, edges, fcure, node_groups

def cureEdges(edges, grouping, full_cure, fcure):

    # edges     - {(node1, node2): [attributes, ...]}
    # fcure     - {(node1, node2): '', ...}

    print('Checking nodes reciprocation...')
    newEdges = {}
    ur_cnt = 0
    fc_cnt = 0
    fc_len = len(fcure)
    edge_cnt = 0
    edge_len = 0 # len(edges)
    for np in edges:
        edge_len += len(edges[np])

    if full_cure:
        # first check all edges for their reciprocation based on full sequence headers
        new_fc_edges = {}
        for edge in fcure: # takes whole headers as node names
            fc_cnt += 1
            if fc_cnt % 10000 == 0 :
                print('  Read  %0.2f%% (%i / %i) edges during full sequence header curation.' % (fc_cnt/fc_len * 100, fc_cnt, fc_len), end='\r', flush=True)
            try:
                x = fcure[(edge[1], edge[0])]
                # node0 = edge[0].split('|')[grouping - 1]
                # node1 = edge[1].split('|')[grouping - 1]
                new_fc_edges[(edge[0], edge[1])] = ''
            except KeyError:
                continue
        print('  Read  %0.2f%% (%i / %i) edges during full sequence header curation.' % (fc_cnt/fc_len * 100, fc_cnt, fc_len))

        # then check all edges between two nodes: either part of sequence header or full sequence header
        for nodes_pair in edges:
            for edge in edges[nodes_pair]:
                edge = edge.split('\t')
                edge_cnt += 1
                if edge_cnt % 10000 == 0 :
                    print('  Processed %0.2f%% (%i / %i) edges for grouped nodes.' % (edge_cnt/edge_len * 100, edge_cnt, edge_len), end='\r', flush=True)
                # check whether the edge has reciprocation based on full sequence header
                try:
                    new_fc_edges[edge[1], edge[0]]
                    try:
                        newEdges[nodes_pair].append('\t'.join(list(nodes_pair) + edge[2:]))
                    except KeyError:
                        newEdges[nodes_pair] = ['\t'.join(list(nodes_pair) + edge[2:])]
                except KeyError:
                    #print('  %s and %s do not reciprocate. (%i/%i)' % (edge[1], edge[0], edge_cnt, edge_len))
                    ur_cnt += 1

    else:
        for nodes_pair in edges:
            for edge in edges[nodes_pair]:
                edge = edge.split('\t')
                edge_cnt += 1
                if edge_cnt % 10000 == 0 :
                    print('  Processed %0.2f%% (%i / %i) edges.' % (edge_cnt/edge_len * 100, edge_cnt, edge_len), end='\r', flush=True)
                try:
                    new_fc_edges[edge[1], edge[0]]
                    try:
                        newEdges[nodes_pair].append('\t'.join(list(nodes_pair) + edges[2:]))
                    except KeyError:
                        newEdges[nodes_pair] = ['\t'.join(list(nodes_pair) + edges[2:])]
                except KeyError:
                    print('  %s and %s do not reciprocate. (%i/%i)' % (edge[1], edge[0], edge_cnt, edge_len))
                    ur_cnt += 1
        # for edge in edges:
        #   edge_cnt += 1
        #   if edge_cnt % 10000 == 0 :
        #       print('  Processed %0.2f%% (%i / %i) edges.' % (edge_cnt/edge_len * 100, edge_cnt, edge_len), end='\r', flush=True)
        #   try:
        #       new_fc_edges[edge[1], edge[0]]
        #       newEdges[edge] = edges[edge]
        #   except KeyError:
        #       print('  %s and %s do not reciprocate. (%i/%i)' % (edge[1], edge[0], edge_cnt, edge_len))
        #       ur_cnt += 1

    print('  Processed %0.2f%% (%i / %i) edges%s.' % (edge_cnt/edge_len * 100, edge_cnt, edge_len, '' if not grouping else ' for grouped nodes'))
    print('  Total number of unreciprocated nodes: %i' % ur_cnt)

    return newEdges


def writeEdges(edges, fields, outfile, merge_edges, addpq):

    edge_len = 0 # len(edges)
    if merge_edges:
        for nodes_pair in edges:
            edge_len += 1 
            first_edge = edges[nodes_pair][0]
            if addpq:
                merged_weigths = [str(len(edges[nodes_pair])), 0.0]
                for e in edges[nodes_pair]:
                    weights = e.split('\t')[-2:]
                    merged_weigths[1] += float(weights[1])
            else:
                merged_weigths = [str(len(edges[nodes_pair]))]
                # merged_weigths = [0.0]
                # for e in edges[nodes_pair]:
                #     weights = e.split('\t')[:-1]
                #     merged_weigths[0] += float(weights[0])

            edges[nodes_pair] = ['\t'.join(first_edge.split('\t')[:2] + ['\t'.join([str(x) for x in merged_weigths])])]
    else:
        for nodes_pair in edges:
            edge_len += len(edges[nodes_pair])

    print('Writing %i edges to file...' % edge_len)
    with open(outfile + '.edges', 'w') as edgesout:
        if merge_edges:
            if not addpq:
                edgesout.write('#' + '\t'.join(fields[:2]) + '\tNumberOfSharedProts=float\n')
            else:
                edgesout.write('#' + '\t'.join(fields[:2]) + '\tNumberOfSharedProts=float\tMergedPQ=float\n')
        else:
            edgesout.write('#' + '\t'.join(fields) + '\n')
        for nodes_pair in edges:
            for edge in edges[nodes_pair]:
                edgesout.write(edge + '\n')
        # for edge in edges:
        #   edgesout.write(edges[edge])


def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(_nsre, s)]   


def writeNodes(nodes, att_names, grouping, extra_att, outfile):

    # read attributes if provided
    print('Writing %i nodes to file...' % len(nodes))
    new_nodes = []
    extra_att_names = []
    if extra_att:
        extra_nodes_atts = {}
        with open(extra_att) as inf:
            extra_att_names = inf.readline().strip().split('\t')[1:] # assumes the first line has attributes names
            for line in inf:
                line = line.strip().split('\t')
                extra_nodes_atts[line[0]] = line[1:]

        if grouping == 0:
            new_nodes = ['\t'.join([node] + node.split('|') + extra_nodes_atts[node.split('|')[0]]) for node in nodes]
        else:
            new_nodes = ['\t'.join(list(node) + extra_nodes_atts[node[0]]) for node in nodes]
    else:
        if grouping == 0:
            new_nodes = ['\t'.join([node] + node.split('|')) for node in nodes]
        else:
            new_nodes = ['\t'.join(list(node)) for node in nodes]


    with open(outfile + '.nodes', 'w') as nodesout:
        nodesout.write('#%s\n' % ('\t'.join(att_names + (extra_att_names if extra_att else []))))
        for node in sorted(new_nodes, key=natural_sort_key):
                nodesout.write('%s\n' % (node))


def get_node_groups(infile, grouping, node_groups):

    gcnt        = 0
    cnt         = 0

    with open(infile) as inf:
        for line in inf:
            if line.startswith('>'):
                gcnt += 1
                header = line[1:].strip()
                node0_header = header.split('|')
                node0 = node0_header[grouping - 1]

                # check wether the node was present in BLAST output
                try:
                    node_groups[node0]
                except KeyError:
                    node_groups[node0] = {header: 1}
                    node_groups[node0]['total'] = 1
                    node_groups[node0]['unique'] = set()
                    node_groups[node0]['unique'].add(header)
                    cnt += 1
                    continue

                # check wether the protein was present in BLAST output
                try:
                    node_groups[node0][header]
                    # node_groups[node0]['total'] += node_groups[node0][header]
                    cnt += 1
                    continue
                except KeyError:
                    node_groups[node0][header] = 1
                    node_groups[node0]['unique'].add(header)
                    # node_groups[node0]['total'] += 1
                    cnt += 1

    print('  Read %i groups with %i proteins in total.' % (cnt, gcnt))

    # print('GCA_900106865.1_genomic--NZ_FNUE01000001.1--pp_8: ',
    #     len(node_groups['GCA_900106865.1_genomic--NZ_FNUE01000001.1--pp_8']), 
    #     node_groups['GCA_900106865.1_genomic--NZ_FNUE01000001.1--pp_8']['total'])
    # print('GCA_001280865.1_genomic--NZ_LGBR01000001.1--pp_7: ', len(node_groups['GCA_001280865.1_genomic--NZ_LGBR01000001.1--pp_7']), node_groups['GCA_001280865.1_genomic--NZ_LGBR01000001.1--pp_7']['total'])

    return node_groups


def normalize_pq(edges, node_groups, mode):

    m = len(edges)
    print(f"  Normalizing PQ value in '{mode}' mode.")
    cnt = 0
    for pair in edges:
        # print(pair)
        n1, n2 = pair
        cnt += 1
        if cnt % 10000 == 0: 
            print('  Processed  %0.2f%% (%i / %i) edges.' % (cnt/m * 100, cnt, m), end='\r', flush=True)
        new_edges = []
        for edge in edges[pair]:
            edge            = edge.split('\t')
            edge_weigth     = float(edge[-1]) # PQ value
            node_group_size = len(node_groups[n1])
            try:
                # normalized by the number of elements in group (0, >100]
                # below calculation just normalizes the total score by the number of encoded proteins which may result in the total score higher than 100 when one protein is shared with more than one from the other node
                if mode == 'elements':
                    pq = edge_weigth / node_group_size 
                    # pq /= 100
                    # print(edge_weigth, node_group_size, pq)

                # normalized by the number of similarities (0, 100]
                # below calculation skews the PQ and penalizes the weigth when one protein is homologous to more than one from the other node, eg. 18/20 shared with 100 should result in total PQ of 90 but having multiple weak similarities increases the divider and influences the total score
                elif mode == 'similarities':
                    total = node_groups[(n1, n2)]['total']
                    pq = edge_weigth / total * 2 # / node_group_size
                    # print(edge_weigth, total, pq)

                # normalized by the number of similarities weighted by the percentage of shared
                elif mode == 'weighted_similarities':
                    shared   = node_groups[(n1, n2)]['total'] / 2
                    unique   = len(node_groups[n1]['unique'])
                    fraction = unique / node_group_size

                    pq = edge_weigth / shared * fraction
                    # print(edge_weigth, shared, unique, node_group_size, fraction, pq)

                # normalize by the number of elements in group and weight it by the percentage of shared
                elif mode == 'weighted_elements':
                    unique   = len(node_groups[n1]['unique'])
                    fraction = unique / node_group_size

                    pq = edge_weigth / node_group_size * fraction
                    # print(edge_weigth, node_group_size, unique, pq)

                elif mode == "both_groups_elements":
                    pq = edge_weigth / (len(node_groups[n1]) + len(node_groups[n2]))

            except KeyError:
                print('WARNING: %s or %s are not present in node_groups.' % (n1, n2))
                pq = 5
            
            edge[-1] = str(pq)
            new_edges.append('\t'.join(edge))
        edges[pair] = new_edges

    print('  Processed  %0.2f%% (%i / %i) edges.' % (cnt/m * 100, cnt, m))
    return edges


def main():
    args = ArgumentParser('Program parsing BLAST output file. To use full headers provide -g 0 and attributes for Node,att_X,att_Y if two parts would result in splitting the sequence header by \'|\'')
    
    args.add_argument('-i', '--infile',
                      type = str,
                      help = 'Path to BLAST outfile.',
                      required = True)
    
    args.add_argument('-o', '--outfile',
                      type = str,
                      help = 'Path to output file - leave without extension.')
    
    args.add_argument('--extra_att',
                      type = str,
                      help = 'Path to input file with new attributes in the first line and the first column matching nodes names. Eveything should be a tab-delimited table.')

    args.add_argument('-a', '--att',
                      type = str,
                      help = 'Names for nodes attributes if grouping selected. [Default: %(default)s]',
                      default = 'Node')

    args.add_argument('-f', '--fields',
                      type = str,
                      help = 'Comma-separated columns names for blast outfmt 6. At least 2 names are required. Ignored '
                           'if --preset chosen.')
    
    args.add_argument('-g', '--grouping',
                      type = int,
                      help = 'Sequence header column after splitting on | sign. If not set, whole sequence headers are used for nodes and edges preparation. If 0 set, the whole sequence header will be used for Node name and the rest will be splitted into separate columns. Add one additional parameter name then. [Default: not set]')

    args.add_argument('-r', '--preset',
                      type = str,
                      help = 'Presets available: '
                           'std == --outfmt 6, '
                           'cqovhsp == --outfmt \'6 std qcovhsp\''
                           'short == --outfmt \'6 qseqid sseqid pident evalue bitscore qcovhsp\'')
    
    args.add_argument('-e', '--evalue',
                      type = float,
                      help = 'Maximum e-value. [Default: %(default)s]',
                      default = 10.0)
    
    args.add_argument('-q', '--qcovhsp',
                      type = float,
                      help = 'Minimum query coverage versus subject. [Default: %(default)s]',
                      default = 0)
    
    args.add_argument('-p', '--pid',
                      type = float,
                      help = 'Minimum %% identity. [Default: %(default)s]',
                      default = 0.0)

    args.add_argument('-ce', '--cure_edges',
                      action = 'store_true',
                      help = 'Cure edges - leave only reciprocated BLAST hits (A -> B and B -> A). Default = None')

    args.add_argument('-fhec', '--full_header_edge_curing',
                      action = 'store_true',
                      help = 'If set, despite grouping parameter, full sequence headers will be used for checking similarity reciprocation. [Default: not set - works for non-grouping run]')

    args.add_argument('--addpq',
                      action = 'store_true',
                      help = 'If set, qcovhsp * pident column for edge will be calculated and added.')

    args.add_argument('--normpq',
                      type = str,
                      help = 'Path to multifasta file used for BLAST search. If set, the PQ value will be divided by the sum of records identified for both nodes it indicated the PQ for. WARNING: works only when --grouping applied.')

    args.add_argument('--norm_mode',
                      type = str,
                      choices = ['elements', 'similarities', 'weighted_similarities', 'weighted_elements', 'both_groups_elements'],
                      help = 'Type of the normalization of the weights: elements - average weight per elements number within a group (0; > 100), similarities - average similarity (0, 100], weighted_similarities - similarities * fraction of shared between groups (0; >100 without --consider_once else 0;100], weighted_elements - elements * fraction of shared (0,100]. [Default: %(default)s]',
                      default = 'both_groups_elements')

    args.add_argument('--consider_once',
                      action = 'store_true',
                      help = 'If set, each protein from one organism (first part of the header) shared with other organisms will be considered only once. This will make sure that weak similarities between multiple proteins will not skew the weight of edges. Then when --normpq applied you can use option weighted_similarities in --norm_mode flag.')

    args.add_argument('--get_full_header',
                      action = 'store_true',
                      help = 'If set, complete header will be added to nodes description during grouping the nodes. Requires the complete set of names for each parameter provided in sequence header. USE ONLY IF the number of sequences and groups after grouping is the same otherwise the nodes number will refer to complete sequence headers. So the numbers of nodes will increase but the number of edges will be the same (not like in grouping 0).')

    args.add_argument('-id', '--identical',
                      action = 'store_true',
                      help = 'Remove if identical query and subject. Default = None')

    args.add_argument('-me', '--merge_edges',
                      action = 'store_true',
                      help = 'If set, edges will be merged and only weight parameter will be written.')
    
    if len(argv[1:]) == 0:
        args.print_usage()
        args.exit()
    
    try:
        args = args.parse_args()
    except:
        args.exit()
    

    if args.preset == 'std':
        args.fields = ['Node', 'Node', 'PercIdentity=float', 'Length', 'Mismatch', 'GapOpen', 'Qstart', 'Qend',
                       'Sstart', 'Send', 'Evalue=float', 'BitScore=float']
    elif args.preset == 'qcovhsp':
        args.fields = ['Node', 'Node', 'PercIdentity=float', 'Length', 'Mismatch', 'GapOpen', 'Qstart', 'Qend',
                       'Sstart', 'Send', 'Evalue=float', 'BitScore=float', 'QueryCoverage=integer']
    elif args.preset == 'short':
        args.fields = ['Node', 'Node', 'PercIdentity=float', 'Evalue=float', 'BitScore=float', 'QueryCoverage=integer']
    else:
        args.fields = args.fields.split(',')
        if len(args.fields) < 2:
            print('\nThere are less than two columns names provided. Quiting...\n')
            exit()

    if args.addpq:
        args.fields += ['PQ=float']

    nodes, edges, fcure, node_groups = reformat(args.infile, args.fields, args.grouping, args.evalue, args.qcovhsp, args.pid, args.identical, args.full_header_edge_curing, args.get_full_header, args.addpq, args.consider_once)

    if args.cure_edges:
        edges = cureEdges(edges, args.grouping, args.full_header_edge_curing, fcure)

    if args.normpq:
        node_groups = get_node_groups(args.normpq, args.grouping, node_groups)
        edges = normalize_pq(edges, node_groups, args.norm_mode)

    writeEdges(edges, args.fields, args.outfile, args.merge_edges, args.addpq)

    att_names = [args.att] if ',' not in args.att else args.att.split(',')
    writeNodes(nodes, att_names, args.grouping, args.extra_att, args.outfile)

    print('Done. Use your output files with anNe.\n')
            
if __name__ == '__main__':
    main()
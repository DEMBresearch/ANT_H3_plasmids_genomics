#!/usr/bin/env python3
__author__ = 'Przemek Decewicz'
__date__ = '17.06.2021'

from argparse import ArgumentParser
from math import ceil, sqrt
from sys import argv, setrecursionlimit
from Graph import Graph

def read_groups_order(groupsOrderFile):

    groupsOrder = dict()
    
    with open(groupsOrderFile) as inf:
        line = inf.readline().strip()
        key = ''

        while line != '':
            if '>' in line[0]:
                key = line[1:]
            else:
                try:
                    groupsOrder[key].append(line)
                except KeyError:
                    groupsOrder[key] = [line]
            line = inf.readline().strip()

    return groupsOrder


def read_orginized_groups_order(groupsOrderFile, G, gAttrs):

    groupsOrder = dict()

    if not gAttrs:
        with open(groupsOrderFile) as inf:
            line = inf.readline().strip().split('\t')
            
            if ';' not in line[1]:
                # only in 1 line
                order = line[1].split(',')
            else:
                order = [tuple(x.split(',')) if ',' in x else (x, ) for x in line[1].split(';')]

            groupsOrder[line[0]] = order
    else:
        for gAtt in gAttrs:
            att_values = sorted(G.groups[gAtt].keys())
            columns = int(ceil(float(sqrt(len(att_values)))))
            indexes = list(range(columns, len(att_values), columns))
            groupsOrder[gAtt] = [tuple(att_values[i - columns:i]) for i in indexes]
            if indexes[-1] < len(att_values):   groupsOrder[gAtt].append(tuple(att_values[indexes[-1]:]))

    return groupsOrder


def main():
    args = ArgumentParser()
    
    args.add_argument('-n', '--nodes',
                      type=str,
                      help='Path to nodes file.',
                      required=True)

    args.add_argument('-e', '--edges',
                      type=str,
                      help='Path to edges file.',
                      required=True)
    
    args.add_argument('-o', '--outfile',
                      type=str,
                      help='Path to output file.',
                      required=True)

    args.add_argument('-t', '--type',
                      type=str,
                      help='Type of layout: simple (all nodes at x=0, y=0), newOrbit, circular or orgAtt (requires --groupsorder !!)  [Default: '
                           'simple]',
                      default='simple',
                      required=False)

    args.add_argument('-w', '--min_weight',
                      type = float,
                      help = 'Filter edges by minimum weight. [Default: 0]',
                      default = 0.0)

    args.add_argument('-d', '--divide_weight',
                      type = float,
                      help = 'If set, the weight will be divided by certain number.')

    args.add_argument('--filter_nodes',
                      type = str,
                      help = 'Path to input file with node names to retain.')

    args.add_argument('--keep_neighbours',
                      action = 'store_true',
                      help = 'If set, the nodes provided with --filter_nodes will be retained with their neighbours.')

    args.add_argument('-a', '--attcore',
                      type=str,
                      help='Attribute core analysis. If provided core analysis summary will be performed. Requires a '
                           'name of an attribute to use during the analysis.'
                           '[Default: None]',
                      default = None,
                      required=False)

    args.add_argument('-cr', '--colrow',
                      type=str,
                      help='Number of columns and rows (comma-separated numbers) for circular layout. These, '
                           'multiplied should give higher or equal number of tiles needed to plot all groups. '
                           'Otherwise it will be calculated automatically.'
                           '[Default: None]',
                      default = None,
                      required=False)

    args.add_argument('-g', '--groups',
                      type=str,
                      help='Comma-separated Node attributes for grouping Nodes when performing circular layout. If '
                           'not set, all attributes will be used. [Default: None]',
                      default=None,
                      required=False)

    args.add_argument('-go', '--groupsorder',
                      type=str,
                      help='If provided, order for given groups will be applied. Requires format: ">Attributename" '
                           'and possible attribute values each in new line. [Default: alphabetical]',
                      default=None,
                      required=False)

    args.add_argument('-ef', '--elipse',
                      type=float,
                      help='Elipse factor - a multiplier float number between 0.25 to 1 which flatens the circle '
                           'from the top and bottom to make an elipse. [Default: 1]',
                      default = 1.0,
                      required=False)

    args.add_argument('-s', '--nodesize',
                      type=int,
                      help='Node size. [Default: 10]',
                      default = 10,
                      required=False)

    args.add_argument('-od', '--orbdist',
                      type=float,
                      help='Stretching of subgraph orbit. Float number between 0.1 to 3 where the higher number the '
                           'smaller the distance between orbits [Default: 0.6]',
                      default = 0.6,
                      required=False)

    args.add_argument('-id', '--interorbdist',
                      type=float,
                      help='Distance between next orbits. Float number between 1.0 to 2.0 where the higher number the '
                           'smaller the distance between orbits [Default: 1.2]',
                      default = 1.2,
                      required=False)

    args.add_argument('-og', '--orbdens',
                      type=float,
                      help='Density on orbits. Float number between 2.0 to 6.0 where the higher number the '
                           'smaller number of subgraphs on orbits [Default: 4.0]',
                      default = 4.0,
                      required=False)

    args.add_argument('-f', '--fasta',
                      type=str,
                      help='Path to FASTA-formatted sequences file corresponding to nodes. Sequence headers have to '
                           'match those in .nodes file. [Default: None]',
                      required=False)

    args.add_argument('-ss', '--skip_summary',
                      action = 'store_true',
                      help = 'If set, finding subgraphs and writting the summary of those will be skipped. Useful for very big graphs.')

    if len(argv[1:]) == 0:
        args.print_help()
        args.exit()
    
    try:
        args = args.parse_args()
    except:
        args.exit()


    if args.elipse < 0.25 or args.elipse > 1.0:
        print('\nWrong elipse factor argument value provided. The number must be in range 0.25 to 1. Quiting...\n')
        exit(1)

    setrecursionlimit(10000)
    G = Graph(args.type, args.nodesize)

    # read nodes
    print('Reading NODES file...')
    G.readNodes(nodeFile=args.nodes)

    # read edges
    print('Reading EDGES file...')
    G.readEdges(edgeFile=args.edges, min_weight=args.min_weight, divide_weight=args.divide_weight)

    if not args.skip_summary:
        # find subgraphs
        print('Finding subgraphs...')
        G.findSubGraphs()

        # write subGraphs summary
        print('Writing subgraphs summary...')
        G.writeSubGraphs(outfile=args.outfile, addseqs=args.fasta)

    # write attribute based core and singleton analysis
    if args.attcore:
        print('Performing attribute core analysis for %s attribute(s).' % args.attcore)
        G.attributeBasedCore(outfile=args.outfile, attributes=args.attcore, addseqs=args.fasta)

    # run layout
    print('Working on %s layout...' % args.type)
    if args.type == 'newOrbit':
        G.newOrbit(args.orbdist, args.orbdens, args.interorbdist)
        # export to GraphML
        G.export2graphml(outfile=args.outfile)
    elif args.type == 'circular':

        if args.colrow:
            colrow = args.colrow.split(',')
            columns = float(colrow[0])
            rows = float(colrow[0])
        else:
            columns = None
            rows = None

        gAttrs = args.groups.split(',') if args.groups else G.getNodeAtts()
        if args.groupsorder: groupsOrder = read_groups_order(args.groupsorder)
        else: groupsOrder = []
        G.groupNodes(gAttrs=gAttrs)
        for group in gAttrs:
            groupOrder = groupsOrder[group] if group in groupsOrder else None
            G.circleLayout(columns=columns, rows=rows, group=group, grouporder=groupOrder, elipse_factor=args.elipse)
            G.export2graphml(outfile=args.outfile, group=group)

    elif args.type == 'orgAtt':
        if not args.groupsorder and not args.groups:
            print('\n\t! Group order file or groups names not provided !\n')
            exit(1)

        gAttrs = args.groups.split(',') if args.groups else G.getNodeAtts()
        G.groupNodes(gAttrs=gAttrs)
        groupsOrder = read_orginized_groups_order(args.groupsorder, G, None if args.groupsorder else gAttrs)

        for group in gAttrs:
            if group in groupsOrder:
                groupOrder = groupsOrder[group]
            else:
                continue
            G.organizedAttributeLayout(group=group, grouporder=groupOrder, orbit_distance=args.orbdist,
                                       inter_orbital_dist=args.interorbdist)
            G.export2graphml(outfile=args.outfile, group=group)

    else:
        G.simple()
        if args.filter_nodes:
            G.filter_nodes(args.filter_nodes, args.keep_neighbours)
        G.export2graphml(outfile=args.outfile)


if __name__ == '__main__':
    main()
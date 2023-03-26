from Node import Node
from Edge import Edge
from math import sin, cos, radians, floor, ceil, pi, sqrt, log10
import numpy as np
from os import path, makedirs

class Graph:
    """This is the graph class."""

    def __init__(self, type, nodesize):
        self.edges        = []     # [Edge class objects]
        self.nodes        = dict() # {nX: Node class object}
        self.nodesConnect = dict() # {nodeName: nX}
        self.directed     = False  # (un)directed
        self.type         = type   # layout type
        self.subGraphs    = {}     # dictionary of subgraphs Size: SET of tuples
        self.groups       = {}     # dictionary of groups based on attributes
        self.nodeSize     = nodesize # the initial size of the nodes
        self.nodeAttNames = ['x=float', 'y=float', 'degree=float'] # nodes attributes
        self.edgeAttNames = []     # edges attributes

    def getNodeAtts(self):
        return self.nodeAttNames[2:]

    def readNodes(self, nodeFile):
        cnt = 0

        with open(nodeFile) as inf:
            line = inf.readline().strip().split('\t')
            # Check whether columns names were provided
            if len(line) > 1:
                if line[0].startswith('#'):
                    line = [line[0][1:]] + line[1:]
                    attNames = [att for att in line]
                    self.nodes
                else:
                    attNames = ['Att_%i' % i for i in range(1, len(line))]
                    cnt += 1
                    self.nodes['n%05d' % cnt] = Node('n%05d' % cnt, line[0], line, attNames)
                    self.nodesConnect[line[0]] = 'n%05d' % cnt
            else:
                #TODO messes when only 1 attribute is provided, why????
                print('Warning!!!! Only 1 column in NODES file provided. It messes up with edges number...')
                attNames = line[1:]
                cnt += 1
                self.nodes['n%05d' % cnt] = Node('n%05d' % cnt, line[0], line[0], attNames)
                self.nodesConnect[line[0]] = 'n%05d' % cnt

            line = inf.readline().split('\t')
            line[-1] = line[-1].strip()
            while line != ['']:
                cnt += 1
                self.nodes['n%05d' % cnt] = Node('n%05d' % cnt, line[0], line, attNames)
                self.nodesConnect[line[0]] = 'n%05d' % cnt
                line = inf.readline().split('\t')
                line[-1] = line[-1].strip()

        print('  Read %i nodes.' % len(self.nodes))
        self.nodeAttNames.extend(attNames)

    def readEdges(self, edgeFile, min_weight, divide_weight):
        cnt = 0

        with open(edgeFile) as inf:
            line = inf.readline().strip().split('\t')

            # Check whether columns names were provided
            if len(line) > 2:
                if '#' in line[0]:
                    attNames = [att for att in line[2:]]
                else:
                    attNames = ['Att_%i' % i for i in range(2, len(line))]
                    cnt += 1
                    # create edge
                    source = self.nodes[self.nodesConnect[line[0]]].getID()
                    target = self.nodes[self.nodesConnect[line[1]]].getID()
                    self.edges.append(
                        Edge(cnt,
                             source,
                             target,
                             line[2:],
                             attNames))

                    # update nodes' informations
                    self.nodes[source].increaseDegree()
                    self.nodes[source].addNeighbour(target)
                    self.nodes[target].increaseDegree()
                    self.nodes[target].addNeighbour(source)

            line = inf.readline().strip().split('\t')
            while line != ['']:
                # filter edges with too small weight
                if float(line[-1]) < min_weight:
                    line = inf.readline().strip().split('\t')
                    continue
                if divide_weight:
                    line[-1] = str(float(line[-1]) / divide_weight)
                cnt += 1
                source = self.nodes[self.nodesConnect[line[0]]].getID()
                target = self.nodes[self.nodesConnect[line[1]]].getID()
                self.edges.append(
                    Edge(cnt,
                         source,
                         target,
                         line[2:],
                         attNames))

                # update nodes' informations
                self.nodes[source].increaseDegree()
                self.nodes[source].addNeighbour(target)
                self.nodes[target].increaseDegree()
                self.nodes[target].addNeighbour(source)

                line = inf.readline().strip().split('\t')

        print('  Read %i edges.' % len(self.edges))
        self.edgeAttNames = attNames

    def checkNeighbours(self, node, subGraph):

        #TODO change to iterative version rather than recursive or when catching the RecursiveError change mode
        #nsize = len(subGraph)
        subGraph.add(node.getID())
        toCheck = node.getNeighbours() - subGraph
        #if nsize == len(subGraph): # means nothing new in subGraph
        #   return subGraph

        for n in toCheck:
            if n in subGraph: continue
            subGraph.add(self.nodes[n].getID())
            subGraph = self.checkNeighbours(self.nodes[n], subGraph)

        return subGraph

    def findSubGraphs(self, ):

        for node in self.nodes:

            # check to what nodes its connected with
            subgraph = self.checkNeighbours(self.nodes[node], set())  # set of nodes that create subgraph

            try:
                self.subGraphs[len(subgraph)].add(tuple(sorted(list(subgraph))))
            except KeyError:
                self.subGraphs[len(subgraph)] = {tuple(sorted(list(subgraph)))}

        # each orbit represents different subgraph_size
        if self.type == 'orbital':
            # leave it as it is
            pass

        # double orbital layout: orbital layouts of orbital layouts of each subgraph
        elif self.type == 'fullorbital':

            if 1 in self.subGraphs:
                tmp = {0: [], 1: self.subGraphs.pop(1)}
            else:
                tmp = {0: []}

            cnt = 0
            for subgraph_size in sorted(self.subGraphs, reverse=True):
                cnt += subgraph_size
                tmp[0].extend(self.subGraphs.pop(subgraph_size))

            tmp[cnt] = tmp.pop(0)

            self.subGraphs = tmp


        #for key in range(2, 4):
        #   print(key, self.subgraphs[key])

    def groupNodes(self, gAttrs):
        """The function groups nodes into subgraphs based on the attributes they have assigned. It changes the
        self.groups attribute of Graph class."""

        # check whether selected attribute is actually in Nodes attributes
        for gAttr in gAttrs:
            if gAttr not in self.nodeAttNames:
                print('Wrong attribute "%s" provided. Check Nodes file and provided attributes.' % gAttr)
                raise KeyError

        for node in self.nodes:
            for gAttr in gAttrs:
                attVal = self.nodes[node].getAttVal(gAttr)

                try:
                    self.groups[gAttr][attVal].append(node)
                except KeyError:
                    try:
                        self.groups[gAttr][attVal] = [node]
                    except KeyError:
                        self.groups[gAttr] = {attVal: [node]}

    def filter_nodes(self, nodes_file, keep_neighbours):
        """Modifies nodes and edges to keep only those indicated in nodes_file with or without neighbours"""

        print('Filtering nodes...')

        keep_nodes = {}
        keep_edges = []

        with open(nodes_file) as inf:
            for line in inf:
                if line.startswith('#'): continue
                keep_nodes[line.strip().split('\t')[0]] = None
        print('  Read %i nodes from file.' % len(keep_nodes))

        # extend the set of nodes to keep with neighbours if set

        if keep_neighbours:
            neighbours = {}
            print('Keeping neighbours for selected nodes...')
            for edge in self.edges:
                source_name = self.nodes[edge.getA()].name
                target_name = self.nodes[edge.getB()].name
                try:
                    keep_nodes[source_name]
                    neighbours[target_name] = None
                    continue
                except KeyError:
                    try:
                        keep_nodes[target_name]
                        neighbours[source_name] = None
                    except KeyError:
                        pass
            keep_nodes.update(neighbours)
            del neighbours
            print('  Keeping %i nodes in total.' % len(keep_nodes))

        # go through edges and keep only those for retained nodes
        print('Filtering edges...')
        for edge in self.edges:
            source_name = self.nodes[edge.getA()].name
            target_name = self.nodes[edge.getB()].name
            try:
                keep_nodes[source_name]
                keep_nodes[target_name]
                keep_edges.append(edge)
            except KeyError:
                pass
        print('  Keeping %i edges.' % len(keep_nodes))

        self.edges = keep_edges
        del keep_edges


    #########################################################################
    ##############################             ##############################
    ##############################   LAYOUTS   ##############################
    ##############################             ##############################
    #########################################################################

    def simple(self, ):

        return

    ###### newOrbit layout

    def newOrbit(self, orbit_distance, orbit_density, inter_orbital_dist):
        """This newOrbit layout will provide a universe of subgraphs lied regularly on orbits with a number of
        subgraphs on each orbit based on its radius and sizes of sorted by size subgraphs."""
        #TODO cure distances between orbit elements and between orbits
        # needs:
        # - self.subGraphs
        # - self.nodeSize
        # default orbit_distance = 0.6  < 0.1 - 3.0 >
        # default orbit_density = 4     < 2.0 - 6.0 >
        orbit_density = 2 * orbit_density

        def newOrbitStart(subgraphs, global_radius, radiuses):
            """Crucial function which prepares each orbit and subgraphs which will be put on it."""

            orbit_len = 2 * pi * global_radius
            subgraphs_diameters_sum = 0.0
            number_of_subgraphs = 0

            # the surver which sums the subgraphs diameters and checks whether the orbit isn't already full
            for radius in radiuses:
                subgraph_space = orbit_density * radius * 2
                if subgraphs_diameters_sum + subgraph_space < orbit_len:
                    subgraphs_diameters_sum += subgraph_space
                    #print(subgraphs_diameters_sum, subgraph_space, orbit_len)
                    number_of_subgraphs += 1
                else:
                    break

            number_of_subgraphs = 1 if global_radius == 0.0 else number_of_subgraphs
            # print(number_of_subgraphs)
            # ceil(global_radius/start_radius/orbit_density)
            to_lay = {}
            laying_order = []
            covered_degrees = radians(0)
            #print(list(range(subgraphs)))
            #print(radiuses)

            # count the percentage share on the orbit of each subgraph that will be layed
            # additional space between subgraphs to share
            dist_space = 0.0 if number_of_subgraphs == 1 else \
                radians(360) * ((1 - subgraphs_diameters_sum / orbit_len) / 2 / number_of_subgraphs)
            degrees_shares = [orbit_density * 2 * radius / subgraphs_diameters_sum
                              for radius in radiuses[:number_of_subgraphs]] if number_of_subgraphs > 1 else [0.0]
            # print(subgraphs_diameters_sum, orbit_len)
            # print('Radians 360 = ', radians(360))
            # print('Radians 360 * share = ', radians(360) * degrees_shares[0])
            # print('Dist space = ', dist_space)
            #print(degrees_shares, sum(degrees_shares))

            for i, subgraph, radius, degree_share in \
                    zip(list(range(number_of_subgraphs)),
                        subgraphs[:number_of_subgraphs],
                        radiuses[:number_of_subgraphs],
                        degrees_shares):

                # wtf is that?
                if len(subgraph) <= 0.65 * number_of_subgraphs:
                    pass #number_of_subgraphs = len(subgraph)

                # prev version
                # subgraph_covered_degrees = (radians(360) / number_of_subgraphs) * orbit_density * radius
                covered_degrees += radians(360) * degree_share * 0.6
                #print(covered_degrees)
                coords = (cos(covered_degrees) * global_radius,
                          sin(covered_degrees) * global_radius)
                covered_degrees += radians(360) * degree_share * 0.4
                #covered_degrees += subgraph_covered_degrees
                to_lay[coords] = [subgraph, radius]
                laying_order.append(coords)


                if i == number_of_subgraphs - 1:
                    break
            #if number_of_subgraphs > 2:
            #   exit()
            return laying_order, to_lay

        subgraphs = []

        # collect subgraphs sizes
        #print(sorted(self.subGraphs))
        for subGraphSize in sorted(self.subGraphs):
            subgraphs.extend(self.subGraphs[subGraphSize])

        # global radius for orbit distance controlling
        global_radius = 0.0
        # R - radius of a certain subgraph
        R = 0.0
        # the sizes of the subgraphs in a decreasing order
        radiuses = []
        for subgraph in subgraphs:
            radiuses.append(((len(subgraph) + 2) * self.nodeSize / orbit_distance + self.nodeSize)/2)

        #print(radiuses)

        while len(subgraphs) > 0: # while there are subgraphs to lay

            laying_order, to_lay = newOrbitStart(subgraphs[::-1], global_radius, radiuses[::-1])
            biggest_to_lay_R = to_lay[laying_order[0]][1]

            for start_coords in laying_order:
                subgraph, R = to_lay[start_coords] # take the subgraph and its radius
                # print('R, subsize: %0.2f, %i' % (R, len(subgraph)))
                #R = ((len(subgraph) + 2) * self.nodeSize / orbit_distance + self.nodeSize)/2

                coords = [(cos(radians(alfa)) * R + start_coords[0],
                           sin(radians(alfa)) * R + start_coords[1])
                           for alfa in np.arange(0, 360, 360 / len(subgraph))]

                for node, coords in zip(subgraph, coords):
                    self.nodes[node].setXY(coords)

            # global radius needs to be increased
            global_radius += (R + inter_orbital_dist * biggest_to_lay_R)
            print('Global radius: ', global_radius)

            tmp = len(to_lay)
            while tmp > 0:
                subgraphs.pop()
                radiuses.pop()
                tmp -= 1


    def organizedAttributeLayout(self, group, grouporder, orbit_distance, inter_orbital_dist):

        """This layout groups nodes based on certain attribute in one line in provided order."""


        print('Preparing layouts for %s attribute.' % group)
        print(group)
        group = self.groups[group]
        #print(group)

        # general radius of each group
        # !!! radius = ((len(subgraph) + 2) * self.nodeSize / orbit_distance + self.nodeSize)/2


        # learn about the dimensions of canvas needed to print the nodes
        height = 0.0
        width = 0.0
        groupsHeights = {}

        for groups in grouporder:
            print(type(groups))
            print('Group:')
            groupHeight = 0.0
            groupWidth = 0.0
            for subGroup in groups:

                print('%6d nodes in \'%s\' subgroup.' % (len(group[subGroup]), subGroup))

                radius = ((len(group[subGroup]) + 2) * self.nodeSize / orbit_distance + self.nodeSize)/2
                if radius > groupHeight: groupHeight = radius
                groupWidth += 2 * radius + inter_orbital_dist * radius

            height += 2 * groupHeight + inter_orbital_dist * groupHeight
            width = groupWidth if groupWidth > width else width
            groupsHeights[groups] = groupHeight

        # determine the starting point
        startPoint = (-width / 2, height / 2)

        # lay the layout
        heightCorrector = 0.0
        for groups in grouporder: # for each row

            widthCorrector = 0.0
            for i, subGroup in enumerate(groups): # for each column in row

                if subGroup not in group:
                    print('Wrong value (%s) for attribute provided in GroupsOrder file. It will not be included in '
                          'results. Moving to the next.' % subGroup)
                    continue

                # correct the starting point of the subgroup
                subGroupSP = (startPoint[0] + widthCorrector, startPoint[1] - groupsHeights[groups] - heightCorrector)

                r = ((len(group[subGroup]) + 2) * self.nodeSize / orbit_distance + self.nodeSize)/2
                R = groupsHeights[groups] - (groupsHeights[groups] - r)

                for node, alfa in zip(sorted(group[subGroup]), np.arange(0, 360, 360 / group[subGroup].__len__())):

                    nR = R - (0.0 if self.nodes[node].degree > 2 else 5 * self.nodeSize)

                    coords = (cos(radians(alfa)) * nR + subGroupSP[0],
                               sin(radians(alfa)) * nR + subGroupSP[1])

                    self.nodes[node].setXY(coords)

                nextr = ((len(group[groups[i + 1]]) + 2) * self.nodeSize / orbit_distance + self.nodeSize)/2 if \
                    i < len(groups) - 1 else 0.0
                widthCorrector += r + inter_orbital_dist * nextr
            heightCorrector += groupsHeights[groups] + inter_orbital_dist * groupsHeights[groups]

        return


    ###### circle/attribute layout

    def circleLayout(self, columns, rows, group, grouporder, elipse_factor):
        """Circle layout is actually an attribute layout which lays the groups of nodes with the same attribute into
        the same circles. """
        #TODO add option to scale the circles and keep the other in line
        #TODO remove regular spaces between circles
        #TODO add option to point how many and which attributes should be in which line
        #TODO make the nodes of lower degree closer/further from the centre - decrease/increase radius by node size

        # group nodes based on certain attribute gattr -> done by groupNodes() function
        # count groups and make a proper starting point for each group
        # for each group create circle / elipse (x/2)^2 + (y/2) ^2 = 1 // x(t) = a cos(t) i y(t) = b cos(t)

        group = self.groups[group]
        groups_number = group.__len__()
        if columns and rows:
            if columns * rows < groups_number:
                print('Provided columns and rows parameters are too small. Changing to proper ones...')
                columns = rows = int(ceil(float(sqrt(groups_number))))
        else:
            columns = rows = int(ceil(float(sqrt(groups_number))))

        biggestSubGroup = 0

        for subGroup in group:
            print('%6d in \'%s\' group.' % (len(group[subGroup]), subGroup))
            if group[subGroup].__len__() > biggestSubGroup: biggestSubGroup = group[subGroup].__len__()

        tileSize = self.nodeSize * biggestSubGroup * 0.5

        y = 0 + (columns / 2 - (0.5 if columns % 2 == 0 else 0)) * tileSize
        x = 0 - (rows / 2 - (0.5 if rows % 2 == 0 else 0)) * tileSize

        print('Biggest group: %i\nTile size: %.2f\nTiles number: %i\nStart point: (%.2f,%.2f)' %
              (biggestSubGroup, tileSize, columns * rows, y, x))

        startPoints = []
        for i in range(columns):
            for j in range(rows):
                startPoints.append((x + j * tileSize, y - i * tileSize))
                #if j < 10:
                #   print((x + j * tileSize, y - i * tileSize))


        for subGroup, tileSP in zip(group if not grouporder else grouporder, startPoints[:len(group) + 1]):
            #nodesPerOrbit = ceil((pi * tileSize) / subGroup.__len__())
            if subGroup not in group:
                print('Wrong value (%s) for attribute provided in GroupsOrder file. It will not be included in '
                      'results. Moving to the next.' % subGroup)
                continue

            coords = [(cos(radians(alfa)) * tileSize / 2.2 + tileSP[0],
                       sin(radians(alfa)) * elipse_factor * tileSize / 2.2 + tileSP[1])
                      for alfa in np.arange(0, 360, 360 / group[subGroup].__len__())]

            #print(group)
            for node, coords in zip(sorted(group[subGroup]), coords):
                self.nodes[node].setXY(coords)
                #print(node[0], self.nodes[node[0]].getXY())



    ########################################################################
    ##############################            ##############################
    ##############################   EXPORT   ##############################
    ##############################            ##############################
    ########################################################################

    ###### Writing graph summary

    def writeSubGraphs(self, outfile, addseqs):
        seqs = {}
        if addseqs:
            with open(addseqs) as inf:
                line = inf.readline()
                while line != '':
                    if '>' in line[0]:
                        header = line[1:].split()[0]
                        seqs[header] = ''
                    else:
                        seqs[header] += line.strip()

                    line = inf.readline()

        #TODO czy zawsze wypisywac orbital summary?
        with open(outfile + '_%s_summary.txt' % self.type, 'w') as outf:
            outf.write('General summary\n'
                       'SubGraphSize\tNumberOfSubGraphs\n')
            for subGraphSize in sorted(self.subGraphs, reverse=True):
                outf.write('%6d\t%6d\n' % (subGraphSize, len(self.subGraphs[subGraphSize])))

            outf.write('\nDetailed summary\n\n')
            for subGraphSize in sorted(self.subGraphs, reverse=True):

                outf.write('----   %i   ----\n\n' % subGraphSize)
                outf.write('NodeID\t%s%s\n' % ('\t'.join(self.nodeAttNames), '' if not addseqs else '\tSequence'))

                separator = len(self.subGraphs[subGraphSize])
                for subGraph in self.subGraphs[subGraphSize]:
                    separator -= 1

                    for node in subGraph:
                        outf.write('%s\t%s%s\n' %
                                   (node,
                                    self.nodes[node].writeSummary(self.nodeAttNames),
                                    '' if not addseqs else '\t' + seqs[self.nodes[node].name]))

                    if separator != 0:
                        outf.write('-------------\n')

                outf.write('\n')

    ###### Export to GraphML format

    def export2graphml(self, outfile, group=None):

        with open(outfile + '_%s%s.graphml' % (self.type, '_' + group if group else ''), 'w') as outf:
            outf.write('<?xml version="1.0" encoding="UTF-8"?>\n'
                       '<graphml xmlns="http://graphml.graphdrawing.org/xmlns">\n')

            # Write header
            keyID = 0
            keyDict = {'node': [], 'edge': []}
            for nodeAtt in self.nodeAttNames:

                if '=' in nodeAtt:
                    nodeAtt = nodeAtt.split('=')
                    attName = nodeAtt[0]
                    attType = nodeAtt[1]
                else:
                    attName = nodeAtt
                    attType = 'string'

                outf.write('<key id="d%i" for="node" attr.name="%s" attr.type="%s"></key>\n' %
                           (keyID, attName, attType))
                keyDict['node'].append(keyID)
                keyID += 1

            for edgeAtt in self.edgeAttNames:
                if '=' in edgeAtt:
                    edgeAtt = edgeAtt.split('=')
                    attName = edgeAtt[0]
                    attType = edgeAtt[1]
                else:
                    attName = edgeAtt
                    attType = 'string'

                outf.write('<key id="d%i" for="edge" attr.name="%s" attr.type="%s"></key>\n' %
                           (keyID, attName, attType))
                keyDict['edge'].append(keyID)
                keyID += 1

            outf.write('<graph id="G" edgedefault="%sdirected">\n' % 'un' if not self.directed else '')

            # Write nodes
            for node in self.nodes:
                outf.write(self.nodes[node].writeGraphlm(keyDict['node'], self.nodeAttNames))

            # Write edges
            for edge in self.edges:
                outf.write(edge.writeGraphlm(keyDict['edge'], self.edgeAttNames))

            # Finish
            outf.write('</graph>\n</graphml>')

            print('Nodes:\t%i' % len(self.nodes))
            print('Edges:\t%i' % len(self.edges))

    ########################################################################
    ##########################                    ##########################
    ##########################   ATTRIBUTE CORE   ##########################
    ##########################                    ##########################
    ########################################################################

    def attributeBasedCore(self, outfile, attributes, addseqs):

        # run through all subgraphs
        # - check the nodes attribute
        # - OPTIONAL check cliqueness
        # - remember the set of its value
        # - #

        # update on 10.11.2017
        # 1. Prepare cores for a single or set of attributes
        # 2. Collect all of the possible attributes variants -> go through the nodes and expand the possible options

        #   !!!!!!!!! only core and singletons !!!!!!!!!!

        # 3. Go through the subgraphs once again and check what variants of each attirbute are in the sugraph
        #    if all of the possible then add this subgraph to the core
        #    elif all of the variants are the same then add this to the singleton
        #    else ommit --- THIS IS THE TODO PART
        # 4. Write the summary for each attribute - make a subdirectory core_analysis/attributeX/core.txt, x_singletons... etc

        outdir =  outfile + '_attribute_core'

        # hold attributes in att_set
        attributes = attributes.split(',') if ',' in attributes else [attributes]
        att_set = {att: set() for att in attributes}
        # get each attributes collection genome: [genome1, genome2, ...]
        print('Getting attributes variants...')
        for subgraph_size in self.subGraphs:
            for subgraph in self.subGraphs[subgraph_size]:
                for node in subgraph:
                    #print('node')
                    for att in attributes:
                        #print(att, self.nodes[node].getAttVal(att))
                        att_set[att].add(self.nodes[node].getAttVal(att))

        # change sets to subdicts
        core = dict()
        for att in att_set:
            core[att] = {val: {} for val in att_set[att]}
            core[att]['core'] = {}

        # analyze the core
        print('Identifying the core set nodes for the attribute(s)...')
        for subgraph_size in sorted(self.subGraphs, reverse=True):
            for subgraph in sorted(self.subGraphs[subgraph_size]):
                # in each subgraph we check if it belongs to the core or singleton in each attribute
                tmp_att_set = {att: set() for att in attributes}

                for node in sorted(subgraph):
                    for att in attributes:
                        #print(self.nodes[node].getAttVal(att))
                        tmp_att_set[att].add(self.nodes[node].getAttVal(att))   

                # check the if sets are att_sets and tmp_att_sets are equal
                for att in attributes:
                    if att_set[att] == tmp_att_set[att]: # core
                        #print('Entered core')
                        try:
                            core[att]['core'][subgraph_size].append(subgraph)
                        except KeyError:
                            #print('Exception core')
                            #print(att, core[att]['core'])
                            core[att]['core'][subgraph_size] = [subgraph]
                            #print(att, core[att]['core'])
                            #input()
                    elif len(tmp_att_set[att]) == 1: # singleton
                        #print('Entered singleton')
                        #print(tmp_att_set[att], subgraph)
                        att_val = tmp_att_set[att].pop()
                        try:
                            core[att][att_val][subgraph_size].append(subgraph)
                        except KeyError:
                            #print('Exception singleton')
                            #print(att, core[att][att_val])
                            core[att][att_val][subgraph_size] = [subgraph]
                            #print(att, core[att][att_val])
                        #input()
                    else:
                        #print('Entered something between core and singleton')
                        pass # other suboptions
                    #print(sorted(list(att_set[att])))
                    #print(att_val)
                    #print(sorted(list(tmp_att_set[att])))


        for att in core:
            print(att)
            for att_val in core[att]:
                print('\t', att_val)
                print('\t\t', '\n\t\t- '.join(['%i: %i' %(key, len(core[att][att_val][key])) for key in  core[att][att_val].keys()]))
                #for size,val in core[att][att_val]:
                #   print('\t\t', size, '-', len(val))
        #exit()
        # read fasta file if needed
        seqs = {}
        if addseqs:
            print('Reading fasta file to add sequences to the output...')
            with open(addseqs) as inf:
                line = inf.readline()
                while line != '':
                    if '>' in line[0]:
                        header = line[1:].split()[0]
                        seqs[header] = ''
                    else:
                        seqs[header] += line.strip()

                    line = inf.readline()

        # make output directory and subdirectories
        print('Writing everything to the files...')
        if not path.isdir(outdir): makedirs(outdir)
        for att in core:
            print(att)
            if not path.isdir(path.join(outdir, att)): makedirs(path.join(outdir, att))

            # write core and singletons for attribute
            for att_val in core[att]:
                #print(subgraph_size)
                #print(core[att][subgraph_size])
                tmp = core[att][att_val]
                with open(path.join(outdir, att, att_val + '.fasta'), 'w') as outfasta:
                    with open(path.join(outdir, att, att_val + '.txt'), 'w') as outf:

                        outf.write('General summary\nSubGraphSize\tNumberOfSubGraphs\n')
                        for subGraphSize in sorted(tmp, reverse=True):
                            outf.write('%6d\t%6d\n' % (subGraphSize, len(tmp[subGraphSize])))

                        outf.write('\nDetailed summary\n\n')
                        for subGraphSize in sorted(tmp, reverse=True):

                            outf.write('----   %i   ----\n\n' % subGraphSize)
                            outf.write('NodeID\t%s%s\n' % ('\t'.join(self.nodeAttNames), '' if not addseqs else '\tSequence'))

                            separator = len(tmp[subGraphSize])
                            for subGraph in tmp[subGraphSize]:
                                separator -= 1

                                for node in subGraph:
                                    outf.write('%s\t%s%s\n' %
                                               (node,
                                                self.nodes[node].writeSummary(self.nodeAttNames),
                                                '' if not addseqs else '\t' + seqs[self.nodes[node].name]))
                                    outfasta.write('>%s\n%s\n' % (self.nodes[node].name, seqs[self.nodes[node].name]))
                                if separator != 0:
                                    outf.write('-------------\n')

                        outf.write('\n')


    def degree_based_sort(self, subgraph):

        """This function sorts nodes within subgraphs based on the degree so that the ones with highest degree could
        be placed either on the outer orbits or in the centre."""

        #TODO add sorting based on nodes degree
        return subgraph
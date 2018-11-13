import os
import glob
import networkx as nx


def get_files_in_directory(input_path, extension):
    # this function is used to return list of files in directory 'input_path'

    files = []
    for infile in glob.glob( os.path.join(str(input_path), u'*')):
        if os.path.isfile(infile) and infile.lower().endswith(extension):
          files.append(infile)
    return files


def load_el_network(
        input_file_path,
        weighted_graph=False,
        directed_graph=False,
        continue_invalid_lines=False
        ):

    # this function is used to load a network file. returns Graph
    # to graph. file content should be in format 'node \t node \t weight'
    # Information:
    # node name should be numeric or alphabetic characters. all spaces will
    # be replaced with 'underscore' character.
    # parameters:
    #   1- input_file_path: network data file. content should be in
    #      format 'node \t node \t weight'
    #   2- weighted_graph: if you want to graph be weighted, this parameter
    #      must be True. default is False.
    #   3- directed_graph: if you want to graph be directed, this parameter
    #      must be True. default is False.
    #   4- continue_invalid_lines: in the case that there are invalid lines
    #      in data file (invalid column count, invalid weights, ...), if this
    #      parameter be True, script will continue to load network. else script
    #      will exit.

    # graph definition

    all_nodes = []
    all_edges = []

    if directed_graph:
        g = nx.DiGraph()
    else:
        g = nx.Graph()

    # setting column count

    if weighted_graph:
        valid_column_count = 3
    else:
        valid_column_count = 2

    # importing data to network

    reading_nodes = False
    reading_edges = False

    inputfile = open(input_file_path, 'r')
    for line_number, l in enumerate(inputfile):
        line = l.strip()
        if line == '':
            continue

        if line == '#nodes':
            reading_nodes = True
            reading_edges = False
            continue
        elif line == '#edges':
            reading_edges = True
            reading_nodes = False
            continue

        ldata = l.split('\t')

        if reading_nodes:
            node = ldata[0].strip().lower().replace(' ','_')
            node = ''.join(ch for ch in node if ch.isalnum() or ch == '_')
            all_nodes.append(node)
            continue
        elif reading_edges:
            if len(ldata) < valid_column_count:
                print 'Error# Invalid line. Line No:' + str(line_number) + ', Line Data: ' + line
                if continue_invalid_lines:
                    continue
                else:
                    print 'Process ended due to previous Error'
                    exit()

            node1 = ldata[0].strip().lower().replace(' ','_')
            node1 = ''.join(ch for ch in node1 if ch.isalnum() or ch == '_')

            node2 = ldata[1].strip().lower().replace(' ','_')
            node2 = ''.join(ch for ch in node2 if ch.isalnum() or ch == '_')

            weight = 1
            if weighted_graph:
                try:
                    weight = float(ldata[2])
                except:
                    print 'Error# Invalid line. Line No:' + str(line_number) + ', Line Data: ' + line
                    if continue_invalid_lines:
                        continue
                    else:
                        print 'Process ended due to previous Error'
                        exit()

            if node1 != '' and node2 != '':
                if weighted_graph:
                    edge = (node1, node2, weight)
                    all_edges.append(edge)
                else:
                    edge = (node1, node2)
                    all_edges.append(edge)
            else:
                print 'Error# Invalid line. nodes are empty. Line No:' + str(line_number) + ', Line Data: ' + line
                if continue_invalid_lines:
                    continue
                else:
                    print 'Process ended due to previous Error'
                    exit()

    g.add_nodes_from(all_nodes)

    if weighted_graph:
        g.add_weighted_edges_from(all_edges)
    else:
        g.add_edges_from(all_edges)

    inputfile.close()
    return g


def read_graph(_file, extension, weighted=False, directed=False):
    # this function is used to read graph data from file and return graph

    if extension == '.graphml':
        graph = nx.read_graphml(_file)
    elif extension == '.gml':
        graph = nx.read_gml(_file)
    elif extension == '.net':
        graph = nx.read_pajek(_file)
    elif extension == '.el':
        graph = load_el_network(_file, weighted, directed)

    else:
        graph = None

    return graph


def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False
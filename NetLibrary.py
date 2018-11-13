import networkx as nx
import matplotlib
import numpy as np
import os
import pandas as pd
from scipy import stats

matplotlib.use('Agg')
import matplotlib.pyplot as plt

global _nn
global _ne
global _lcc_size
global _adc
global _abc
global _acc
global _alc
global _acoc
global _ncc
global _accs
global _trans
global _dens
global _max_cs
global _dac
global _adcon
global _acfcc
global _nec
global _akc
global _nocc
global _nd
global _ae
global _rad
global _ap
global _asp


def draw_network(graph, network_name):
    # this function is used to save network visualization to a file with
    # name 'network_name'

    nx.draw_networkx(
        graph,
        pos=nx.shell_layout(graph),
        with_labels=False,
        node_color='red',
        font_size=8,
        node_size=50,
        edge_color='gray')

    plt.savefig(network_name)


def compute_network_features(graph, network_name):
    # this function is used to calculate network features and returns result
    # as a dictionary: Dict<feature_name,feature_value>
    # --------------------------------------------------------------------------------

    network_features = dict()
    node_features_list = list()

    netclass = network_name.split('___')

    if len(netclass) > 1:
        network_features['group'] = netclass[0]
        network_features['Network Name'] = netclass[1]
    else:
        network_features['group'] = '_unknown_'
        network_features['Network Name'] = network_name

    if graph.is_directed():
        network_features['Is Directed?'] = True
    else:
        network_features['Is Directed?'] = False

    if graph.is_multigraph():
        network_features['Is MultiGraph?'] = True
    else:
        network_features['Is MultiGraph?'] = False

    # Global Attributes
    # --------------------------------------------------------------------------------
    # number of nodes
    # --------------------------------------------------------------------------------
    if _nn:
        try:
            nodes_count = nx.number_of_nodes(graph)
            network_features['Number of Nodes'] = nodes_count
        except:
            network_features['Number of Nodes'] = 'NA'

    # number of edges
    # --------------------------------------------------------------------------------
    if _ne:
        try:
            edges_count = nx.number_of_edges(graph)
            network_features['Number of Edges'] = edges_count
        except:
            network_features['Number of Edges'] = 'NA'

    # network density
    # --------------------------------------------------------------------------------
    if _dens:
        try:
            density = nx.density(graph)
            network_features['Density'] = density
        except:
            network_features['Density'] = 'NA'

    # graph degree assortativity
    # --------------------------------------------------------------------------------
    if _dac:
        try:
            graph_degree_assortativity = nx.degree_assortativity_coefficient(graph)
            network_features['Graph Degree Assortativity'] = graph_degree_assortativity
        except:
            network_features['Graph Degree Assortativity'] = 'NA'

    # avg. closeness centrality
    # --------------------------------------------------------------------------------
    if _acc:
        try:
            ccn = nx.closeness_centrality(graph)
            mccn = np.mean(ccn.values())
            network_features['Avg. Closeness Centrality'] = mccn
        except:
            network_features['Avg. Closeness Centrality'] = 'NA'

    # avg. betweenness centrality
    # --------------------------------------------------------------------------------
    if _abc:
        try:
            bcn = nx.betweenness_centrality(graph)
            mbcn = np.mean(bcn.values())
            network_features['Avg. Betweenness Centrality'] = mbcn
        except:
            network_features['Avg. Betweenness Centrality'] = 'NA'

    # avg. degree centrality
    # --------------------------------------------------------------------------------
    if _adc:
        try:
            dcn = nx.degree_centrality(graph)
            mdcn = np.mean(dcn.values())
            network_features['Avg. Degree Centrality'] = mdcn
        except:
            network_features['Avg. Degree Centrality'] = 'NA'

    # avg. degree connectivity
    # --------------------------------------------------------------------------------
    if _adcon:
        try:
            dc = nx.average_degree_connectivity(graph)
            adc = np.mean(dc.values())
            network_features['Avg. Degree Connectivity'] = adc
        except:
            network_features['Avg. Degree Connectivity'] = 'NA'

    # avg. load centrality
    # --------------------------------------------------------------------------------
    if _alc:
        try:
            lc = nx.load_centrality(graph)
            mlc = np.mean(lc.values())
            network_features['Avg. Load Centrality'] = mlc
        except:
            network_features['Avg. Load Centrality'] = 'NA'

    # avg. edge betweenness centrality
    # --------------------------------------------------------------------------------

    # try:
    #     ebc = nx.edge_betweenness_centrality(graph)
    #     mebc = np.mean(ebc.values())
    #     network_features['Avg. Edge Betweenness centrality'] = mebc
    # except:
    #     network_features['Avg. Edge Betweenness centrality'] = 'NA'

    # edge connectivity
    # --------------------------------------------------------------------------------
    # try:
    #     ec = nx.edge_connectivity(graph)
    #     network_features['Edge Connectivity'] = ec
    # except:
    #     network_features['Edge Connectivity'] = 'NA'

    # diameter
    # --------------------------------------------------------------------------------
    if _nd:
        try:
            diameter = nx.diameter(graph)
            network_features['Diameter'] = diameter
        except:
            network_features['Diameter'] = 'NA'

    # eccentricity
    # --------------------------------------------------------------------------------
    if _ae:
        try:
            eccentricity = nx.eccentricity(graph)
            network_features['Avg. Eccentricity'] = np.mean(eccentricity.values())
        except:
            network_features['Eccentricity'] = 'NA'

    # radius
    # --------------------------------------------------------------------------------
    if _rad:
        try:
            radius = nx.radius(graph)
            network_features['Radius'] = radius
        except:
            network_features['Radius'] = 'NA'

    # Non MultiGraph Features
    # --------------------------------------------------------------------------------
    if not graph.is_multigraph():

        # transitivity
        # ----------------------------------------------------------------------------
        if _trans:
            try:
                transitivity = nx.transitivity(graph)
                network_features['Transitivity'] = transitivity
            except:
                network_features['Transitivity'] = 'NA'

        # Katz centrality
        # ----------------------------------------------------------------------------
        if _akc:
            try:
                katz = nx.katz_centrality(graph)
                mean_katz = np.mean(katz.values())
                network_features['Avg. Katz Centrality'] = mean_katz
            except:
                network_features['Avg. Katz Centrality'] = 'NA'

        # PageRank
        # ----------------------------------------------------------------------------
        if _ap:
            try:
                pagerank = nx.pagerank(graph)
                mean_pagerank = np.mean(pagerank.values())
                network_features['Avg. PageRank'] = mean_pagerank
            except:
                network_features['Avg. PageRank'] = 'NA'

    # Undirected Graphs
    # --------------------------------------------------------------------------------

    if not nx.is_directed(graph):

        # Degree
        # ----------------------------------------------------------------------------
        #
        # try:
        #     all_degrees = nx.degree(graph)
        #     mean_degrees = np.mean(all_degrees.values())
        #     network_features['Avg. Degree'] = mean_degrees
        # except:
        #     network_features['Avg. Degree'] = 'NA'

        # connected components
        # ----------------------------------------------------------------------------
        if _nocc:
            try:
                cc_number = nx.number_connected_components(graph)
                network_features['Number of Connected Components'] = cc_number
            except:
                network_features['Number of Connected Components'] = 'NA'

        # lcc size fraction && avg. cc size
        # ----------------------------------------------------------------------------
        if _accs or _lcc_size:
            try:
                cc_list = list(nx.connected_components(graph))
                cc_sizes = []
                for cc in cc_list:
                    cc_sizes.append(len(cc))

                lcc_size = np.max(cc_sizes)
                if _accs:
                    network_features['lcc_size_fraction'] = lcc_size / float(nodes_count)
                if _lcc_size:
                    mean_cc_sizes = np.mean(cc_sizes)
                    network_features['Avg. Connected Component Size'] = mean_cc_sizes
            except:
                if _accs:
                    network_features['lcc_size_fraction'] = 'NA'
                if _lcc_size:
                    network_features['Avg. Connected Component Size'] = 'NA'

        # communicability centrality for Undirected networks
        # ----------------------------------------------------------------------------
        if not graph.is_multigraph():
            if _acoc:
                try:
                    cc = nx.communicability_centrality(graph)
                    mcc = np.mean(cc.values())
                    network_features['Avg. Communicability Centrality'] = mcc
                except:
                    network_features['Avg. Communicability Centrality'] = 'NA'

            # clustering coefficient
            # -------------------------------------------------------------------------
            if _ncc:
                try:
                    clustering_coefficient = nx.average_clustering(graph)
                    network_features['Network Clustering Coefficient'] = clustering_coefficient
                except:
                    network_features['Network Clustering Coefficient'] = 'NA'

        # clique analysis for Undirected networks
        # ----------------------------------------------------------------------------
        if _max_cs:
            try:
                cliques_obj = nx.find_cliques(graph)
                cliques = [clq for clq in cliques_obj]

                clique_sizes = []
                for c in cliques:
                    clique_sizes.append(len(c))

                # user_clique_size = 5
                if len(clique_sizes) > 0:
                    # network_features['No of Cliques with size ' + str(user_clique_size)] \
                    # = clique_sizes.count(user_clique_size)
                    network_features['Avg. Clique Size'] = np.mean(clique_sizes)
                    network_features['Max Clique Size'] = np.max(clique_sizes)
                else:
                    # network_features['No of Cliques with size ' + str(user_clique_size)] = 0
                    network_features['Avg. Clique Size'] = 0
                    network_features['Max Clique Size'] = 0
            except:
                # network_features['No of Cliques with size ' + str(user_clique_size)] = 'NA'
                network_features['Avg. Clique Size'] = 'NA'
                network_features['Max Clique Size'] = 'NA'

                # else:
                # try:
                #     all_in_degrees = nx.DiGraph.in_degree(graph)
                #     all_out_degrees = nx.DiGraph.out_degree(graph)
                #
                #     mean_in_degrees = np.mean(all_in_degrees.values())
                #     mean_out_degrees = np.mean(all_out_degrees.values())
                #
                #     network_features['Avg. In Degree'] = mean_in_degrees
                #     network_features['Ave. Out Degree'] = mean_out_degrees
                # except:
                #     network_features['Avg. In Degree'] = 'NA'
                #     network_features['Ave. Out Degree'] = 'NA'

    # Nodes Features Calculation

    for node in graph.nodes():

        node_features = dict()

        try:
            node_features['group'] = network_name
        except:
            node_features['group'] = 'NA'

        if _abc:
            try:
                node_features['Betweenness Centrality'] = bcn[node]
            except:
                node_features['Betweenness Centrality'] = 'NA'

        if _acc:
            try:
                node_features['Closeness Centrality'] = ccn[node]
            except:
                node_features['Closeness Centrality'] = 'NA'

        if _adc:
            try:
                node_features['Degree Centrality'] = dcn[node]
            except:
                node_features['Degree Centrality'] = 'NA'

        if _alc:
            try:
                node_features['Load Centrality'] = lc[node]
            except:
                node_features['Load Centrality'] = 'NA'

        if _ae:
            try:
                node_features['Eccentricity'] = eccentricity[node]
            except:
                node_features['Eccentricity'] = 'NA'

        if not graph.is_multigraph():
            if _akc:
                try:
                    node_features['Katz Centrality'] = katz[node]
                except:
                    node_features['Katz Centrality'] = 'NA'

            if _ap:
                try:
                    node_features['PageRank'] = pagerank[node]
                except:
                    node_features['PageRank'] = 'NA'

        if not nx.is_directed(graph):
            # try:
            #     node_features['Degree'] = all_degrees[node]
            # except:
            #     node_features['Degree'] = 'NA'

            if not graph.is_multigraph():
                if _acoc:
                    try:
                        node_features['Communicability Centrality'] = cc[node]
                    except:
                        node_features['Communicability Centrality'] = 'NA'
                        # else:
                        # try:
                        #     node_features['In Degree'] = all_in_degrees[node]
                        # except:
                        #     node_features['In Degree'] = 'NA'
                        #
                        # try:
                        #     node_features['Out Degree'] = all_out_degrees[node]
                        # except:
                        #     node_features['Out Degree'] = 'NA'

        node_features_list.append(node_features)

    return network_features, node_features_list


def statistic_analysis(input_tsv_file, output_directory, draw_boxplot=True):
    _file = open(input_tsv_file, 'r')
    header_data = _file.readline().split('\t')
    _file.close()

    try:
        header_data.remove('group')
    except:
        print 'Specified file does not contain \'group\' attribute'
        return

    header_items = [item.strip() for item in header_data]

    data = pd.read_csv(input_tsv_file, sep='\t')

    if draw_boxplot:
        plt.rcParams['font.size'] = 18

        for item in header_items:
            data.boxplot(item, by='group', figsize=(20, 13))

            locs, labels = plt.xticks()
            plt.setp(labels, rotation=90, fontsize=14)
            plt.savefig(os.getcwd()
                        + os.path.sep
                        + output_directory
                        + os.path.sep
                        + item
                        + '.png')
            plt.close('all')

    result_file = open(output_directory + os.path.sep + 'statistic_result.txt', 'w')

    grps = pd.unique(data.group.values)

    for item in header_items:
        result_file.write('\n=============================================================================\n')
        result_file.write('============================== basic statistics =============================\n')
        for i, g in enumerate(grps):
            d = pd.Series(data[item][data.group == g].values).describe()
            result_file.write(str('Group: ' + g + ', Feature: ' + item + '\nDescription:\n' + str(d) + '\n'))
            if i < (len(grps) - 1):
                result_file.write('\n-----------------------------------------------------------------------------\n')

        result_file.write('\n=============================================================================\n')
        result_file.write('============================== ANOVA Test result ============================\n\n')

        d_data = [data[item][data.group == grp].values for grp in grps]
        result_file.write('For feature: ' + item + '\n')
        f, p = stats.f_oneway(*d_data)
        result_file.write('F-value: ' + str(f) + '\nP-value: ' + str(p) + '\n')

    result_file.close()

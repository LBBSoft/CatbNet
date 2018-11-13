from __future__ import division
import NetLibrary
from PyQt4.QtGui import QMainWindow, QFileDialog
from UI import Window
from PyQt4 import QtGui
import Utilities
import shutil, time, os, math, threading, sys


class Application(QMainWindow, Window):
    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.setupUi(self)
        self.pbar.setValue(0)
        self.LblWorking.hide()
        self.LblWorking.setText("Working ...")

    def on_pushButton_released(self):
        self.LblNetworksDirectory.setText(str(QFileDialog.getExistingDirectory(self, "Select Directory")))

    def close_application(self):
        self.close()

    def on_BtnStart_released(self):
        self.pbar.setValue(0)
        # check app conditions
        if not self.RdoDirected.isChecked() and not self.RdoUndirected.isChecked():
            w = QtGui.QWidget()
            QtGui.QMessageBox.information(w, "Message", "Choose your networks are directed or undirected, please")
            return
        if self.LblNetworksDirectory.text() == '':
            w = QtGui.QWidget()
            QtGui.QMessageBox.information(w, "Message", "Choose network files directory, please")
            return
        if not self.RdoEdgeList.isChecked() and not self.RdoGML.isChecked() and not self.RdoGraphML.isChecked() and not self.RdoPajek.isChecked():
            w = QtGui.QWidget()
            QtGui.QMessageBox.information(w, "Message", "Choose network files format, please")
            return

        NetLibrary._nn = self.chkNodeCount.isChecked()
        NetLibrary._ne = self.chkEdgeCount.isChecked()
        NetLibrary._lcc_size = self.ChkLccSize.isChecked()
        NetLibrary._adc = self.ChkMeanDegreeCentrality.isChecked()
        NetLibrary._abc = self.ChkMeanBetCentrality.isChecked()
        NetLibrary._acc = self.ChkMeanCloCentrality.isChecked()
        NetLibrary._alc = self.ChkMeanLoadCentrality.isChecked()
        NetLibrary._acoc = self.ChkMeanCommCentrality.isChecked()
        NetLibrary._ncc = self.ChkClusteringCoefficient.isChecked()
        NetLibrary._accs = self.ChkMeanCCSize.isChecked()
        NetLibrary._trans = self.ChkTransitivity.isChecked()
        NetLibrary._dens = self.ChkDenisty.isChecked()
        NetLibrary._max_cs = self.ChkMaxCliqueSize.isChecked()
        NetLibrary._dac = self.ChkDegreeAssortativityCoefficient.isChecked()
        NetLibrary._adcon = self.ChkAvgDegreeConnectivity.isChecked()
        NetLibrary._acfcc = self.ChkAvgCFClosCentrality.isChecked()
        NetLibrary._nec = self.ChkEigenvectorCentrality.isChecked()
        NetLibrary._akc = self.ChkKatzCentrality.isChecked()
        NetLibrary._nocc = self.ChkNumberOfConnectedComponents.isChecked()
        NetLibrary._nd = self.ChkNetDiameter.isChecked()
        NetLibrary._ae = self.ChkAvgEccentricity.isChecked()
        NetLibrary._rad = self.ChkNetworkRadius.isChecked()
        NetLibrary._ap = self.ChkAvgPagerank.isChecked()
        NetLibrary._asp = self.ChkAvgShortestPath.isChecked()


        self.BtnStart.setDisabled(True)
        self.LblWorking.show()
        self.LblWorking.setText("Working ...")
        t = threading.Thread(target=self.worker)
        t.start()

    def worker(self):
        # region getting network format

        if self.RdoEdgeList.isChecked():
            ext = '.el'
        elif self.RdoGraphML.isChecked():
            ext = '.graphml'
        elif self.RdoGML.isChecked():
            ext = '.gml'
        elif self.RdoPajek.isChecked():
            ext = '.net'
        # endregion

        try:
            if os.path.exists('Output_Directory'):
                shutil.rmtree('Output_Directory')

            if not os.path.exists('Output_Directory'):
                os.makedirs('Output_Directory')
            if not os.path.exists('Output_Directory' + os.path.sep + 'bp_networks'):
                os.makedirs('Output_Directory' + os.path.sep + 'bp_networks')
            if not os.path.exists('Output_Directory' + os.path.sep + 'bp_nodes'):
                os.makedirs('Output_Directory' + os.path.sep + 'bp_nodes')

        except:
            print("Please remove 'Output_Directory' manually and try again.")
            QtGui.QApplication.processEvents()
            self.pbar.setValue(100)
            self.BtnStart.setDisabled(False)
            self.LblWorking.setText("Failed!")
            sys.exit(-1)


        # region resolving network settings

        directed = self.RdoDirected.isChecked()
        weighted = self.ChkWeighted.isChecked()
        # endregion

        # region loading files & networks

        files = Utilities.get_files_in_directory(self.LblNetworksDirectory.text(), ext)
        if len(files) == 0:
            print("not found file with extension:" + str(ext))
            QtGui.QApplication.processEvents()
            self.pbar.setValue(100)
            self.BtnStart.setDisabled(False)
            self.LblWorking.setText("Failed!")
            sys.exit(-1)

        all_results = []
        all_nodes_results = []
        for _file in files:
            graph = Utilities.read_graph(_file, ext, weighted, directed)
            if graph is None:
                print "for file: " + str(_file) + "graph in None."
                continue

            network_name = os.path.basename(_file)

            # region feature calculation
            # -----------------------------------------------------------------------------------
            results, nodes_results = NetLibrary.compute_network_features(graph, network_name)
            all_results.append(results)
            all_nodes_results.append(nodes_results)
            # endregion

            # region doing progress
            # -----------------------------------------------------------------------------------

            QtGui.QApplication.processEvents()
            val = math.ceil(self.pbar.value() + 1 / len(files) * 100)
            if val >= 100:
                val = 95
            self.pbar.setValue(val)
            # endregion

        # endregion

        # region writing full network features

        if len(all_results) > 0:
            common_items = []
            for i, lst in enumerate(all_results):
                if i is 0:
                    common_items = lst.keys()
                else:
                    tmp1 = set(common_items)
                    tmp2 = set(lst.keys())
                    common_items = list(tmp1.intersection(tmp2))

        common_items.sort()
        tsv = open('Output_Directory' + os.path.sep + 'network_results.tsv', 'w')

        predefined_headers = ['group', 'Network Name', 'Is Directed?', 'Is MultiGraph?', 'Number of Nodes',
                              'Number of Edges']
        QtGui.QApplication.processEvents()

        for it, item in enumerate(all_results):
            if it is 0:
                tsv.write('group')
                tsv.write('\tNetwork Name')
                tsv.write('\tIs Directed?')
                tsv.write('\tIs MultiGraph?')
                tsv.write('\tNumber of Nodes')
                tsv.write('\tNumber of Edges')

                for header in common_items:
                    if header not in predefined_headers:
                        tsv.write(str('\t' + header))

                tsv.write('\n')

                tsv.write(str(item['group']))
                tsv.write('\t' + str(item['Network Name']))
                tsv.write('\t' + str(item['Is Directed?']))
                tsv.write('\t' + str(item['Is MultiGraph?']))
                tsv.write('\t' + str(item['Number of Nodes']))
                tsv.write('\t' + str(item['Number of Edges']))

                for header in common_items:
                    if header not in predefined_headers:
                        tsv.write(str('\t' + str(item[header])))

                tsv.write('\n')

            else:
                tsv.write(str(item['group']))
                tsv.write('\t' + str(item['Network Name']))
                tsv.write('\t' + str(item['Is Directed?']))
                tsv.write('\t' + str(item['Is MultiGraph?']))
                tsv.write('\t' + str(item['Number of Nodes']))
                tsv.write('\t' + str(item['Number of Edges']))

                for header in common_items:
                    if header not in predefined_headers:
                        tsv.write(str('\t' + str(item[header])))

                tsv.write('\n')
        tsv.close()
        # endregion

        # region writing network features for analysis

        bef = open('Output_Directory' + os.path.sep + 'network_results.tsv', 'r')
        new = open('Output_Directory' + os.path.sep + 'temp.tsv', 'w')

        for line in bef:
            data = line.split('\t')
            data[1:4] = []
            newline = '\t'.join(data)
            new.write(newline)
        bef.close()
        new.close()

        # endregion

        # region writing node features
        # ---------------------------------------------------------------------------------------

        if len(all_nodes_results) > 0:
            common_items = []
            for l, net_item in enumerate(all_nodes_results):
                for j, node_item in enumerate(net_item):
                    if j is 0 and l is 0:
                        common_items = node_item.keys()
                    else:
                        tmp1 = set(common_items)
                        tmp2 = set(node_item.keys())
                        common_items = list(tmp1.intersection(tmp2))

        common_items.sort()

        tsv = open('Output_Directory' + os.path.sep + 'nodes_result.tsv', 'w')

        predefined_headers = ['group']

        for k, net_item in enumerate(all_nodes_results):
            for it, item in enumerate(net_item):
                if it is 0 and k is 0:
                    tsv.write('group')

                    for header in common_items:
                        if header not in predefined_headers:
                            tsv.write(str('\t' + header))

                    tsv.write('\n')

                    tsv.write(str(item['group']))

                    for header in common_items:
                        if header not in predefined_headers:
                            tsv.write(str('\t' + str(item[header])))

                    tsv.write('\n')

                else:
                    tsv.write(str(item['group']))

                    for header in common_items:
                        if header not in predefined_headers:
                            tsv.write(str('\t' + str(item[header])))

                    tsv.write('\n')
        tsv.close()
        # endregion

        # region drawing bax plots

        # NetLibrary.draw_box_plots_for_groups('Output_Directory/temp.tsv',
        #                                      'bp_networks')
        # NetLibrary.draw_box_plots_for_groups('Output_Directory/nodes_result.tsv',
        #                                      'bp_nodes')
        # endregion

        # region statistical progress

        if self.ChkStatisticAnalysisOfAllNetworks.isChecked():
            NetLibrary.statistic_analysis('Output_Directory' + os.path.sep + 'nodes_result.tsv'
                                          , 'Output_Directory' + os.path.sep + 'bp_nodes',
                                          True)

        if self.ChkStatisticAnalysisOfGroups.isChecked() and self.ChkFilesAreGrouped.isChecked():
            NetLibrary.statistic_analysis('Output_Directory' + os.path.sep + 'temp.tsv'
                                          , 'Output_Directory' + os.path.sep + 'bp_networks',
                                          True)

        # endregion

        # region doing progress
        QtGui.QApplication.processEvents()
        self.pbar.setValue(100)
        self.BtnStart.setDisabled(False)
        os.remove('Output_Directory' + os.path.sep + 'temp.tsv')

        self.LblWorking.setText("Success!")

        # endregion

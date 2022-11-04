#! python3

import numpy as np
import os
import pandas as pd
from ast import literal_eval
from sympy import Symbol
from plot_network_better import *
from pathlib import Path

class GetNetworkData:
    
    def __init__(self, name, path=None):
        if not path:
            self.path = r"..\data\{}".format(name)

        CWD = os.getcwd()
        os.chdir(self.path)
        files = os.listdir()
        self.net_file = [Path(file).absolute() for file in files if file.endswith('_net.tntp')][0]
        try:
            self.node_file = [Path(file).absolute() for file in files if file.endswith('_node.tntp')][0]
        except IndexError:
            self.node_file = None
        self.trips_file = [Path(file).absolute() for file in files if file.endswith('_trips.tntp')][0]
        os.chdir(CWD)
        
        self.net_name = name
        self.net = self.__import_net()
        self.net.index = pd.MultiIndex.from_tuples(list(zip(self.net.init_node, self.net.term_node)))

    def __import_net(self):
        """Returns the entire network."""
        with open(self.net_file) as netfile:
            net = pd.read_csv(netfile, skiprows=8, sep='\t')
            net.columns = [s.strip().lower() for s in net.columns]
            net.drop(['~', ';'], axis=1, inplace=True)
        return net
    
    @property
    def links(self):
        """Returns all the links present in the network."""
        return list(zip(self.net.init_node, self.net.term_node))

    @property
    def lp_func(self):
        """Returns the link performance functions for the links in the network."""
        lp_funcs = {}
        for link in self.links:
            links_data = self.net.loc[link]
            t_0 = links_data['free_flow_time']
            b = links_data['b']
            x = Symbol('x')
            capacity = links_data['capacity']
            power = links_data['power']
            # Link travel time = free flow time * (1 + B*(flow/capacity)**Power).
            lpf = f'{t_0}*(1 + {b}*((x/{capacity})**{power}))'
            lp_funcs[link] = lpf
        return lp_funcs
    
    @property
    def demands(self):
        """Returns the demands of the from origin nodes to destination nodes
        for the network in pandas dataframe object."""
        demand = {}
        with open(self.trips_file) as file:
            text = file.read()
            origins_data = text.split('Origin')[1:]    
            for data in origins_data:
                row = data.replace(' ', '').replace(';', ',').strip().split('\n')
                origin = literal_eval(row[0])
                origin_demand = literal_eval('{' + ''.join(row[1:]) + '}')
                demand[origin] = origin_demand
            df = pd.DataFrame(demand)
        return df

    @property
    def node_positions(self):
        """Returns the (x, y) coordinates of the Network."""
        pos = {}
        try:
            with open(self.node_file) as file:
                for line in file:
                    node, x, y = line.strip().replace(';', '').strip('\t').split('\t')
                    if node.isdigit() and self.is_number(x) and self.is_number(y):
                        pos[literal_eval(node)] = (literal_eval(x), literal_eval(y))
        except (FileNotFoundError, TypeError) as e:
            for node in self.net.init_node:
                pos[node] = tuple(np.random.randint(low=-990, high=991, size=2))
        return pos

    @staticmethod
    def is_number(string):
        try:
            return True if float(string) else False
        except ValueError:
            return False

    def plotNet(self):
        f, ax = plt.subplots(figsize=(16,9))
        ax = plot_net(links=self.links, xx=None, tt=None, title=self.net_name, pos=self.node_positions)
        return ax

if __name__ == '__main__':

    NETWORK = 'SiouxFalls'
    SiouxFalls = GetNetworkData(NETWORK)
    demands = SiouxFalls.demands
    positions = SiouxFalls.node_positions
    links = SiouxFalls.links
    lp_func = SiouxFalls.lp_func

##    SiouxFalls.plotNet()
##    plt.show()

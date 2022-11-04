#! python3

import numpy as np
import os
import pandas as pd
##import openmatrix as omx
import sys
from pprint import pprint
from ast import literal_eval
from sympy import Symbol

class GetNetworkData:
    
    def __init__(self, name, path=None):
        if not path:
            self.path = r"H:\BOOKS\CLASS\Master's\Code\TransportationNetworks-master\{}".format(name)
        self.net_name = name
        self.net = self.__import_net()
        self.net.index = pd.MultiIndex.from_tuples(list(zip(self.net.init_node, self.net.term_node)))

    def __import_net(self):
        """Returns the entire network."""
        
        with open(os.path.join(self.path, f'{self.net_name}_net.tntp')) as netfile:
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
            lpf = str(t_0 * (1 + b*(x/capacity)**power))
            lp_funcs[link] = lpf
            
        # Link generalized cost = Link travel time + toll_factor * toll + distance_factor * distance
        # Not Implemented
        
        return lp_funcs
    
    @property
    def demands(self):
        """Returns the demands of the from origin nodes to destination nodes
        for the network in pandas dataframe object."""
        demand = {}
        with open(os.path.join(self.path, f'{self.net_name}_trips.tntp')) as file:
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
        with open(os.path.join(self.path, f'{self.net_name}_node.tntp')) as file:
            for line in file:
                node, x, y = line.strip().replace(';', '').strip('\t').split('\t')
                if node.isdigit() and self.is_number(x) and self.is_number(y):
                    pos[literal_eval(node)] = (literal_eval(x), literal_eval(y))
        return pos

    @staticmethod
    def is_number(string):
        try:
            return True if float(string) else False
        except ValueError:
            return False


if __name__ == '__main__':
    NETWORK = 'SiouxFalls'
    SiouxFalls = GetNetworkData(NETWORK)
    demands = SiouxFalls.demands
    positions = SiouxFalls.node_positions
    links = SiouxFalls.links
    lp_func = SiouxFalls.lp_func
    
    

    

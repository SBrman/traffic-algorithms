#! python3 

__author__ = "Simanta Barman"
__email__ = "barma017@umn.edu"

from Network import *
from Shortest_path import *
import math
from pprint import pprint, pformat

logging.disable(logging.INFO)
logging.basicConfig(level=logging.DEBUG, format="%(message)s")

class STOCH_loader:
    def __init__(self, link_labels, demand, theta=None):
        self.A = list(link_labels.keys())
        self.link_labels = link_labels
        self.demand = demand
        self.theta = theta if theta else 1

    def __get_node_labels(self, origin):
        """Returns node labels, backnode vector as Dicts using dijkstra\'s shortest path algorithm.
        Inputs: origin, link_labels"""
        labels, backnode = Shortest_path(labels=self.link_labels, start=origin).dijkstra()
        labels = {node: label for node, label in labels.items() if label < math.inf}
        
        return labels

    def __get_reasonable_links(self, origin, node_labels=None):
        if not node_labels:
            node_labels, backnode = self.__get_node_labels(origin)
        reasonable_links = {(i, j) for (i, j) in self.A if {i, j} <= set(node_labels.keys())
                                                        and node_labels[i] < node_labels[j]}
        return reasonable_links

    def __create_bush(self, links):
        """Returns Bush which is created using the reasonable links only."""
        lp_funcs = {link: link_label for link, link_label in self.link_labels.items() if link in links}
        G = Network()
        G.add_links(lp_funcs.keys())
        G.add_lpf(lp_funcs)
        return G
    
    def load_STOCH(self, origin):
        """Returns the link flows for a given origin to the destination."""
        # Step 1: Getting reasonable links
        node_L = self.__get_node_labels(origin)
        reasonable_links = self.__get_reasonable_links(origin, node_labels=node_L)
        Bush = self.__create_bush(links=reasonable_links)

        link_likelihood = {(i, j): round(math.exp(self.theta * (node_L[j] - node_L[i] - t_ij)), 6)
                  for (i, j), t_ij in self.link_labels.items() if (i, j) in reasonable_links}

        # Step 2: Calculating node and link weights in forward topological order
        # Getting topsort from node_labels from dijkstra
        _, topsorted_nodes = zip(*sorted(zip(node_L.values(), node_L.keys())))

        node_weight = {}
        link_weight = {}
        for i in topsorted_nodes:
            # 2.a Setting the current node i in topological order and setting Wi = 1 
            node_weight[i] = 1
            for (i, j) in Bush.forward_star(i):
                link_weight[(i, j)] = node_weight[i] * link_likelihood[(i, j)]
            node_weight[i] = sum(link_weight[(k, i)] for (k, i) in Bush.reverse_star(i))
                
        # Step 3: Calculating node and link weights in reverse topological order after loading demand
        link_flows = {}
        node_link_flow = {}
        for j in topsorted_nodes[::-1]:
            ########################################################################################
            ## Shouldn't this demand be the sum of all demand that comes to the destination node? ##
            ## No, because STOCH loading is done for every origin then the link_flows sum is used ##
            ########################################################################################
            node_link_flow[j] = self.demand.get((origin, j), 0) + sum(link_flows.get((j, k), 0) 
                                                            for (j, k) in Bush.forward_star(j))

            for (i, j) in Bush.reverse_star(j):
                link_flows[(i, j)] = round(node_link_flow[j] * (link_weight[(i, j)] / node_weight[j]), 6)

        return link_flows
    
    def load(self):
        link_flows = {link: 0 for link in self.A}
        for origin, _ in self.demand.keys():
            x_ij = self.load_STOCH(origin=origin)
            for link, x in x_ij.items():
                link_flows[link] = link_flows[link] + x

        return link_flows
    

if __name__ == "__main__":
    lp_funcs = {(1, 2): 2, (1, 3): 1, (2, 3): 2, (3, 2): 2, (2, 4): 1, (3, 4): 3}
    demand = {(1, 4): 100}
    
    STOCH_loader = STOCH_loader(lp_funcs, demand)
    xij = STOCH_loader.load()
    print(f'Link flows = \n{pformat(xij)}')

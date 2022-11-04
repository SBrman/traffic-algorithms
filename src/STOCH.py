#! python3 

__author__ = "Simanta Barman"
__email__ = "barma017@umn.edu"

from Network import *
from Shortest_path import *
import math
from pprint import pprint, pformat

logging.disable(logging.INFO)
logging.basicConfig(level=logging.DEBUG, format="%(message)s")

class Tap_SUE:
    def __init__(self, Graph, demand):
        self.G = Graph
        self.demand = demand
        self.origins = tuple({node for node, _ in demand.keys()})
        self.destinations = tuple({node for _, node in demand.keys()})

    def __get_node_labels(self, origin, link_labels):
        """Returns node labels, backnode vector as Dicts using dijkstra\'s shortest path algorithm.
        Inputs: origin, link_labels"""

        labels, backnode = Shortest_path(labels=link_labels, start=origin).dijkstra()
        return labels

    def __get_reasonable_links(self, origin, node_labels=None):
        if not node_labels:
            link_labels = {link: tij(0) for link, tij in self.G.lp_func.items()}
            node_labels, backnode = Shortest_path(labels=link_labels, start=origin).dijkstra()
        reasonable_links = {(i, j) for (i, j) in self.G.A if node_labels[i] < node_labels[j]}
        return reasonable_links

    def __create_bush(self, links, link_labels):
        """Returns Bush which is created using the reasonable links only."""
        
        lp_funcs = {link: lp_func for link, lp_func in link_labels.items() if link in links}
        return create_network(links=links, lp_funcs=lp_funcs)
    
    def load_STOCH(self, origin, fixed_link_labels=None, theta=1):

        if not fixed_link_labels:
            fixed_link_labels = {link: tij(0) for link, tij in self.G.lp_func.items()}
            
        # Step 1: Getting reasonable links
        node_L = self.__get_node_labels(origin, fixed_link_labels)
        reasonable_links = self.__get_reasonable_links(origin, node_labels=node_L)
        Bush = self.__create_bush(links=reasonable_links, link_labels=fixed_link_labels)

        logging.warning(f'Step 1:\n{"-"*79}\n{"-"*79}')
        logging.warning(f'L = {node_L}\n\nReasonable Links = {reasonable_links}')

        link_likelihood = {(i, j): round(math.exp(theta * (node_L[j] - node_L[i] - t_ij)), 6)
                  for (i, j), t_ij in fixed_link_labels.items() if (i, j) in reasonable_links}

        logging.warning(f'\nLink Likelihood = \n{pformat(link_likelihood, indent=4, width=20)}\n')
        logging.warning(f'{"-"*79}\n')


        # Step 2: Calculating node and link weights in forward topological order
        # Getting topsort from node_labels from dijkstra
        _, topsorted_node = zip(*sorted(zip(node_L.values(), node_L.keys())))
        logging.warning(f'\n\nStep 2:\n{"-"*79}\n{"-"*79}')
        logging.warning(f'Topological Sort = {topsorted_node}\n')

        node_weight = {}
        link_weight = {}
        for i in topsorted_node:
            # 2.a Setting the current node i in topological order and setting Wi = 1 
            node_weight[i] = 1
            logging.warning(f'Node-{i}: \n{" "*8}W_{i} = {node_weight[i]}{" "*8}(Initializing)\n')
            for (i, j) in Bush.forward_star(i):
                link_weight[(i, j)] = node_weight[i] * link_likelihood[(i, j)]
                logging.warning(f'{" "*8}W_{i}{j} = W_{i} * L_{i}{j} = {node_weight[i]} * {link_likelihood[(i, j)]}'
                                f' = {link_weight[(i, j)]}')
            if len(list(Bush.forward_star(i))) == 0:
                logging.warning(f'{" "*8}Outdegrees = 0')
            node_weight[i] = sum(link_weight[(k, i)] for (k, i) in Bush.reverse_star(i))

            node_weight_calc = " + ".join([str(link_weight[(k, i)]) for (k, i) in Bush.reverse_star(i)])
            node_weight_calc_var = ' + '.join([f'W_{k}{i}' for (k, i) in Bush.reverse_star(i)])

            if len(list(Bush.reverse_star(i))) == 0:
                logging.warning(f'\n{" "*8}Indegrees = 0')
                
            if node_weight_calc:
                if len(node_weight_calc) > 1:
                    logging.warning(f'\n{" "*8}W_{i} = {node_weight_calc_var} = {node_weight_calc} = {node_weight[i]}\n')
                else:
                    logging.warning(f'\n{" "*8}W_{i} = {node_weight_calc_var} = {node_weight[i]}\n')
            else:
                logging.warning(f'\n{" "*8}W_{i} = {node_weight[i]}\n')
        logging.warning(f'Node {i} is the destinatio. So, Stopping.')
        logging.warning(f'{"-"*79}\n')
                
        # Step 3: Calculating node and link weights in reverse topological order after loading demand
        logging.warning(f'\n\nStep 3:\n{"-"*79}\n{"-"*79}')
        logging.warning(f'Reverse Topological Sort = {topsorted_node[::-1]}\n')

        link_flow = {}
        node_link_flow = {}
        for j in topsorted_node[::-1]:
            #########################################################################################
            ## Shouldn't this demand be the sum of all demand that comes to the destination node?  ##
            ## If this STOCH loading is used for every origin then sum of demand is not required.  ##
            #########################################################################################
            node_link_flow[j] = self.demand.get((origin, j), 0) + sum(link_flow.get((j, k), 0) for (j, k) in Bush.forward_star(j))

            d = str(self.demand.get((origin, j), 0))
            s = ' + '.join([str(link_flow.get((j, k), 0)) for (j, k) in Bush.forward_star(j)])
            if not s: s=0
            
            logging.warning(f'Node-{j}: \n{" "*8}X_{j} = {d} + ({s}) = {node_link_flow[j]}{" "*8}(Initializing)\n')
            
            for (i, j) in Bush.reverse_star(j):
                link_flow[(i, j)] = round(node_link_flow[j] * (link_weight[(i, j)] / node_weight[j]), 6)
                logging.warning(f'{" "*8}x_{i}{j} = {node_link_flow[j]} * ({link_weight[(i, j)]} / {node_weight[j]})'
                                f' = {link_flow[(i, j)]}\n')
        logging.warning(f'Node {j} is the origin. So, Stopping.')

        logging.warning(f'{"-"*79}\nSTOCH Finished. \n\nLink flows = \n{pformat(link_flow, indent=4)}')

        

def create_network(links=None, lp_funcs=None):
    G = Network()
    if not links:
        links = lp_funcs.keys()
    G.add_links(links)
    if lp_funcs:
        G.add_lpf(lp_funcs)
    G.topological_sort()
    G.TON = list(G.ordered_nodes.values())
    return G


if __name__ == "__main__":
    def net():
        lp_funcs = {(1, 2): 2, (1, 3): 1, (2, 3): 2, (3, 2): 2, (2, 4): 1, (3, 4): 3}
        links = list(lp_funcs.keys())
        G = create_network(links, lp_funcs)
        demand = {(1, 4): 100}
        return G, demand

    G, demand = net()

    Tap = Tap_SUE(Graph=G, demand=demand)
    Tap.load_STOCH(origin=1)

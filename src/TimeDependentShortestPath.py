#! python3

"""
Time Dependent Shortest Paths, TDSP

Algorithms
------------
A1: Fifo must hold, cost = travel time. Time can be continuous or discrete.

B: FIFO does not need to hold. Cost can be anything. Time must be discrete. 
	
	B1. Single origin, single departure times --> All destination.
	B2. Single origin, all departure times --> A single destination. (NOT IMPLEMENTED FOR NOW)
"""

__author__ = "simanta barman"
__email__ = "barma017@umn.edu"


import logging
import numpy as np
from pprint import pprint, pformat
from math import inf
from Network import Network
from network_reader import data_reader

logging.disable(logging.CRITICAL)
logging.basicConfig(level=logging.DEBUG, format='%(message)s')


class callable_dict(dict):
    def __call__(self, key):
        return self.get(key, self[len(self) - 1])
        

class TimeDependentShortestPath:

    def __init__(self, Graph, origin_node, departure_time, time_horizon, delta_t=1):
        self.G = Graph
        self.r = origin_node
        self.t0 = departure_time

        self.T = time_horizon + delta_t
        self.delta_t = delta_t
        self.__fifo = ""

    @property
    def fifo_holds(self):
        """Returns a boolean indicating whether network is FIFO or not"""
        if self.__fifo == "FIFO":
            return True
        elif self.__fifo == "NO FIFO":
            return False
        
        for t in np.arange(self.t0, self.T, self.delta_t):
            for (i, j), tau in self.G.link_cost.items():
                if tau(t + self.delta_t) - tau(t) < -self.delta_t:
                    self.__fifo = "NO FIFO"
                    return False
            else:
                self.__fifo = "FIFO"
        return True
                
    def tdsp_fifo(self):
        """Returns the time dependent shortest path labels and backnodes only for FIFO networks."""

        assert self.fifo_holds, "Network is not a FIFO network. Use tdsp_non_fifo"
        
        # Step 1: Intiialiazing the labels
        label = {i: inf for i in self.G.N}
        label[self.r] = self.t0

        # Step 2: Initializing finalized nodes and backnode vector
        finalized = set()
        backnode = {i: -1 for i in self.G.N}

        while True:
            # Step 3: Choosing an unfinalized node with the lowest label
            _, node = min([(label_val, node) for node, label_val in label.items() if node not in finalized])

            # Step 4: Finalizing selected node
            finalized = finalized.union({node}) 

            # Step 5: label and backnode vector updating
            for i, j in self.G.forward_star(node):
                old_label_j = label[j]
                temp_new_label_j = label[i] + self.G.link_cost[i, j](label[i])
                
                if temp_new_label_j > self.T:
                    continue

                # Step 5a: Updating label
                label[j] = min(temp_new_label_j, old_label_j)

                # Step 5b: Updating backnode vector
                if old_label_j != label[j]:
                    backnode[j] = i

            # Step 6: Terminating Criteria
            if len(finalized) == len(self.G.N):
                return (label, backnode)

    def tdsp_non_fifo(self):
        """Returns the time dependent shortest path labels and backnodes for any discrete time network."""

        # Step 1: Initialize backnode vector and labels
        backnode = {i: {t: -1 for t in np.arange(self.t0, self.T, self.delta_t)} for i in self.G.N}
        label = {i: {t: inf for t in np.arange(self.t0, self.T, self.delta_t)} for i in self.G.N}
        logging.debug(f'q = {backnode}')
        

        # Step 2: Set label of L_r^{t0} to be zero
        label[self.r][self.t0] = 0
        logging.debug(f'L = {label}')
        
        # Step 3: Initialize the current time to departure time
        t = self.t0
        logging.debug(f't = {t}')

        # Step 4: Update label and backnode vector
        while t < self.T: 
            for i in self.G.N:
                logging.debug(f'Node = {i} at time = {t}')
                if label[i][t] >= inf:
                    continue                
                for (i, j) in self.G.forward_star(i):
                    for t_prime in np.arange(t + self.delta_t, self.T, self.delta_t):
                        if label[j][t_prime] <= label[i][t] + self.G.tau[i, j](t) and t_prime != label[i][t] + self.G.tau[i, j](t):
                            continue
                        logging.debug(f'Link({i}:{t}, {j}:{t_prime})')
                        
                        logging.critical(f"link ({i}:{t}, {j}:{t_prime})")
                        old_label = label[j][t_prime]
                        label[j][t_prime] = min(label[i][t] + self.G.link_cost[i, j](t), label[j][t_prime])
                        logging.debug(f"Old Label = {old_label}")
                        logging.debug(f'New label = min(L[{i}][{t}] + tau[{i}, {j}]({t}), L[{j}][{t_prime}])'
                              f'          = min({label[i][t]} + {self.G.link_cost[i, j](t)}, {label[j][t_prime]})')

                        if old_label != label[j][t_prime]:
                            backnode[j][t_prime] = (i, t)
                            logging.debug(f'Label changed, q[{j}:{t_prime}] = {i}')

            t += self.delta_t
            logging.debug(f't = {t}')
        return (label, backnode)

    def shortest_path(self, end_time=None):
        """TODO: Not working."""
        
        if end_time is None:   
            _, backnodes = self.tdsp_fifo()
        else:
            _, backnodes = self.tdsp_non_fifo()
            backnodes = {node: backnodes[end_time] for node, backnodes in backnodes.items()}
            
        path = {}
        for s in self.G.N:
            node = s
            path[s] = [node]
            while node != self.r:
                node = backnodes[node]
                if node == -1:
                    break
                path[s].append(node)
        path = {s: pp[::-1] for s, pp in path.items()}
        return path

        
if __name__ == "__main__":

    from problem_nets import *
        
    Tdsp = TimeDependentShortestPath(Graph=e2(), origin_node='a', departure_time=0, time_horizon=7)
    L, q = Tdsp.tdsp_non_fifo()

    print(f'Labels = {L}\nBacknodes = {q}')

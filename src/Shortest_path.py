#! python 3

"""Chapter 2 -- Shortest paths dijkstra and A* only"""

from math import inf
import logging
from copy import deepcopy
from collections import OrderedDict
from time import time

logging.disable(logging.CRITICAL)
logging.basicConfig(level=logging.DEBUG, format='%(message)s')


class Shortest_path:

    def __init__(self, start, labels=None, heuristic=None):
        
        self.L = labels
        self.A = {link for link in labels}
        nodes = {i for i, _ in labels}.union({j for _, j in labels})
        self.N = {node: node for node in nodes}
        self.start = start
        if heuristic == None:
            self.heuristic = {node: 0 for node in self.N}
        else:
            self.heuristic = heuristic

    def forward_star(self, node):
        """Returns the Forward star for the node"""
        return self.__star(node, link_inner_index=0)

    def reverse_star(self, node):
        """Returns the Reverse star for the node"""
        return self.__star(node, link_inner_index=1)

    def __star(self, node, link_inner_index):
        """Returns the forward star or reverse star for the node"""
        for link in self.A:
            if self.N[node] == link[link_inner_index]:
                yield link

    def dijkstra(self):
        """Returns Shortest path Labels and Backnode vector for Cyclic network 
        using Dijksta's Algorithm (label setting method)"""
        return self.a_star('dijkstra')

    def a_star(self, method='a_star'):
        """Returns Shortest path Labels and Backnode vector for Cyclic network 
        using A* Algorithm (One origin to one destination)"""
        
        if method=='a_star':
            g_is = self.heuristic
            for g in g_is.values():
                if g < 0:
                    raise Exception('Heuristic values must be non-negative')
                elif g > 0:
                    logging.debug('Using A* Algorithm: ')
                    break
            else:
                logging.debug('Can\'t use A* (Invalid heuristics). Using Dijksta\'s Algorithm: ')
                g_is = {node: 0 for node in self.N}

        elif method=='dijkstra':
            logging.debug('Using Dijksta\'s Algorithm: ')
            g_is = {node: 0 for node in self.N}

        # Step 1, Initializing Label vector
        L_r = {i: inf for i in self.N}
        L_r[self.start] = 0

        # Step 2, initializing Finalized set and Backnode vector
        F = set()
        q_r = {i: -1 for i in self.N}

        iterations = 0
        while True:
            # Step 3, Selecting node based on min L_ir + g_is if it's not in F
            unfinalized = {node: (L_ir + g_is[node]) for node, L_ir in L_r.items()}
            
            for node in F:
                del unfinalized[node]
            
            min_value, i = min(zip(unfinalized.values(), unfinalized.keys()))

            # Step 4, Finalize node j by adding it to F
            F.add(i)
            if len(F) == len(self.N):
                break

            # Step 5, Update the labels for the outgoing links from node j
            for (i, j) in self.forward_star(i):
                c_ij = self.L[(i, j)]
                L_ir = min(L_r[j], (L_r[i] + c_ij))

                if L_ir < L_r[j]:
                    L_r[j] = L_ir

                    # Step 6
                    q_r[j] = i
                    
            iterations += 1
            
            # Step 7
            if len(F) == len(self.N):
                break
            
        logging.debug(f'Total iterations - {iterations}')

        return (L_r, q_r)


def main():

    # Acyclic Network
    A = [[1, 2], [1, 3], [2, 3], [2, 4], [3, 4]]
    L = {(1, 2): 2, (1, 3): 4, (2, 3): 1, (2, 4): 5, (3, 4): 2}
    heuristic = {1: 3, 2: 2, 3: 1, 4: 0}
    
    As = [['1', '2'], ['1', '3'], ['2', '3'], ['2', '4'], ['3', '4']]
    Ls = {('1', '2'): 2, ('1', '3'): 4, ('2', '3'): 1, ('2', '4'): 5, ('3', '4'): 2}
    heuristic_s = {'1': 3, '2': 2, '3': 1, '4': 0}
    
    SP = Shortest_path_dijkstra(start='1', labels=Ls, heuristic=heuristic_s)

    for algorithm in [SP.dijkstra, SP.a_star]:
        L, q = algorithm()
        print(f'L_r: {L}\nq_r: {q}\n')

    return SP

if __name__ == '__main__':
    main()

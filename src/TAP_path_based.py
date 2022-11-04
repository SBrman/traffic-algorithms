#! python3

__author__ = "Simanta Barman"
__email__ = "barma017@umn.edu"

from Network import *
from Shortest_path import *
from GetNetworkData import *
import logging
import numpy as np
import inspect
from time import time
from collections import deque

#logging.disable(logging.DEBUG)
##logging.disable(logging.ERROR)
logging.disable(logging.CRITICAL)
logging.basicConfig(level=logging.DEBUG, format='%(message)s')

def create_network(links, lp_funcs=None, alpha_labels=None):
    G = Network()
    G.add_links(links)
    if lp_funcs:
        G.add_lpf(lp_funcs)
    if alpha_labels: 
        G.add_alpha_labels(alpha_labels)
    G.topological_sort()
    G.TON = list(G.ordered_nodes.values())
    return G

def ddx(f, x):
    h = 1e-6
    return (f(x+h) - f(x-h)) / (2*h)


class TAP:
    def __init__(self, Graph, demand, routes=None):
        # Getting routes from demand matrix if not given.
        if routes == None:
            routes = (tuple({r for r, _ in demand.keys()}), \
                        tuple({s for _, s in demand.keys()}))
        self.G = Graph
        self.demand = demand
        self.routes = routes
        self.origin = routes[0]
        self.destination = routes[1]

    def __get_x_dict(self, h_dict):
        x_dict = {link: 0 for link in self.G.A}
        for path, path_flow in h_dict.items():
            if path == None:
                continue
            for link in path:
                x_dict[link] += path_flow
        return x_dict

    def get_gamma_and_aec(self, t, x, k, d_rs):
        """Returns the convergence measures, Relative Gap (γ1) and Average
        Excess Cost (AEC)"""
        t_dot_x = sum((t[link]*x[link] for link in t))
        k_dot_d = sum((k.get(path, 0) * d_rs.get(path, 0) for path in k))
        gamma1 = ((t_dot_x / k_dot_d) - 1)
        aec = (t_dot_x - k_dot_d) / sum(self.demand.values())
        
        print(f'Rel Gap = {gamma1}\t\t\t AEC = {aec}')
        return gamma1, aec

    def __get_A3(self, longest_paths):
        """Generator for links in non-basic path only"""
        for path in longest_paths:
            for link in path:
                yield link

    def __get_A4(self, shortest_paths):
        """Generator for links in basic path only"""
        for path in shortest_paths:
            for link in path:
                yield link

    def __get_ddx_tij_for_links_in_A3_U_A4(self, longest_paths, shortest_paths, h_dict):
        x_dict = self.__get_x_dict(h_dict)
        A3_U_A4 = set(self.__get_A3(longest_paths)).union(set(self.__get_A4(shortest_paths)))

        h = 1e-6
        ddx = lambda f, x: ((f(x + h) - f(x - h)) / (2*h))

        ddx_tij_wrt_xij = sum((ddx(self.G.lp_func[link], x_dict[link]) for link in A3_U_A4))

        return np.round(ddx_tij_wrt_xij, 6)
    
    def __get_SPTT(self, start, end, t_ij=None):
        """Returns the shortest path (all the links in the path) and shortest paths travel time."""

        # Getting corrected_t_ij which contains only the labels for links that are only in path from start to end.
        # Otherwise links in other paths are included in the shortest path calculation
        label_nodes = set()
        queue = [start]
        links = set(t_ij)
        while queue:
            node = queue.pop(0)
            for i, j in links:
                if node == i:
                    if node not in queue and node not in label_nodes:
                        queue.append(j)
            label_nodes.add(node)
        corrected_t_ij = {(i, j): label for (i, j), label in t_ij.items() if i in label_nodes and j in label_nodes}

        P = Shortest_path(start=start, labels=corrected_t_ij)
        Lr, backnode = P.dijkstra()
        path = [end]
        while True:
            if path[-1] == start:
                break
            elif path[-1] in backnode:
                path.append(backnode[path[-1]])
            else:
                path = None
                break

        if path:
            shortest_path = path[::-1]
            shortest_path_links = tuple(zip(shortest_path[:-1], shortest_path[1:]))
            shortest_path_time = sum((t_ij[(i, j)] for i, j in shortest_path_links))
            return shortest_path_links, shortest_path_time
        else:
            return None, None

    def __get_links_from_paths(self, path, path_type):
        """Generator for all the links from all the paths"""
        for link in path:
            yield link
            
    def __get_A3(self, longest_path):
        """Generator for links in non-basic path only"""
        return self.__get_links_from_paths(longest_path, 'lp')

    def __get_A4(self, shortest_path):
        """Generator for links in basic path only"""
        return self.__get_links_from_paths(shortest_path, 'sp')
    
    def __get_delta_h(self, sps, lps, c_pi_star, c_pi):
        """Returns delta_h"""
        
        c_pi_delta = np.round((c_pi - c_pi_star), 6)
        
        x_dict = self.__get_x_dict(self.h_dict)
        A3_U_A4 = set(self.__get_A3(lps)).union(set(self.__get_A4(sps)))
        A3_cap_A4 = set(self.__get_A3(lps)).intersection(set(self.__get_A4(sps)))

        ddx_tij_wrt_xij = sum((ddx(self.G.lp_func[link], x_dict[link]) for link in A3_U_A4 if link not in A3_cap_A4))
        print(f'ddx = {list((ddx(self.G.lp_func[link], x_dict[link]), link) for link in A3_U_A4 if link not in A3_cap_A4)}')
        delta_h = np.round(c_pi_delta / ddx_tij_wrt_xij, 1)

        return delta_h

    def __get_path_travel_time(self, path):
        time = 0
        for link in path:
            time += self.t_dict[link]
        return time

    def __shift_path_flows(self, r, s, basic_path):
        if len(self.used_paths[(r, s)]) == 1:
            self.h_dict[list(self.used_paths[(r, s)])[0]] = self.demand[(r, s)]
        else:
            c_pi_star = self.__get_path_travel_time(basic_path)
            print(f"c_pi_star {basic_path} = {c_pi_star}")
            non_basic_paths = {path for path in self.used_paths[(r, s)] if path != basic_path}
            for non_basic_path in non_basic_paths:
                c_pi = self.__get_path_travel_time(non_basic_path)
                print(f"c_pi {non_basic_path} = {c_pi}")
                delta_h = self.__get_delta_h(basic_path, non_basic_path, c_pi_star, c_pi)
                print(f'delta_h = {delta_h}')
                self.h_dict[non_basic_path] -= delta_h
            self.h_dict[basic_path] = self.demand[(r, s)] - sum(self.h_dict[non_basic_path] for non_basic_path in non_basic_paths)

    def gradient_projection(self, GAMMA_1=0.02, AEC=1e-5, max_iter=30):
        
        self.h_dict = {path: 0 for (r, s) in self.demand for path in self.G.find_all_paths(start=r, end=s)[0]}
        x_dict = self.__get_x_dict(self.h_dict)
        self.t_dict = {(i, j): lpf(0) for (i, j), lpf in self.G.lp_func.items()}     # freeflow travel times
        kappa = {}

        # 1. Initializing used_paths
        self.used_paths = {(r, s): set() for r in self.origin for s in self.destination if self.demand[(r,s)] > 0}
        print(len(self.used_paths))        
        
        # 2.
        i = 0
        while i <= max_iter:
            i += 1
            print(f'Iteration - {i}\n\nUsed_paths = {self.used_paths}')
                
            for (r, s), used_path in self.used_paths.items():
                print(f'Inner iter for {(r, s)}')
                print(f'\n\nUsed path {(r, s)} = {used_path}:\n')
                print(x_dict, self.t_dict)
                
                # 2.a Shortest path
                shortest_path, shortest_path_time = self.__get_SPTT(start=r, end=s, t_ij=self.t_dict)
                if shortest_path == None:
                    continue
                
                kappa[(r, s)] = shortest_path_time
                
                # Adding shortest path to the set of used paths
                self.used_paths[(r, s)].add(shortest_path)
                print(f'\nAdding SP {shortest_path} to used_paths: \nnew UP = {self.used_paths}\n\n')

                # 2.b Shift flow
                self.__shift_path_flows(r, s, shortest_path)
                
                # 2.c update travel times
                x_dict = self.__get_x_dict(self.h_dict)
                self.t_dict = {(i, j): lpf(x_dict[(i, j)]) for (i, j), lpf in G.lp_func.items()}

                # 3. drop unused paths
                for (r, s), paths in self.used_paths.copy().items():
                    for path in paths:
                        if not self.h_dict[path] != 0:
                            self.used_paths[(r, s)].discard(path)

            kappa = {}
            for (r, s), used_path in self.used_paths.items():
                shortest_path, shortest_path_time = self.__get_SPTT(start=r, end=s, t_ij=self.t_dict)
                kappa[(r, s)] = shortest_path_time if shortest_path_time < kappa.get((r, s), 100000) else 100000
                            
            print(f'Path flows = {self.h_dict}\nLink flows = {x_dict}\nTravel_times = {self.t_dict}\nkappa = {kappa}\nd = {self.demand}\n')

            gamma_1, aec = self.get_gamma_and_aec(t=self.t_dict, x=x_dict, k=kappa, d_rs=self.demand)

            if gamma_1 <= GAMMA_1 or aec <= AEC or i >= max_iter:
                return x_dict, self.t_dict

            
        
if __name__ == "__main__":
    from problem_nets import *

    # Get the Network and OD matrix
    G, d, routes = neww()
    
    Tap = TAP(Graph=G, demand=d, routes=None)
    x = Tap.gradient_projection()
    print(x)

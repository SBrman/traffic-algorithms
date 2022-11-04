#! python3

from Network import *
from Shortest_path import *
import logging
from pprint import pprint

logging.disable(logging.INFO)
logging.basicConfig(level=logging.DEBUG, format="%(message)s")


class TAP:
    def __init__(self, Graph, demand):
        self.G = Graph
        self.demand = demand
        self.origin, self.destination =  zip(*list(demand.keys()))

    def __t(self, x):
        return {link: tij(x) for link, tij in self.G.lp_funcs_vect_vals.items()}

    def __get_t_dot_x(self, t, x):
        # t here is a dictionary not functions
        return sum(t[link]*x[link] for link in self.G.A)
        
    def _restricted_aec(self, x, fancy_x):
        """Returns restricted Average access cost.
        x = current solution, dict(links, linflows)
        fancy_x = set of target solutions {dict(links, linflows), ...}"""
        
        t = self.__t(x)
        min_i_tx = min(self.__get_t_dot_x(t, x_star) for x_star in fancy_x)
        
        return (self.__get_t_dot_x(t, x) - min_i_tx) / sum(self.demand.values())
    
    def _get_SPTT(self, start, end, t_ij=None):
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

    def get_kappa_and_sps(self, t_vector_dict):
        """Returns kappa and shortest paths,
        kappa = {(r, s): path_travel_time}
        shortest_paths = {(r, s): [links_in_shortest_path_from_r_to_s]}
        """

        kappa, shortest_paths = {}, {}
        for r in self.origin:
            for s in self.destination:
                if self.demand[(r, s)] == 0:
                    continue
                shortest_path, shortest_path_time = self._get_SPTT(start=r, end=s, t_ij=t_vector_dict)
                if not shortest_path:
                    continue
                shortest_paths[shortest_path] = shortest_path_time
                kappa[(r, s)] = shortest_path_time
##        logging.info(f'shortest_paths = {shortest_paths}')
        return kappa, shortest_paths
        
    def __get_all_or_nothing_x(self, shortest_paths):
        """Returns all or nothing assignments link flow"""
        x_aon = {link: 0 for link in self.G.A}
        for path in shortest_paths:
            for link in path:
                x_aon[link] += self.demand[(path[0][0], path[-1][-1])]
        return x_aon

    def __get_delta_x(self, t, x, fancy_x):
        # numerator_first_term
        nft = []
        for x_star in fancy_x:
            x_minus_x_star = {link: x[link] - x_star[link] for link in x}
            nft.append(max(self.__get_t_dot_x(t, x_minus_x_star), 0))
            logging.debug(nft)
   
        if not (denominator := sum(nft)):
            return {link: 0 for link in self.G.A}

        delta_x = {}
        for i, x_star in enumerate(fancy_x):
            
            for link in self.G.A:
                delta_x[link] = delta_x.get(link, 0) + (nft[i] / denominator) * (x_star[link] - x[link])

        return delta_x

            
    def _shift_flows(self, x, mu, delta_x):
        return {link: x[link] + mu * delta_x[link] for link in self.G.A}
    
    def simplical_decomposition(self, AEC_LIMIT=1e3):

        # Step 1: Initialize fancy_x
        fancy_x = []

        # Step 2: Find shortest paths for all OD pairs
        # Getting shortes_paths for all od pairs using freeflow travel times
        x = {link: 0 for link in self.G.A}
        t = {link: t_ij(x) for link, t_ij in self.G.lp_funcs_vect_vals.items()}

        i = 1
        while True:
            logging.warning(f"\n\nIteration-{i}:\nTravel times = {t}")
            logging.warning(f"x = {x}")
            
            _, shortest_paths = self.get_kappa_and_sps(t) # {(r, s): path}
            logging.warning(f"Shortest path = {shortest_paths}")

            # Step 3: Form all or nothing assignment x* based on the shortest paths
            x_star = self.__get_all_or_nothing_x(shortest_paths)
            logging.warning(f"x* = {x_star}")

            # Step 4: if x* in fancy_x stop (Terminating criteria)
            if x_star in fancy_x:
                break
            
            # Step 5: add x* to fancy_x
            fancy_x.append(x_star)
            logging.warning(f"Fancy x = {fancy_x}")

            # Step 6: Subproblem: Find a restricted equilibrium x using only the vectors in fancy_x
            # Start of Subproblem:
            j = 1
            while True:
                # (a) Find the improvement direction delta_x using equation (9.34).
                delta_x = self.__get_delta_x(t, x, fancy_x)
                logging.info(delta_x)
                
                # (b) Update x = x + mu * delta_x with mu suffiencienty small (to reduce AEC')
                mu = 1/j
                logging.info(f"mu = {mu}")

                if i == 1:
                    x = x_star

                x = self._shift_flows(x, mu, delta_x) # x + mu * delta_x
                logging.info(f"New x = {x}")
                
                # (c) Update travel times
                t = {link: t_ij(x) for link, t_ij in self.G.lp_funcs_vect_vals.items()}
                logging.info(f"New t = {t}")
                
                # (d) Return to step 1 of subproblem unless AEC' is small enough
                if (aec := self._restricted_aec(x, fancy_x)) < AEC_LIMIT:
                    break
                
                j += 1

            # Return to step 2, (while loop takes care of that)
            i += 1

        return x
            
if __name__ == '__main__':

    def net():
        # Defining the network
        G = Network()

        # Following dict must have this format, (all the key value pairs must be in new and separate lines)
        # Dict = {
        #           key: value pairs
        #        }
        G.lp_funcs_vect_vals = {
                                (1, 2): lambda x: 10 * x.get((1, 2), 0),
                                (1, 3): lambda x: 50 + x.get((1, 3), 0) + 0.5 * x.get((2, 3), 0),
                                (2, 3): lambda x: 10 + x.get((2, 3), 0) + 0.5 * x.get((1, 3), 0),
                                (2, 4): lambda x: 50 + x.get((2, 4), 0) + 0.5 * x.get((3, 4), 0),
                                (3, 4): lambda x: 10 * x.get((3, 4), 0) + 5 * x.get((2, 4), 0)
                                }
        
        G.add_links(G.lp_funcs_vect_vals.keys())
        
        # Defining demand for origin to destination (OD matrix)
        d = {(1, 4): 6}
        routes = None

        return G, d, routes
    

    def sheffi():
        # Defining the network
        G = Network()

        G.lp_funcs_vect_vals = {
                                 (1, 2): lambda x: 1.5 + 0.01 * x.get((1, 2), 0) **2 ,
                                 (1, 4): lambda x: 2.5 + 0.03 * x.get((1, 4), 0) **2 ,
                                 (2, 3): lambda x: 2.5 + 0.01 * x.get((2, 3), 0) **2 ,
                                 (2, 5): lambda x: 3.5 + 0.03 * x.get((2, 5), 0) **2 + 0.015 * x.get((4, 5), 0) ** 2,
                                 (3, 6): lambda x: 4.5 + 0.03 * x.get((3, 6), 0) **2 + 0.015 * x.get((5, 6), 0) ** 2,
                                 (4, 5): lambda x: 4.5 + 0.01 * x.get((4, 5), 0) **2 + 0.0025 * x.get((2, 5), 0) ** 2,
                                 (4, 7): lambda x: 5.5 + 0.03 * x.get((4, 7), 0) **2 ,
                                 (5, 6): lambda x: 5.5 + 0.01 * x.get((5, 6), 0) **2 + 0.0025 * x.get((3, 6), 0) ** 2,
                                 (5, 8): lambda x: 6.5 + 0.03 * x.get((5, 8), 0) **2 + 0.015 * x.get((7, 8), 0) ** 2,
                                 (6, 9): lambda x: 7.5 + 0.03 * x.get((6, 9), 0) **2 + 0.015 * x.get((8, 9), 0) ** 2,
                                 (7, 8): lambda x: 7.5 + 0.01 * x.get((7, 8), 0) **2 + 0.0025 * x.get((5, 8), 0) ** 2,
                                 (8, 9): lambda x: 8.5 + 0.01 * x.get((8, 9), 0) **2 + 0.0025 * x.get((6, 9), 0) ** 2
                                }
        
        G.add_links(G.lp_funcs_vect_vals.keys())
        
        # Defining demand for origin to destination (OD matrix)
        d = {(1, 9): 12, (4, 9): 4}
        routes = None

        return G, d, routes
    
    # Selecting a network
    G, d, routes = sheffi()
    
    # Initializing the Traffic assignment problem
    Tap = TAP(Graph=G, demand=d)
    x = Tap.simplical_decomposition()
    print(x)

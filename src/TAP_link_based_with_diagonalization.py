#! python3

__author__ = "Simanta Barman"
__email__ = "barma017@umn.edu"

from Network import *
from GetNetworkData import *
from Shortest_path import *
import logging
from pprint import pformat
import numpy as np
import sympy
import time
from Root_finder import Root_finder
import inspect
import re

logging.disable(logging.WARNING)
logging.disable(logging.CRITICAL)
logging.basicConfig(level=logging.DEBUG, format='%(message)s')


class TAP:
    def __init__(self, Graph, demand, routes=None):
        self.G = Graph
        self.demand = demand
        self.origin = tuple({r for r, _ in demand.keys()})
        self.destination = tuple({s for _, s in demand.keys()})

    def get_gamma_and_aec(self, t, x, k, d_rs):
        """Returns the convergence measures, Relative Gap (γ1) and Average
        Excess Cost (AEC)"""
        t_dot_x = sum(t[link]*x[link] for link in t)
        k_dot_d = sum(k.get(path, 0) * d_rs.get(path, 0) for path in d_rs)
##        lpf = lambda i, j: inspect.getsource(self.G.lp_func[(i,j)]).split("return")[1].strip().strip('(').strip(')')
##        beckmann = sum(sympy.integrate(lpf(i, j), (Symbol("x"), 0, x[(i, j)])) for (i, j) in self.G.A)
        
        gamma1 = np.round(np.round(t_dot_x / k_dot_d, 6) - 1, 6)
        aec = (t_dot_x - k_dot_d) / sum(self.demand.values())
        print(f'Rel Gap = {gamma1}\t\t AEC = {aec}')
        return gamma1, aec

    def msa(self, gamma1_limit=1e-4, aec_limit=0.02):
        """Returns the final x, t, iteration, γ1, AEC using the
        Method of Successive Averages (MSA)"""
        return self.__tap_algo(gamma1_limit, aec_limit, 'MSA')
    
    def frank_wolfe(self, gamma1_limit=1e-4, aec_limit=0.01, root_finder='Newtons method'):
        """Returns the final x, t, iteration, γ1, AEC using the
        Frank Wolfe (FW) Algorithm"""
        return self.__tap_algo(gamma1_limit, aec_limit, 'FW', root_finder)

    def __adaptive_step_size(self, x_dict, x_star_dict, root_finder):
        """Returns adaptive step size, λ by minimizing the ζ(λ) function."""
        
        ζ_prime_λ_str = ''
        x, λ = sympy.symbols('x λ')
        for (i, j), lp_func in self.G.lp_func.items():
            new_x = λ*x_star_dict[(i,j)] + (1-λ)*x_dict[(i,j)]
            ζ_prime_λ_str += ' + ' + str(lp_func(new_x) * (x_star_dict[(i,j)] - x_dict[(i,j)]))
        logging.critical(f'ζ\'(λ) = {ζ_prime_λ_str}\n')
        ζ_prime_λ_str = sympy.simplify(ζ_prime_λ_str)      # simplifies the equation but this is very slow
        logging.critical(f'       = {ζ_prime_λ_str}\n')
        ζ_prime_λ = sympy.lambdify(λ, ζ_prime_λ_str)
        
        # Getting the root the zeta function using bisection or newton's method
        ζ_prime_roots = Root_finder(ζ_prime_λ)
        λ = ζ_prime_roots.bisect() if root_finder == 'Bisection' else ζ_prime_roots.newton_raphson()
        logging.critical(f'λ = {λ}')
        return λ

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
        logging.warning(f'shortest_paths = {shortest_paths}')
        return kappa, shortest_paths
        
    def __get_all_or_nothing_x(self, shortest_paths):
        """Returns all or nothing assignments link flow"""
        x_aon = {link: 0 for link in self.G.A}
        for path in shortest_paths:
            for link in path:
                x_aon[link] += self.demand[(path[0][0], path[-1][-1])]

        return x_aon

    def __diagonalize(self, link_flows):
        t = {}
        for link, lpf in self.G.lp_funcs_vect_vals.items():
            lp_func_str_source = inspect.getsourcelines(G.lp_funcs_vect_vals[link])[0][0].strip().split('lambda x: ')[1].split('+')
            lp_func_str = ''
            for ii in lp_func_str_source:
                try:
                    if float(ii):
                        lp_func_str += ' + ' + str(ii)
                except:
                    pass
                if str(link) in ii:
                    lp_func_str += ' + ' + ii.replace(f'x.get({link}, 0)', 'x').strip(',')
                elif str(link) not in ii:
                    num, iii = m if len(m := ii.split('*')) == 2 else (1, m[0])
                    if other_link := re.findall(r'\((\(.+?, .+?\))', ii):
                        lp_func_str += f' + {num} * ' + str(link_flows[literal_eval(other_link[0])])
                
                t[link] = lp_func_str
##        print(f'\nDiagonalized link performance functions: \n{pformat(t)}\n')
        self.G.add_lpf(t)
        
    def __tap_algo(self, gamma1_limit, aec_limit, method, root_finder='Newtons method', epsilon=0.01):

        # Freeflow travel times
        t = {link: lp_func({link: 0 for link in self.G.A})
             for link, lp_func in self.G.lp_funcs_vect_vals.items()}
        
        i = 1
        while True:

            print(f'Iteration - {i}:\t', end=' ')
            # Step 1: Find the shortest path between each origin and destination.
            # Then check termination criteria if iteration is not 1.
            kappa, shortest_paths = self.get_kappa_and_sps(t)
            
            if i != 1:
                # Getting Relative gap and AEC
                gamma1, aec = self.get_gamma_and_aec(t, x, kappa, self.demand)
                # Terminating Criteria
                if abs(gamma1) == gamma1_limit or abs(aec) <= aec_limit or i > 30:
                    x = {link: np.round(flow, 6) for link, flow in x.items()}
                    break
                
            # Step 2. Shift travelers onto shortest paths.
            # 2a. Find the link flows if everybody were traveling on the
            # shortest paths found in step 1, store these in x*.
            x_star = self.__get_all_or_nothing_x(shortest_paths)

            # 2b. If this is the first iteration, set x = x* and move to step 3. Otherwise,
            # continue with step c.
            if i == 1:
                x = x_star
            else:
                # 2c. Using the current solution x, form the diagonalized link performance
                # functions for each links
                diagonalized_lp_funcs = self.__diagonalize(x)

                # 2d Find lambda
                _lambda = self.__adaptive_step_size(x, x_star, "Newton-Raphson")

                # Update x
                x = {link: _lambda*x_star[link] + (1-_lambda)*x[link] for link in self.G.A}

                
            t = {link: tt(x) for link, tt in self.G.lp_funcs_vect_vals.items()}

##            print(f'\nx = {pformat(x)}\n')
##            print(f'x_star = {pformat(x_star)}\n')
##            print(f't = {pformat(t)}')
            
            i += 1

        return x, t, i, gamma1, aec

def print_result(start, x, t, i, gamma, aec):
    print(f'Total Iterations = {i}\t(Took {time.time() - start} seconds)\nγ1 = {gamma}\tAEC = {aec}\nx⃗ = '
          f'\n{pformat(x, indent=8)}\nTravel times for the links:\nt = \n{pformat(t, indent=8)}')
    

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

        # Following dict must have this format, (all the key value pairs must be in new and separate lines)
        # Dict = {
        #           key: value pairs
        #        }
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
        
        
    # Parameters for determining terminating criteria
    LIMITS = {'gamma1_limit': 0.001, 'aec_limit': 0.01}
    
    # Selecting a network
    G, d, routes = sheffi()
    
    # Initializing the Traffic assignment problem
    Tap = TAP(Graph=G, demand=d, routes=None)

    # Method of Successive Averages
    start = time.time()
    x, t, i, gamma, aec = Tap.msa(**LIMITS)
    print_result(start, x, t, i, gamma, aec)
    
    # Frank Wolfe (Bisection)
    start = time.time()
    x, t, i, gamma, aec = Tap.frank_wolfe(**LIMITS, root_finder='Bisection')
    print_result(start, x, t, i, gamma, aec)
    
    # Frank Wolfe (Newtons Method)
    start = time.time()
    x, t, i, gamma, aec = Tap.frank_wolfe(**LIMITS, root_finder='Newtons method')
    print_result(start, x, t, i, gamma, aec)


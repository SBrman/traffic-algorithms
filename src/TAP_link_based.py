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

#logging.disable(logging.NOTSET)
logging.disable(logging.DEBUG)
#logging.disable(logging.INFO)
#logging.disable(logging.WARNING)
logging.disable(logging.CRITICAL)
logging.basicConfig(level=logging.DEBUG, format='%(message)s')


class TAP:
    def __init__(self, Graph, demand, routes=None):

        # Getting routes from demand matrix if not given.
        if routes == None:
            routes = (tuple({r for r, _ in demand.keys()}), tuple({s for _, s in demand.keys()}))
        self.G = Graph
        self.demand = demand
        self.routes = routes

    def get_gamma_and_aec(self, t, x, k, d_rs):
        """Returns the convergence measures, Relative Gap (γ1) and Average
        Excess Cost (AEC)"""
        t_dot_x = sum((t[link]*x[link] for link in t))
        k_dot_d = sum((k.get(path, 0) * d_rs.get(path, 0) for path in d_rs))
        lpf = lambda i, j: inspect.getsource(self.G.lp_func[(i,j)]).split("return")[1].strip().strip('(').strip(')')
        beckmann = sum(sympy.integrate(lpf(i, j), (Symbol("x"), 0, x[(i, j)])) for (i, j) in self.G.A)
        print(f'TSTT = {t_dot_x}\t\tSPTT = {k_dot_d}\t\tBeckmann objective value = {beckmann}')
        
        gamma1 = ((t_dot_x / k_dot_d) - 1)
        aec = (t_dot_x - k_dot_d) / sum(self.demand.values())
        print(f'Rel Gap = {gamma1}\t\t\t AEC = {aec}\n')
        return gamma1, aec

    def msa(self, gamma1_limit=1e-4, aec_limit=0.02):
        """Returns the final x, t, iteration, γ1, AEC using the
        Method of Successive Averages (MSA)"""
        return self.__tap_algo(gamma1_limit, aec_limit, 'MSA')
    
    def frank_wolfe(self, gamma1_limit=1e-4, aec_limit=0.01, root_finder='Newtons method'):
        """Returns the final x, t, iteration, γ1, AEC using the
        Frank Wolfe (FW) Algorithm"""
        return self.__tap_algo(gamma1_limit, aec_limit, 'FW', root_finder)

    def conjugate_frank_wolfe(self, gamma1_limit=1e-4, aec_limit=0.01, epsilon=0.01, root_finder='Newtons method'):
        """Returns the final x, t, iteration, γ1, AEC using the
        Frank Wolfe (FW) Algorithm"""
        return self.__tap_algo(gamma1_limit, aec_limit, 'CFW', root_finder, epsilon)

    def __adaptive_step_size(self, x_dict, x_star_dict, root_finder):
        """Returns adaptive step size, λ by minimizing the ζ(λ) function."""
        
        ζ_prime_λ_str = ''
        x, λ = sympy.symbols('x λ')
        for (i, j), lp_func in self.G.lp_func.items():
            new_x = λ*x_star_dict[(i,j)] + (1-λ)*x_dict[(i,j)]
            ζ_prime_λ_str += ' + ' + str(lp_func(new_x) * (x_star_dict[(i,j)] - x_dict[(i,j)]))
        #ζ_prime_λ_str = sympy.simplify(ζ_prime_λ_str)      # simplifies the equation but this is very slow
        logging.info(f'ζ\'(λ) = {ζ_prime_λ_str}\n')
        ζ_prime_λ = sympy.lambdify(λ, ζ_prime_λ_str)
        
        # Getting the root the zeta function using bisection or newton's method
        ζ_prime_roots = Root_finder(ζ_prime_λ)
        λ = ζ_prime_roots.bisect() if root_finder == 'Bisection' else ζ_prime_roots.newton_raphson()
        return λ

    def __get_alpha(self, x, x_star, x_aon, epsilon=0.01):
        """Returns the alpha value."""
        ddx = lambda f, x: ((x+1e-6)-(x-1e-6)) / (2*1e-6)
        t_prime_dict = {link: ddx(lp_func, x[link]) for link, lp_func in self.G.lp_func.items()}

        x_star_old_minus_x = {link: x_star[link] - x[link] for link in self.G.A}
        x_aon_minus_x = {link: x_aon[link] - x[link] for link in self.G.A}
        x_aon_minus_x_star = {link: x_aon[link] - x_star[link] for link in self.G.A}

        numerator = sum((x_star_old_minus_x[link]*x_aon_minus_x[link]*t_prime_dict[link] for link in self.G.A))
        denominator = sum((x_star_old_minus_x[link]*x_aon_minus_x_star[link]*t_prime_dict[link] for link in self.G.A))

        if denominator == 0:
            return 0
        else:
            alpha = (numerator / denominator)
            return min([alpha, 1-epsilon]) if alpha > 0 else 0
        
        #return max([min([(numerator / denominator), 1-epsilon]), 0]) if denominator != 0 else 0
    
    def __get_title(self, method, root_finder):
        """Returns the title of the method being used."""
        title = '\n' + ('Conjugate Frank-Wolfe' if method == 'CFW' else ('Frank-Wolfe' if method == 'FW' else 'MSA'))
        if method != 'MSA': title += f' ({root_finder})'
        title += ': \n' + '-'*79
        return title

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
        for r in self.routes[0]:
            for s in self.routes[1]:
                if self.demand[(r, s)] == 0:
                    continue
                shortest_path, shortest_path_time = self._get_SPTT(start=r, end=s, t_ij=t_vector_dict)
                if not shortest_path:
                    continue
                shortest_paths[shortest_path] = shortest_path_time
                kappa[(r, s)] = shortest_path_time
        logging.warning(f'shortest_paths = {shortest_paths}')
        return kappa, shortest_paths
        
    
    def __tap_algo(self, gamma1_limit, aec_limit, method, root_finder='Newtons method', epsilon=0.01):
        print(self.__get_title(method, root_finder))
        
        # Initialization: set i = 1
        i, λ = 1, 1

        # Iteration 1: Finding the shortest path by setting both path's link flows to zero.
        x_dict = {link: 0 for link in self.G.A}
        
        while True:
            i_text = i-1 if method == 'CFW' else i
            logging.warning(f'Iteration-{i_text}'.center(79) +'\n'+ '-'*79)

            
            logging.warning(f'x⃗ = {x_dict}')

            # Getting the t vector by substituting x in the link performance functions.
            t_vector_dict = {link: lp_func(x_dict[link]) for link, lp_func in self.G.lp_func.items()}
            logging.warning(f't⃗ = {t_vector_dict}')

            # Getting Kappa (shortest path time) values
            kappa, shortest_paths = self.get_kappa_and_sps(t_vector_dict)

            # Getting the All or nothing target x vector
            x_aon_dict = {link: 0 for link in self.G.A}
            for path in shortest_paths:
                for link in path:
                    x_aon_dict[link] += self.demand[(path[0][0], path[-1][-1])]         # (r, s) = path[0][0], path[-1][-1]
            logging.warning(f'x_aon⃗ = {x_aon_dict}')
            
            # Getting ther target x* vector
            if method != 'CFW' or i <= 1:
                x_star_dict = x_aon_dict
            else:
                # Get the α value and x*
                α = self.__get_alpha(x_dict, x_star_dict, x_aon_dict, epsilon)
                x_star_dict = {link: α*x_star_dict[link] + (1-α)*x_aon_dict[link] for link in self.G.A}

            logging.warning(f'x*⃗ = {x_star_dict}')

            if i != 1:
                # d_rs = {(r, s): self.demand[(r, s)] for path in shortest_paths if self.demand[(r := path[0][0], s := path[-1][-1])] != 0}
                d_rs = self.demand
                logging.warning('\nγ1 and AEC calculation:')
                parameters = {'t⃗': t_vector_dict, 'x⃗': x_dict, 'κʳˢ⃗': kappa, 'dʳˢ⃗': d_rs}
                for name, parameter in parameters.items():
                    logging.info(f'{name} = {parameter}')

                # Getting Relative gap and AEC
                gamma1, aec = self.get_gamma_and_aec(*parameters.values())
                logging.warning(f'\nγ1 = {gamma1}\tAEC = {aec}\n\nx⃗ = {x_dict}\n\n')

                # Terminating Criteria
                if abs(gamma1) < 0.02:#gamma1_limit or abs(aec) <= aec_limit or i == 31:
                    x_dict = {link: np.round(flow, 2) for link, flow in x_dict.items()}
                    break
                
                # Getting the step size λ
                λ = self.__adaptive_step_size(x_dict, x_star_dict, root_finder) if method in ['FW', 'CFW'] else 1/(i+1)
            logging.warning(f'λ = {λ}')

            # Updating the x vector by shifting travelers
            x_dict = {link: λ*x_star_dict[link] + (1-λ)*x_dict[link] for link in self.G.A}
            logging.warning(f'x⃗ = {x_dict}\n')

            # Without rounding MSA takes a lot of time
            if method == 'MSA':
                x_dict = {link: np.round(flow, 4) for link, flow in x_dict.items()}
                x_star_dict = {link: np.round(flow, 4) for link, flow in x_star_dict.items()}

            # Next Iteration
            i += 1

        return x_dict, t_vector_dict, i, np.round(float(gamma1), 5), np.round(float(aec), 4)

def print_result(start, x, t, i, gamma, aec):
    print(f'Total Iterations = {i}\t(Took {time.time() - start} seconds)\nγ1 = {gamma}\tAEC = {aec}\nx⃗ = '
          f'\n{pformat(x, indent=8)}\nTravel times for the links:\nt = \n{pformat(t, indent=8)}')
    

if __name__ == '__main__':

    from problem_nets import *
        
    # Parameters for determining terminating criteria
    LIMITS = {'gamma1_limit': 0.02, 'aec_limit': 0.01}
    
    # Selecting a network
    G, d, routes = neww()
    
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
    
    # Conjugate Frank Wolfe (Bisection)
    start = time.time()
    x, t, i, gamma, aec = Tap.conjugate_frank_wolfe(**LIMITS, epsilon=0.01, root_finder='Bisection')
    print_result(start, x, t, i, gamma, aec)

    # Conjugate Frank Wolfe (Newtons Method)
    start = time.time()
    x, t, i, gamma, aec = Tap.conjugate_frank_wolfe(**LIMITS, epsilon=0.01, root_finder='Newtons method')
    print_result(start, x, t, i, gamma, aec)

#! python3

__author__ = "Simanta Barman"
__email__ = "barma017@umn.edu"

import numpy as np
from math import inf
from sympy import Symbol, lambdify, Rational
from fractions import Fraction
from ast import literal_eval
import pyinputplus as pyip
import logging
from pprint import pprint, pformat
from itertools import product, combinations
from copy import deepcopy
import inspect
from collections import OrderedDict

logging.disable(logging.CRITICAL)
logging.basicConfig(level=logging.DEBUG, format='%(message)s')


class Network:
    
    def __init__(self):
        # Links
        self.A = set()
        # Nodes
        self.N = {}
        # link Performance Function
        self.lp_func = {}
        self.alphas = {}

    def __update_nodes(self):

        self.N = self.__get_nodes_from_links()

    def __get_nodes_from_links(self, links=None):
        if links == None:
            links = self.A
        # Key and value are the same
        return {node: node for link in self.A for node in link}
        
    def add_links(self, links=None):
        if links != None: 
            for link in links:
                try:
                    if len(link) == 2:
                        self.A.add(tuple(link))
                    else:
                        print('Enter a valid link')
                except TypeError:
                    print('Enter a valid link')
        else:
            while True:
                try:
                    inp = literal_eval(input('Enter link: '))
                    if type(inp) == tuple:
                        self.A.add(inp)
                except:
                    break
        self.__update_nodes()
        self.topological_sort()

    def del_links(self, links):
        for link in links:
            if link == None:
                continue     
            if link in self.A.copy():
                self.A.discard(link)
            else:
                print(f'Link {link} does not exist in the set of links in this Network.')
            if link in self.lp_func.copy():
                del self.lp_func[link]
        self.__update_nodes()

    def add_lpf(self, lp_funcs=None, link=None):
        x = Symbol('x')
        if lp_funcs != None:
            if type(lp_funcs) == dict:
                for link, func in lp_funcs.items():
                    self.A.add(link)
                    self.lp_func[link] = lambdify(x, func, 'numpy') 
            elif type(lp_funcs) == str and link != None:
                self.A.add(link)
                self.lp_func[link] = lambdify(x, lp_funcs, 'numpy')
            elif link == None:
                    print('Did not specify a link')

        else:
            inp = pyip.inputYesNo('Manually enter Link Performance Functions for each link \n(Enter Yes/No): ')
            if inp != 'no':
                links = list(self.A)
                links.sort()
                for link in links:
                    while True:
                        func_inp = input(f'Enter the Link Performance Function for link {link}: ')
                        if func_inp.isprintable():
                            break
                        print('Invalid Expression, Try again.')
                    self.lp_func[link] = lambdify(x, func_inp, 'numpy')
            else:
                print('Could not modify Link Performance Functions, Try again.') 

    def get_lp_funcs(self):
        lp_funcs_dict = {}
        
        logging.info('Link Performance functions:')
        for link, lp_func in self.lp_func.items():
            lp_funcs_dict[link] = inspect.getsource(lp_func).partition("return")[2].strip()
            logging.info(f'\t\tLink{link}: f(x) = {lp_funcs_dict[link]}')

        return lp_funcs_dict

    def add_alpha_labels(self, alphas=None, link=None):
        """Adds alpha labels to the given links"""
        self.__add_labels('alpha' , labels=alphas, link=link)

    def __add_labels(self, labelname, labels=None, link=None):
        label_dict = {}
        
        if labels != None:
            if type(labels) == dict:
                for link, label in labels.items():
                    self.A.add(link)
                    label_dict[link] = Rational(str(Fraction(float(label)).limit_denominator()))
            
            elif type(labels) in [int, float] and link != None:
                self.A.add(link)
                label_dict[link] = Rational(str(Fraction(float(labels)).limit_denominator()))
            
            elif link == None:
                print('Did not specify a link')
        else:
            inp = pyip.inputYesNo(f'Manually enter {labelname} labels for each link \n(Enter Yes/No): ')
            if inp != 'no':
                links = list(self.A)
                links.sort()
                for link in links:
                    while True:
                        label_inp = input(f'Enter the {labelname} label for the link {link}: ')
                        if label_inp.isprintable():
                            break
                        print('Invalid label, try again.')
                    label_dict[link] = Rational(Rational(str(Fraction(float(label_inp)).limit_denominator())))
            else:
                print('Could not modify alpha labels, Try again.')

        if labelname == 'alpha':
            self.alphas.update(label_dict)

    def __get_all_degrees(self, links, nodes, degree_type='in'):
        if links == None:
           links = self.A
        if nodes == None:
           nodes = self.__get_nodes_from_links(links=links)
        if degree_type == 'in':
            pos = 1
        elif degree_type == 'out':
            pos = 0
        else:
            raise Exception('Invalid degree_type. Enter \'in\' or \'out\'.')

        degrees = {node: 0 for node in nodes.values()}
        for node in nodes:
            for link in links:
                if node == link[pos]:
                    degrees[node] += 1
        degrees = OrderedDict(sorted(degrees.items(), key=lambda x: x[1]))
        return degrees

    def get_indegrees(self, links=None, nodes=None):

        return self.__get_all_degrees(links, nodes, degree_type='in')

    def get_outdegrees(self, links=None, nodes=None):

        return self.__get_all_degrees(nodes, links, degree_type='out')
    
    def topological_sort(self):
        """Return network type and if the network is an acyclic network then this function
        returns nodes sorted in a topological order"""
        ordered_nodes, i = {}, 0
        nodes, links = deepcopy(self.N), deepcopy(self.A)
        original_nodes = {i: node for i, node in enumerate(self.N, 1)}

        while True:
            indegrees = self.get_indegrees(links, nodes)
            # outdegrees = self.get_outdegrees(links, nodes)

            logging.debug(f'nodes = {nodes}\nlinks = {links}\nindegrees = {indegrees}\n')#outdegrees = {outdegrees}')
            
            for node, indegree in indegrees.items():
                if indegree == 0:
                    i += 1
                    ordered_nodes[i] = nodes.pop(node)
                    for link in links.copy():
                        if node in link:
                            links.discard(link)

            if len(indegrees) == 0:
                self.type, self.ordered_nodes = 'Acyclic Network', OrderedDict(ordered_nodes)
                break
            elif 0 not in set(indegrees.values()) and len(indegrees) > 1:
                self.type, self.ordered_nodes = 'Cyclic Network', original_nodes
                break

        return self.type, self.ordered_nodes
        
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

    def find_all_paths(self, start, end, max_node_repeat=1):
        """
        Returns all the possible paths (in 2 formats) from start Node to end Node.
        Path Formats:    0. Paths with links,    1. Compact paths.
        Usage--
        >>> Network.find_all_paths(start=start Node, end=end Node [, max_node_repeat=1])[index_of_path_format]
        """

        self.pi = set()
        b_v = {node: [] for node in self.N}
        for node in self.N:
            for (i, j) in self.forward_star(node):
                b_v[j].append(i)
        #logging.debug(b_v)
        
        # BFS search
        queue = [[end]]
        while queue:
            temp_path = queue.pop(0)
            if temp_path[-1] == start:
                new_path = tuple(temp_path[::-1])
                #logging.debug(f'Adding to path: {new_path}')
                self.pi.add(new_path)
            else:
                for node in b_v[temp_path[-1]]:
                    temp_path.append(node)
                    #logging.debug(temp_path)
                    if temp_path.count(node) > max_node_repeat:
                        #logging.debug('Cyclic Network, so discard this path.')
                        temp_path.pop()
                        continue
                    queue.append(temp_path.copy())
                    temp_path.remove(node)

        alt_paths = set()
        for path in self.pi:
            temp_alt_path = []
            for i in range(len(path)-1):
                temp_alt_path.append((path[i], path[i+1]))
            alt_paths.add(tuple(temp_alt_path))

        return alt_paths, self.pi

    def all_routes(self):
        """Returns a generator of all possible combinations of routes. Example:
        >>> G = Network(4)
        >>> list(G.all_routes())
        [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
        """
        return combinations(self.N.values(), 2)

    def link_path_adjacency_matrix(self, repeat_node=1, routes=None):
        """Returns Link Path Adjacency matrix, rows (list of links) and
        columns (list of paths) of the matrix.
        Usage:
        >>> matrix, link_rows, path_cols = G.link_path_adjacency_matrix()
        >>> matrix
        array([[1, 0, 0, 1, 0, 0],
               [0, 1, 1, 0, 0, 0],
               [1, 0, 0, 0, 0, 1],
               [0, 1, 0, 1, 1, 0],
               [0, 1, 0, 0, 0, 0],
               [1, 0, 1, 0, 0, 1]])
        >>> link_rows
        [(1, 2), (1, 3), (2, 3), (2, 4), (3, 2), (3, 4)]
        >>> path_cols
        [(1, 3, 2), (1, 2), (1, 2, 3), (1, 3), (1, 2, 3, 4), (1, 3, 2, 4),
        (1, 2, 4), (1, 3, 4), (2, 3), (2, 4), (2, 3, 4), (3, 4), (3, 2, 4)]
        """
        paths_all = {}
        p_all = {}
        if routes == None:
            routes = self.all_routes()
        else:
            routes = product(*routes)
            
        # G.all_routes = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
        for route in routes:
            paths_all[route], p_all[route] = self.find_all_paths(*route, max_node_repeat=repeat_node)
            # example: 1st paths_i2j = paths_1to2 = {((1, 2),), ((1, 3), (3, 2))}

        path_col = []

        lp_mat = []
        link_row = list(self.A)
        link_row.sort()

        for n, link in enumerate(link_row):
            row = []
            for pathVal in paths_all.values():
                for p in pathVal:
                    if n == 0:
                        path_col.append(p)
                    if link in p:
                        #logging.debug(f'link= {link}     path={p}')
                        row.append(1)
                    else:
                        row.append(0)
            lp_mat.append(row)
            
        lp_mat = np.array(lp_mat)

        return lp_mat, link_row, path_col


def topological_sort(links):
    """Return network type and if the network is an acyclic network then this function
        returns nodes sorted in a topological order"""
    
    nodes = {node for link in links for node in link}
    original_nodes = {i: node for i, node in enumerate(nodes, 1)}

    def __get_all_degrees(nodes, links, degree_type='in'):
        if degree_type == 'in':
            pos = 1
        elif degree_type == 'out':
            pos = 0
        else:
            raise Exception('Invalid degree_type. Enter \'in\' or \'out\'.')

        degrees = {node: 0 for node in nodes.values()}
        for node in nodes:
            for link in links:
                if node == link[pos]:
                    degrees[node] += 1
        degrees = OrderedDict(sorted(degrees.items(), key=lambda x: x[1]))
        return degrees

    def get_indegrees(nodes, links):

        return __get_all_degrees(nodes, links, degree_type='in')

    def get_outdegrees(nodes, links):

        return __get_all_degrees(nodes, links, degree_type='out')
        
    def main_sorting(links):    
        ordered_nodes, i = {}, 0
        nodes = {node: node for link in links for node in link}

        while True:
            indegrees = get_indegrees(nodes, links)
            outdegrees = get_outdegrees(nodes, links)
            #logging.debug(f'indegrees = {indegrees}\noutdegrees = {outdegrees}')
            
            for node, indegree in indegrees.items():
                if indegree == 0:
                    i += 1
                    node_key = list(nodes.keys())[list(nodes.values()).index(node)]
                    ordered_nodes[i] = nodes.pop(node_key)
                    for link in links.copy():
                        if node in link:
                            links.discard(link)

            if len(indegrees) == 0:
                return 'Acyclic Network', OrderedDict(ordered_nodes)
            elif 0 not in set(indegrees.values()) and len(indegrees) > 1:
                return 'Cyclic Network', original_nodes
                
    return main_sorting(links)


def get_h(paths, d):
    divider = {}        # Assuming all the link flows are divided equally between links
    for path in paths:
        for k in d:
            if (path[0][0], path[-1][-1]) == k:
                if k not in divider:
                    divider.setdefault(k, 1)
                else: 
                    divider[k] += 1

    h_values = {}
    for path in paths:
        for k in d:
            if (path[0][0], path[-1][-1]) == k:
                h_values[path] = d[k] / divider[k]

    return h_values

def prettify_matrix(matrix, indent=8, distance=5):
    """Returns pretty matrix"""

    matrixx = (' '*indent+ '+- ' + ' '*(distance*(len(matrix[0]))) + ' -+\n')
    for inner_list in matrix:
        matrixx +=  (' '*indent+'|')
        for elem in inner_list:
            matrixx += (str(elem).rjust(distance))
        matrixx += ('|'.rjust(distance)) + '\n'
    matrixx += (' '*indent+ '+- ' + ' '*(distance*(len(matrix[0]))) + ' -+')

    return matrixx

def get_xhtc_and_lpaMatrix(G, path_flows, routes, h_vector_dict=None):
    """Returns x\u20D7, t\u20D7, c\u20D7 and the link path adjacency matrix
    Example:
    >>> x, h, t, c, lpa_mat, link_row, path_column, report = get_xhtc_and_lpaMatrix(Graph,\
                                                                    path_flows, routes)
    >>> lpa_mat
    array([[1, 0, 0, 1, 0, 0],
           [0, 1, 1, 0, 0, 0],
           [1, 0, 0, 0, 0, 1],
           [0, 1, 0, 1, 1, 0],
           [0, 1, 0, 0, 0, 0],
           [1, 0, 1, 0, 0, 1]])
    >>> link_row
    [(1, 2), (1, 3), (2, 3), (2, 4), (3, 2), (3, 4)]
    >>> path_column
    [(1, 2, 3, 4), (1, 3, 2, 4), (1, 3, 4), (1, 2, 4), (2, 4), (2, 3, 4)]
    >>> x
    {(1, 2): 20.0, (1, 3): 20.0, (2, 3): 40.0, (2, 4): 50.0, (3, 2): 10.0, (3, 4): 50.0}
    >>> h
    {((1, 2), (2, 3), (3, 4)): 10.0, ((1, 3), (3, 2), (2, 4)): 10.0, ((1, 3), (3, 4)): 10.0,
    ((1, 2), (2, 4)): 10.0, ((2, 4),): 30.0, ((2, 3), (3, 4)): 30.0}
    >>> t
    {(1, 2): 200.0, (1, 3): 70.0, (2, 3): 50.0, (3, 2): 20.0, (2, 4): 100.0, (3, 4): 500.0}
    >>> c
    {(1, 2, 3, 4): 750.0, (1, 3, 2, 4): 190.0, (1, 3, 4): 570.0, (1, 2, 4): 300.0,
    (2, 4): 100.0, (2, 3, 4): 550.0}
    """
    # Start of Calculations:
    lpa_matrix, link_row, path_col = G.link_path_adjacency_matrix(routes=routes)
    alt_path_col = [tuple([path[0][0]] + list(np.array(path)[:,1])) for path in path_col]
    report = '\u0394 =\n'
    report += prettify_matrix(lpa_matrix)

    # Getting link flows for each links as h_vector
    if h_vector_dict == None:
        h_vector_dict = get_h(path_col, path_flows)
        
    h_vector = [h_vector_dict[path] for path in path_col]
        
    # Dot product of delta(lpa_matrix) and h_vector
    report += '\n\nDot product of \u0394 and h gives X. So,\n'
    report += f'    X = \u0394 . h\n==> X = \n'
    x_vector = lpa_matrix.dot(h_vector)
    x_vector_dict = dict(zip(link_row, x_vector))
    report += pformat(x_vector_dict, indent=8)

    # Getting t_vector from the link performance functions
    report += '\n\nSubstituting the X values in each link\'s link performance function ==>\n'
    report += '    t = \n' 
    t_vector_dict = {link: float(lp_func(x_vector_dict[link])) for link, lp_func in G.lp_func.items()}
    report += pformat(t_vector_dict, indent=8)
    t_vector = [t_vector_dict[link] for link in link_row]

    # Getting the c_vector
    report += '\n\nDot product of \u0394^T and t gives c.\n'
    report += f'    c = \u0394^T . t\n==> c = \n'
    c_vector = lpa_matrix.transpose().dot(np.array(t_vector))
    
    c_vector_dict = dict(zip(alt_path_col, c_vector))
    report += pformat(c_vector_dict, indent=8)

    logging.info(report)

    return x_vector_dict, h_vector_dict, t_vector_dict, c_vector_dict, lpa_matrix, link_row, alt_path_col, report

if __name__ == '__main__':

    # Define the Network
    Graph = Network()          # Network of 4 nodes
    Graph.add_links({(2, 4), (1, 2), (3, 4), (2, 3), (3, 2), (1, 3)})
    Graph.add_lpf({(1, 2): '10*x', (1, 3): 'x + 50', (2, 3): 'x + 10', (3, 2): 'x + 10',
               (2, 4): 'x + 50', (3, 4): '10*x'})
    
    Graph.topological_sort()
    print(f'Network type: {Graph.type}\n'
          f'Topological Order of the Nodes: {Graph.ordered_nodes}') 
    

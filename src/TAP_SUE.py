#! python3 

__author__ = "Simanta Barman"
__email__ = "barma017@umn.edu"

from Network import *
from Shortest_path import *
from STOCH_loader import *
from pprint import pprint, pformat

logging.disable(logging.DEBUG)
logging.basicConfig(level=logging.DEBUG, format="%(message)s")


class Tap_SUE:
    def __init__(self, Graph, demand, theta=1):
        self.G = G
        self.demand = demand
        self.theta = theta

    def __load_demand_using_STOCH(self, link_labels, theta):
        """Returns Link flow dictionary after loading the network
        with demand using Dial's STOCH algorithm."""
        return STOCH_loader(link_labels, self.demand, theta).load()
    
    def __get_travel_times(self, link_flows):
        return {link: t(link_flows[link]) for link, t in self.G.lp_func.items()}

    def msa(self, precision=0.02):
        iteration = 1
        # Choosing a initial feasible link assignment 
        x = {link: 0 for link in self.G.A}
            
        while True:
            # Update link travel times
            t = self.__get_travel_times(link_flows=x)
            
            # Calculating target link flows, x* using STOCH algo
            x_star = self.__load_demand_using_STOCH(t, theta=self.theta)
            
            # Convex combination using MSA
            _lambda = 1/(iteration+1)
            x = {link: round(_lambda*x_star[link] + (1 - _lambda) * x[link], 6) for link in self.G.A}

            # Termination criteria
            if abs(sum(x[link] - x_star[link] for link in self.G.A) / len(x)) <= precision:
                print(f'{iteration = }')
                print(f'TSTT = {sum(t[ij]*x[ij] for ij in self.G.A)}')
                return x
            
            iteration += 1

if __name__ == "__main__":
    def sheffi_11_1():
        # Defining the network
        G = Network()

        links = {(1, 2), (5, 8), (1, 4), (2, 3), (4, 5), (8, 9), (5, 6), (3, 6), (2, 5), (6, 9), (4, 7), (7, 8)}
        lpfs = {(i, j): f"{((i+j)/2)} + {((j-i)/100)} * x**2" for (i, j) in links}
        
        G.add_links(links)
        G.add_lpf(lpfs)

        # Defining demand for origin to destination
        d = {(1, 9): 12, (4, 9): 4}

        return G, d

    G, d = sheffi_11_1()
    
    Tap_SUE = Tap_SUE(Graph=G, demand=d, theta=1e10)
    x_ij = Tap_SUE.msa()
    pprint(x_ij)

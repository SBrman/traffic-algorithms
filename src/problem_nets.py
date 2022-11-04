from Network import *


def small_network():
    # Defining the network
    G = Network()
    G.add_links({(1, 2), (1, 3), (3, 2)})
    G.add_lpf({(1, 2): '10+x', (1, 3): '20+x', (3, 2): '0'})

    # Defining demand for origin to destination
    d = {(1, 2): 50}
    routes = [(1,), (2,)]

    return G, d, routes

def large_network():
    # Defining the network
    G = Network()
    G.add_links({(2, 4), (1, 5), (6, 4), (6, 3), (5, 6), (2, 5), (1, 3)})
    G.add_lpf({(1, 3): 'x/100 + 10', (1, 5): 'x/100 + 10', (2, 4): 'x/100 + 10', (2, 5): 'x/100 + 10',
               (5, 6): 'x/100 + 10', (6, 3): 'x/100 + 10', (6, 4): 'x/100 + 10'})
    
    # Defining demand for origin to destination (OD matrix)
    d = {(1, 3): 5000, (1, 4): 0, (2, 3): 0, (2, 4): 10000}
    routes = [(1,2), (3,4)]

    return G, d, routes

def new_net():
    
    links = {(1, 2), (5, 8), (1, 4), (2, 3), (4, 5), (8, 9), (5, 6), 
             (3, 6), (2, 5), (6, 9), (4, 7), (7, 8)}
    
    lp_funcs = {(1, 2): '((1/40000)*x**2 + 3)', (1, 4): '((1/40000)*x**2 + 3)', 
                (2, 3): '((1/40000)*x**2 + 3)', (2, 5): '((1/10000)*x**2 + 5)', 
                (3, 6): '((1/40000)*x**2 + 3)', (4, 5): '((1/10000)*x**2 + 5)', 
                (4, 7): '((1/40000)*x**2 + 3)', (5, 6): '((1/10000)*x**2 + 5)', 
                (5, 8): '((1/10000)*x**2 + 5)', (6, 9): '((1/40000)*x**2 + 3)', 
                (7, 8): '((1/40000)*x**2 + 3)', (8, 9): '((1/10000)*x**2 + 5)'}
    
    alpha_labels = None

    # Create the network
    G = Network()
    G.add_links(links)
    G.add_lpf(lp_funcs)
    G.pos = {7: (0, 2), 8: (1, 2), 9: (2, 2),
             4: (0, 1), 5: (1, 1), 6: (2, 1),
             1: (0, 0), 2: (1, 0), 3: (2, 0)}
    
    # Defining demand for origin to destination (OD matrix)
    d = {(4, 9): 1000, (1, 9): 1000}
    routes = [(4,1), (9,)]
    return G, d, routes

def fig7():
    links = {(7, 4), (1, 2), (9, 6), (4, 2), (2, 3), (4, 5),
             (8, 9), (6, 3), (4, 1), (7, 8), (5, 2)}

    lp_funcs = {(1, 2): '(2*x**2 + 4)', (2, 3): '(x**2 + 2)',
                (4, 1): '(2*x**2 + 4)', (4, 2): '(2*x**2 + 40)',
                (4, 5): '(x**2 + 2)', (5, 2): '(x**2 + 2)',
                (6, 3): '(x**2 + 2)', (7, 4): '(x**2 + 2)',
                (7, 8): '(x**2 + 2)', (8, 9): '(x**2 + 2)',
                (9, 6): '(x**2 + 2)', (8, 5): '(x**2 + 2)',
                (5, 6): '(x**2 + 2)'}

    # Create the network
    G = Network()
    G.add_links(links)
    G.add_lpf(lp_funcs)
    G.pos = {1: (0, 2), 2: (1, 2), 3: (2, 2),
             4: (0, 1), 5: (1, 1), 6: (2, 1),
             7: (0, 0), 8: (1, 0), 9: (2, 0)}
    
    # Defining demand for origin to destination (OD matrix)
    d = {(7, 3): 10}
    routes = [(7,), (3,)]
    
    return G, d, routes

def SiouxFalls():
    G = Network()
    SF = GetNetworkData(name='SiouxFalls')
    G.add_links(SF.links)
    G.add_lpf(SF.lp_func)
    
    unformatted_demands = SF.demands
    d = {}
    for i in G.N:
        for j in G.N:
            demand_ij = unformatted_demands[i][j]
            d[(i, j)] = demand_ij
                
    routes = [tuple({r for r, _ in d}), tuple({s for _, s in d})]
        
    return G, d, routes

def sheffi_11_1():
    # Defining the network
    G = Network()

    links = {(1, 2), (5, 8), (1, 4), (2, 3), (4, 5), (8, 9), (5, 6), (3, 6), (2, 5), (6, 9), (4, 7), (7, 8)}
    lpfs = {(i, j): f"{((i+j)/2)} + {((j-i)/100)} * x**2" for (i, j) in links}
    
    G.add_links(links)
    G.add_lpf(lpfs)

    # Defining demand for origin to destination
    d = {(1, 9): 12, (4, 9): 4}

    return G, d, None
    
def sheffi_ex():
    # Defining the network
    G = Network()
    
    G.add_links({(1, 2), (2, 3), (1, 12), (12, 2), (1, 3)})
    G.add_lpf({(1, 2): '2+x**2', (1, 12): '1+3*x', (12, 2): '0',
               (2, 3): '3+x', (1, 3): '3.5 + 0.2*x'})

    # Defining demand for origin to destination
    d = {(1, 3): 4}
    routes = [(1,), (3,)]

    return G, d, routes

def sheffi_ex_1():
    # Defining the network
    G = Network()
    G.add_links({(1, 2), (2, 3), (1, 12), (12, 2), (1, 3)})
    G.add_lpf({(1, 2): '2+3*x**2', (1, 12): '1+6*x', (12, 2): '0',
               (2, 3): '3+2*x', (1, 3): '3.5 + 0.4*x'})

    # Defining demand for origin to destination
    d = {(1, 3): 4}
    routes = [(1,), (3,)]

    return G, d, routes

def neww():
    # Defining the network
    G = Network()

    lpfs = {(1, 2): '10+0.5*x', (2, 23): '20+0.1*x**2', (23, 3): '0',
               (2, 3): '30+0.05*x**2', (1, 3): '40+0.002*x**3'}
    G.add_lpf(lpfs)
    G.add_links(lpfs.keys())

    # Defining demand for origin to destination
    d = {(1, 3): 25, (2, 3): 8}

    return G, d, None

def p1():
    # Defining the network
    G = Network()
    
    # tdsp1 BLU book
    link_cost = {
                    (1, 3): lambda t: 10,
                    (1, 2): lambda t: 4 + t,
                    (2, 3): lambda t: max(5 - t/2, 1),
                    (2, 4): lambda t: 5,
                    (3, 4): lambda t: t / 2
                }
    
    G.link_cost = link_cost
    G.add_links(list(link_cost.keys()))

    return G


def p2():
    # Defining the network
    G = Network()
    # tdsp2 BLU book
    link_cost = {
                    (1, 2): lambda t: 1,
                    (1, 3): lambda t: 1,
                    (2, 3): lambda t: 1,
                    (3, 4): lambda t: 5 if t == 2 else 10 
                }

    G.tau = {
                    (1, 2): lambda t: 1,
                    (1, 3): lambda t: 3,
                    (2, 3): lambda t: 1,
                    (3, 4): lambda t: 1
            }
    
    G.link_cost = link_cost
    G.add_links(list(link_cost.keys()))

    return G

def e2():
    # Defining the network
    G = Network()
    # tdsp2
    link_cost = {
                    ('a', 'b'): lambda t: [1, 2, 1, 3, 2, 4, 3][t],
                    ('a', 'ab'): lambda t: [0, 0, 0, 0, 0, 0, 0][t],
                    ('ab', 'b'): lambda t: [2, 1, 3, 2, 1, 2, 3][t],
                    ('a', 'c'): lambda t: [5, 6, 5, 4, 5, 7, 6][t],
                    ('b', 'c'): lambda t: [4, 3, 5, 3, 2, 4, 3][t] 
                }

    G.tau = link_cost.copy()
    
    G.link_cost = link_cost
    G.add_links(list(link_cost.keys()))

    return G

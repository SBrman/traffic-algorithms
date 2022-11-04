#! python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx

def plot_net(links, xx=None, tt=None, title=None, pos=None):

    lxpos = 0.35 if tt else 0.5
    ltpos = 0.65 if xx else 0.5
    
    f, ax = plt.subplots(figsize=(16, 9))

    G = nx.DiGraph()
    G.add_edges_from(links)

    nodes = {i for i, _ in links}.union({j for _, j in links})
    if not pos:
        pos = {node: np.random.randint(0, 10000, size=(1, 2)) for node in nodes}
            
    # edges
    nx.draw_networkx(G, pos=pos, with_labels=True, node_color='lightblue', font_size=15, node_size=600)
    if xx:
        nx.draw_networkx_edge_labels(G, pos, edge_labels={node: 'x='+str(x) for node, x in xx.items()}, label_pos=lxpos, font_size=14, font_color='b', bbox=dict(facecolor='w', edgecolor='b'))
    if tt:
        nx.draw_networkx_edge_labels(G, pos, edge_labels={node: 't='+str(t) for node, t in tt.items()}, label_pos=ltpos, font_size=14, font_color='purple', bbox=dict(facecolor='w', edgecolor='purple'))
    if title: 
        plt.suptitle(f'Figure: {title}', bbox=dict(facecolor='None', alpha=0.5), y=0.01, va='bottom')
    #plt.title('Link Flows and Travel times for each link.')
    plt.subplots_adjust(left=0.03, bottom=0.05, right=0.97, top=0.95, wspace=None, hspace=None)
    return ax

if __name__ == '__main__':
    f, ax = plt.subplots(figsize=(16,9))
    A = {(3, 4), (4, 3), (3, 1), (18, 20), (5, 4), (16, 17), (12, 13), (21, 22), (22, 20), (22, 23), (9, 5), (23, 22), (8, 9), (9, 8), (8, 6), (10, 9), (11, 14), (1, 3), (19, 15), (10, 15), (6, 2), (16, 10), (18, 7), (6, 5), (15, 14), (24, 23), (6, 8), (18, 16), (12, 3), (4, 5), (20, 19), (5, 6), (21, 24), (3, 12), (5, 9), (4, 11), (20, 22), (14, 15), (11, 4), (23, 24), (9, 10), (1, 2), (10, 11), (2, 1), (11, 10), (19, 17), (19, 20), (10, 17), (24, 13), (15, 10), (15, 19), (16, 18), (15, 22), (12, 11), (20, 18), (7, 18), (23, 14), (14, 11), (20, 21), (21, 20), (22, 15), (22, 21), (17, 10), (8, 7), (14, 23), (17, 16), (17, 19), (8, 16), (11, 12), (10, 16), (2, 6), (13, 12), (16, 8), (24, 21), (7, 8), (13, 24)}
    pos = {1: (-96.77041974, 43.61282792), 2: (-96.71125063, 43.60581298), 3: (-96.77430341, 43.5729616), 4: (-96.74716843, 43.56365362), 5: (-96.73156909, 43.56403357), 6: (-96.71164389, 43.58758553), 7: (-96.69342281, 43.5638436), 8: (-96.71138171, 43.56232379), 9: (-96.73124137, 43.54859634), 10: (-96.73143801, 43.54527088), 11: (-96.74684071, 43.54413068), 12: (-96.78013678, 43.54394065), 13: (-96.79337655, 43.49070718), 14: (-96.75103549, 43.52930613), 15: (-96.73150355, 43.52940117), 16: (-96.71138171, 43.54674361), 17: (-96.71138171, 43.54128009), 18: (-96.69407825, 43.54674361), 19: (-96.71131617, 43.52959125), 20: (-96.71118508, 43.5153335), 21: (-96.7309792, 43.51048509), 22: (-96.73124137, 43.51485818), 23: (-96.75090441, 43.51485818), 24: (-96.74920028, 43.50316422)}
    plot_net(A, None, None, title='Network', pos=pos)
    plt.savefig(r'.\\Net_x_t_initial.png', dpi=120)

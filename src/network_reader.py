#! python3

def data_reader(path):    
    links = {}
    link_cost = {}
    with open(path, 'r') as file:
        for idx, line in enumerate(file.readlines()):
            if idx == 0:
                continue
            link_id, i, j, *t = line.strip('\n').split('\t')
            i, j = (int(i), int(j))
            links[int(link_id)] = (i, j)
            link_cost[i, j] = {i: int(tt) for i, tt in enumerate(t)}

    return links, link_cost
    

import json
from networkx.readwrite import json_graph
from pandas import read_csv
import csv

def read_graph_from_json(json_filename, update_population=True, set_edge_weights=True):
    print(f"Reading graph from {json_filename}...")
    with open(json_filename) as f:
        data = json.load(f)
    G = json_graph.adjacency_graph(data) 
    if update_population:
        for i in G.nodes:
            G.nodes[i]['TOTPOP'] = G.nodes[i]['P0010001'] 
    if set_edge_weights:
        for e in G.edges:
            G.edges[e]['weight'] = 1
    return G

def read_plan_from_csv(csv_filename):
    geoid_to_district = dict()
    unassigned = list()
    f = read_csv(csv_filename, skipinitialspace=True)
    for index, row in f.iterrows():
        g = str( row['GEOID20'] ) # get GEOID
        if len(g) == 14: # fix issue with leading zeros in block-level plans
            g = '0' + g
        d = str( row['District'] ) # get district
        if d in { 'ZZ', 'ZZZ' }:
            unassigned.append(g)
        else:
            geoid_to_district[g] = int(d)-1
    if unassigned:
        print("The following GEOIDs are unassigned:",unassigned)
    return geoid_to_district

def write_plan_to_csv(geoid_to_district, csv_filename):
    print(f"Writing plan to {csv_filename}...")
    with open(csv_filename, 'w', newline="") as f: 
        csvwriter = csv.writer(f) 
        for geoid, district in geoid_to_district.items():
            row = [ str(geoid), district+1 ]
            csvwriter.writerow(row) 
    return

def append_dict_as_row(filename, dict_of_elem, field_names):
    with open(filename, 'a+', newline='') as write_obj:
        dict_writer = csv.DictWriter(write_obj, fieldnames=field_names)
        dict_writer.writerow(dict_of_elem)
        
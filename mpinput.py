# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 17:46:32 2020

@author: rundx
"""

import numpy as np
import networkx as nx



#Input matpower instance file
filename = "14-bus.txt"

def load_data(filename):
    #Read file to learn appropriate number of Buses to use and Rows to Skip
    #Finds tags for bus/gen/branch data, counts numbers, and breaks early if all data is found
    with open(filename) as f:
        read_data = f.readlines()
    
    
    
    #Cuts through and generates bus/gen/branch data from the text file
    #Outputs are as basic numpy arrays 
    line_count = 0
    for line in read_data:
        line_count += 1
        s = line.split()
        if len(s) > 1 and s[0] == 'mpc.bus':
            bus_start_line = line_count
            bus_skip = bus_start_line
        if len(s) > 1 and s[0] == 'mpc.gen':
            bus_end_line = line_count
            num_buses = bus_end_line - bus_start_line - 5
            gen_start_line = bus_end_line
            gen_skip = gen_start_line
        if len(s) > 2 and s[0] == 'mpc.gencost':
            gen_cost_start = line_count
            gen_end_line = gen_cost_start
            num_gens = gen_end_line - gen_start_line - 5
        if len(s) > 2 and s[0] == 'mpc.branch':
            gen_cost_end_line = line_count
            branch_start_line = gen_cost_end_line
            branch_skip = branch_start_line
        if len(s) > 2 and s[1] == 'INFO':
            branch_end_line = line_count
            num_branches = branch_end_line - branch_start_line - 3
            break


    
    
    #Comprehension tool to remove semicolons from last data entry in each row of the file
    convert = lambda x: float(x.decode('ascii').strip(';') or 0)
    
    
    #Read bus, gen, branch data into numpy arrays
    
    bus_dtype = [('bus_i', int),('bus_type', int),('bus_pd',float),('bus_qd',float),
                 ('bus_Gs',float),('bus_Bs', float),('bus_area',int),('bus_Vm',float),
                 ('bus_Va',float),('bus_baseKV',float),('bus_zone',int),
                 ('bus_Vmax',float),('bus_Vmin',float)]
    gen_dtype = [('bus_i',int),('Pg',float),('Qg',float),('Qmax',float),
                 ('Qmin',float),('Vg',float),('mBase',float),('status',int),('Pmax',float),
                 ('Pmin',float),('Pc1',float)]
    branch_dtype = [('fbus',int),('tbus',int),('r',float),('x',float),('b',float),
                    ('rateA',float),('rateB',float),('rateC',float),('ratio',float),
                    ('angle',float),('status',float),('angmin',float),('angmax',float)]
    gen_cost_dtype = [('gen_type', int),('startup',float),('shutdown',float),('something',int),
                      ('quad_cost',float),('lin_cost',float),('flat_cost',float)]
    bus_data = np.genfromtxt(filename, skip_header = bus_skip, 
                             autostrip= True, max_rows = num_buses, dtype = bus_dtype,
                             converters = {-1: convert})
    gen_data = np.genfromtxt(filename, skip_header = gen_skip, 
                             autostrip= True, max_rows = num_gens, dtype = gen_dtype,
                             converters = {-1: convert}, filling_values = 0)
    branch_data = np.genfromtxt(filename, skip_header = branch_skip,
                                autostrip= True, max_rows = num_branches, dtype = branch_dtype,
                                converters = {-1: convert})
    gen_cost_data = np.genfromtxt(filename, skip_header = gen_cost_start, dtype = gen_cost_dtype,
                                  autostrip = True, max_rows = num_gens, 
                                  converters = {-1:convert}, filling_values = 0)
    idx = np.argsort(bus_data['bus_i'])
    bus_data['bus_i'] = np.array(bus_data['bus_i'])[idx]
    bus_data['bus_pd'] = abs(np.array(bus_data['bus_pd'])[idx])
    #Generator costs are typically embedded in the OPF data at the end of the file
    '''gen_cost = np.genfromtxt(filename, skip_header = gen_cost_start, 
                             autostrip = True, max_rows = num_gens, converters = {5:convert},
                             usecols = (5))
    
    
    gen_data['Pc1'] = gen_cost'''
    gen_data['Pc1'] = np.array(gen_cost_data['lin_cost'])
    if not np.any(gen_data['Pc1']):
        gen_data['Pc1'] = 1

    return bus_data,gen_data,branch_data


#Modifies existing networkx object, graph, to assign nodes, edges, and all
#commonly used attributes:
#bus_i, bus_pd
def encode_graph(graph, bus_data, gen_data, branch_data, demand_scale, capacity_scale, gen_scale, cost_scale):
    
    #add nodes to the graph from this list of bus data names
    #create two dictionaries, each with keys from the list of nodes
    #but with values from list of nodes and from demand values
    #add these named attributes to networkx node objects
    graph.add_nodes_from(bus_data['bus_i'])
    bus_i_dict = dict(zip(bus_data['bus_i'],bus_data['bus_i']))
    bus_pd_dict = dict(zip(bus_data['bus_i'],bus_data['bus_pd']))
    bus_zone_dict = dict(zip(bus_data['bus_i'],bus_data['bus_area']))
    #bus_pd_dict = {key: .1*value for (key,value) in nx.get_node_attributes(graph,'bus_pd').items()}
    nx.set_node_attributes(graph, bus_i_dict, 'bus_i')
    #Demand Modifier scales input demand values by multiplying
    demand_modifier = 1*demand_scale
    nx.set_node_attributes(graph, {key: demand_modifier*value for (key,value) in bus_pd_dict.items()}, 'bus_pd')
    nx.set_node_attributes(graph, {key: value for (key,value) in bus_zone_dict.items()}, 'bus_zone')
    
    
    #Simple loop through generators to add generator properties as data to nodes for which 
    #generators exist; always sets minimum generation to 0
    gen_modifier = gen_scale
    #gen_modifier = 1
    j = 0
    gen_cost_scale = cost_scale
    for i in gen_data['bus_i']:
        if not graph.nodes[i].get('gen_Pmax'):
            graph.nodes[i]['gen_Pmax'] = gen_data['Pmax'][j]*gen_modifier
            graph.nodes[i]['gen_Pmin'] = gen_data['Pmin'][j]
            graph.nodes[i]['gen_cost'] = gen_data['Pc1'][j]*gen_cost_scale
        else:
            graph.nodes[i]['gen_Pmax'] += gen_data['Pmax'][j]*gen_modifier
            graph.nodes[i]['gen_Pmin'] += gen_data['Pmin'][j]
            graph.nodes[i]['gen_cost'] += gen_data['Pc1'][j]*gen_cost_scale
        #cost default to 1 since most of these datasets are missing this information
        j += 1
    for i in bus_data['bus_i']:
        if not graph.nodes[i].get('gen_Pmax'):
            graph.nodes[i]['gen_Pmax'] = 0
            graph.nodes[i]['gen_Pmin'] = 0
            graph.nodes[i]['gen_cost'] = 0

    
    #add edges to the graph based on fbus and tbus data
    ### NOTE NOTE NOTE ###
    #to reuse a zip: save list(zip(whatever,whatever)) as a variable and use that to 
    #be able to iterate through repeatedly
    '''
    for i in range(len(branch_data['rateC'])):
        if branch_data['rateC'][i] == 0:
            branch_data['rateC'][i] = 660
    '''      
            
    #Order all edges properly:
    for i in range(len(branch_data['fbus'])):
        if branch_data['fbus'][i] > branch_data['tbus'][i]:
            branch_data['fbus'][i], branch_data['tbus'][i]= branch_data['tbus'][i], branch_data['fbus'][i]
            
    
    
    graph.add_edges_from(zip(branch_data['fbus'],branch_data['tbus']))
    branch_x_dict = dict(zip(zip(branch_data['fbus'],branch_data['tbus']),abs(branch_data['x'])))
    nx.set_edge_attributes(graph, branch_x_dict, 'branch_x')
    branch_b_dict = dict(zip(zip(branch_data['fbus'],branch_data['tbus']),abs(branch_data['b'])))
    nx.set_edge_attributes(graph, branch_b_dict, 'branch_b')
    for edge in graph.edges:
        if graph.edges[edge]['branch_b'] == 0:
            graph.edges[edge]['branch_b'] = graph.edges[edge]['branch_x']
    cap_modifier = 1*capacity_scale
    branch_cap_dict = dict(zip(zip(branch_data['fbus'],branch_data['tbus']), cap_modifier*branch_data['rateC']))
    nx.set_edge_attributes(graph, branch_cap_dict, 'branch_cap')
    branch_CR_dict = dict(zip(zip(branch_data['fbus'],branch_data['tbus']), cap_modifier*branch_data['rateC']*abs(branch_data['x'])))
    nx.set_edge_attributes(graph, branch_CR_dict, 'branch_CR')
    branch_neg_CR_dict = dict(zip(zip(branch_data['fbus'],branch_data['tbus']), -1*cap_modifier*branch_data['rateC']*abs(branch_data['x'])))
    nx.set_edge_attributes(graph, branch_neg_CR_dict, 'branch_neg_CR')
    branch_cap_dict = dict(zip(zip(branch_data['fbus'],branch_data['tbus']), branch_data['status']))
    nx.set_edge_attributes(graph, branch_cap_dict, 'branch_status')
    branch_exp_cost_dict = dict(zip(zip(branch_data['fbus'],branch_data['tbus']), branch_data['rateA']))
    nx.set_edge_attributes(graph,branch_exp_cost_dict,'branch_exp_cost')
    branch_cand_cost_dict = dict(zip(zip(branch_data['fbus'],branch_data['tbus']), branch_data['rateB']))
    nx.set_edge_attributes(graph,branch_cand_cost_dict,'branch_cand_cost')

'''

filename = 'pglib-opf-master/pglib_opf_case14_ieee.m'

bus_data = np.array([])
gen_data = np.array([])
branch_data = np.array([])
bus_data,gen_data,branch_data = load_data(filename)

graph = nx.Graph()
encode_graph(graph, bus_data, gen_data, branch_data)'''
# -*- coding: utf-8 -*-
"""
Created on Fri May 28 12:09:38 2021

@author: Kyle Skolfield
"""


import gurobipy as gp
import networkx as nx
import numpy as np
import mpinput as mp
import time
#import cProfile
import sys
from itertools import islice
#from itertools import combinations
import random
import math
#import nicify_tree_decomp as ntd
#import os


########################################
########### HELPER FUNCTIONS ###########
########################################


def path_lengths(graph,path_set):
    lengths = []
    for path in path_set:
        length = 0
        for edge in path:
            length += graph.edges[edge]['branch_CR']
        lengths.append(length)
    return lengths

def k_shortest_paths(G, source, target, k, weight=None):
    return list(islice(nx.shortest_simple_paths(G, source, target, weight=weight), k))

#returns edge list of a decently long path between start node and end node
#may not actually reach end node
#what to do if it doesn't reach?
def path_expander(graph,start_node,end_node):
    used_v = []
    used_v.append(start_node)
    curr_node = start_node
    main_subgraph = graph.copy()
    while len(main_subgraph.nodes) > 1:
        #print(curr_node, list(main_subgraph.nodes))
        curr_weight = 0
        neighbors = nx.neighbors(main_subgraph,curr_node)
        main_subgraph.remove_node(curr_node)
        for adj_node in neighbors:
            reachable_nodes = nx.dfs_preorder_nodes(main_subgraph,adj_node)
            temp_subgraph = nx.subgraph(main_subgraph,reachable_nodes)
            temp_weight = temp_subgraph.size(weight = 'branch_CR')
            if  temp_weight >= curr_weight:
                curr_weight = temp_weight
                new_node = adj_node
        used_v.append(new_node)
        if curr_node == new_node or curr_node == end_node:
            break
        curr_node = new_node
    return list(nx.utils.pairwise(used_v))

def induced_path(start_node,graph):
    temp_graph = graph.copy()
    neighbors = list(temp_graph.neighbors(start_node))
    path_nodes = [start_node]
    while len(neighbors) > 0:
        next_node = random.choice(neighbors)
        #print('traveling to node', next_node)
        path_nodes.append(next_node)
        #print('path so far', path_nodes)
        neighbors.remove(next_node)
        for node in neighbors:
            temp_graph.remove_node(node)
            #print('removing node', node)
            neighbors.remove(node)
        neighbors = list(set(temp_graph.neighbors(next_node)) - set(path_nodes))
    path_edges = list(nx.utils.pairwise(path_nodes))
    temp_graph = temp_graph.edge_subgraph(path_edges).copy()
        
    return temp_graph

def long_induced_path(graph, repeats):
    best_length = 0
    #can specificy subset of nodes
    subset = graph.nodes
    for node in subset:
        #temp_best = 0
        #print(node)
        for i in range(repeats):
            path_graph = induced_path(node,graph)
            path_length = path_graph.size(weight='branch_CR')
            if path_length > best_length:
                best_length = path_length
    return best_length

###################################
########### END HELPERS ###########
###################################



########################################
########### BASE DCOPF MODEL ###########
########################################

#For full DCOPF, use flow capacity constraints from OTS and force all switch variables to be constant 1

#Load flow equations
def load_cons_flow_eq(mod,graph):
    M = 2*.6*max({key: 1/value for  (key,value) in nx.get_edge_attributes(graph,'branch_x').items()}.values())
    #M = 2*6.8
    edge_x = nx.get_edge_attributes(graph, 'branch_x')
    mod.addConstrs((-1/edge_x[i,j] * (mod._bus_angle[i] - mod._bus_angle[j] ) - mod._corr_flow[i,j] + (1 - mod._switch[i,j])*M >= 0 for (i,j) in mod._corr_flow),
                   name = "F_eq_pos")
    mod.addConstrs((-1/edge_x[i,j] * (mod._bus_angle[i] - mod._bus_angle[j] ) - mod._corr_flow[i,j] - (1 - mod._switch[i,j])*M <= 0 for (i,j) in mod._corr_flow),
                   name = "F_eq_neg")

#Load bus angle difference constraints
def load_cons_bus_ang_diff(mod,graph,bus_ang):
    mod.addConstrs((mod._bus_angle[i] - mod._bus_angle[j] <= bus_ang for (i,j) in mod._corr_flow))
    mod.addConstrs((mod._bus_angle[j] - mod._bus_angle[i] <= bus_ang for (i,j) in mod._corr_flow))
  

#Load flow balance constraints
def load_cons_balance(mod,graph):
    for i in graph.nodes:
        mod.addConstr(gp.quicksum([mod._corr_flow[j,i] for j in graph.neighbors(i) if j < i]) - 
                      gp.quicksum([mod._corr_flow[i,j] for j in graph.neighbors(i) if j > i]) + 
                      (mod._gen.get(i) or 0) == nx.get_node_attributes(graph, 'bus_pd')[i])

#######################################
########### END DCOPF MODEL ###########
#######################################




#######################################
########### OTS CONSTRAINTS ###########
#######################################

#Load flow capacity constraints
def load_cons_OTS_flow_cap(mod,graph):
    edge_caps = nx.get_edge_attributes(graph,'branch_cap')
    mod.addConstrs((mod._corr_flow[i,j] <= edge_caps[i,j] * mod._switch[i,j] for (i,j) in mod._corr_flow), name="F_cap_pos")
    mod.addConstrs((mod._corr_flow[i,j] >= -edge_caps[i,j] * mod._switch[i,j] for (i,j) in mod._corr_flow), name="F_cap_neg")

#Load maximum number of switched off lines constraint
def load_cons_switch_cap(mod,graph,cap):
    mod.addConstr(len(mod._corr_flow) - gp.quicksum(mod._switch) <=cap, name='switch_cap')
     
###########################################
########### END OTS CONSTRAINTS ###########
###########################################



#######################################
########### TEP CONSTRAINTS ###########
#######################################

def load_cons_flow_eq_TEP(mod,graph):
    M = 2*.6*max({key: 1/value for  (key,value) in nx.get_edge_attributes(graph,'branch_x').items()}.values())
    #M = 2*6.8
    edge_x = nx.get_edge_attributes(graph, 'branch_x')
    mod.addConstrs((-1/edge_x[i,j] * (mod._bus_angle[i] - mod._bus_angle[j] ) - mod._corr_flow[i,j] + 
                    (1 - mod._new_line[i,j])*M >= 0 for (i,j) in mod._corr_flow if edge_status[i,j] == 0),
                   name = "F_eq_pos")
    mod.addConstrs((-1/edge_x[i,j] * (mod._bus_angle[i] - mod._bus_angle[j] ) - mod._corr_flow[i,j] - 
                    (1 - mod._new_line[i,j])*M <= 0 for (i,j) in mod._corr_flow if edge_status[i,j] == 0),
                   name = "F_eq_neg")
    mod.addConstrs((-1/edge_x[i,j] * (mod._bus_angle[i] - mod._bus_angle[j] ) ==  mod._corr_flow[i,j] for (i,j) in mod._corr_flow if edge_status[i,j] == 1),
                   name = "F_eq")


def load_TEP_cons_flow_cap(mod,graph):
    edge_caps = nx.get_edge_attributes(graph,'branch_cap')
    edge_status = nx.get_edge_attributes(graph,'branch_status')
    mod.addConstrs((mod._corr_flow[i,j] <= edge_caps[i,j] * mod._new_line[i,j] for (i,j) in mod._corr_flow if edge_status[i,j] == 0), name="F_cap_pos")
    mod.addConstrs((mod._corr_flow[i,j] >= -edge_caps[i,j] * mod._new_line[i,j] for (i,j) in mod._corr_flow if edge_status[i,j] == 0), name="F_cap_neg")
    mod.addConstrs((mod._corr_flow[i,j] <= edge_caps[i,j] + edge_caps[i,j] * mod._cap_exp[i,j]  for (i,j) in mod._corr_flow if edge_status[i,j] == 1), name="F_cap_pos")
    mod.addConstrs((mod._corr_flow[i,j] >= -edge_caps[i,j] - edge_caps[i,j] * mod._cap_exp[i,j] for (i,j) in mod._corr_flow if edge_status[i,j] == 1), name="F_cap_neg")
    
    
def load_cons_flow_eq_TEP_hyb(mod,graph):
    edge_x = nx.get_edge_attributes(graph, 'branch_x')
    mod.addConstrs((-1/edge_x[i,j] * (mod._bus_angle[i] - mod._bus_angle[j] ) ==  mod._corr_flow[i,j] for (i,j) in mod._corr_flow if edge_status[i,j] == 1),
                   name = "F_eq")

###########################################
########### END TEP CONSTRAINTS ###########
###########################################



def fused_graph_pos(graph, model_list, lin_mod, hyb_mod, trans_mod):
    pos_flow = graph.copy()
    
    for edge in pos_flow.edges:
        if model_list == ['lin']:
            if lin_mod._corr_flow[edge].x <= 0:
                pos_flow.remove_edge(edge[0],edge[1])
        elif model_list == ['hyb']:
            if hyb_mod._corr_flow[edge].x <= 0:
                pos_flow.remove_edge(edge[0],edge[1])
        elif model_list == ['trans']:
            if trans_mod._corr_flow[edge].x <= 0:
                pos_flow.remove_edge(edge[0],edge[1])
        elif 'lin' in model_list and 'hyb' in model_list:
            if lin_mod._corr_flow[edge].x <= 0 or hyb_mod._corr_flow[edge].x <=0:
                pos_flow.remove_edge(edge[0],edge[1])
        elif 'lin' in model_list and 'trans' in model_list:
            if lin_mod._corr_flow[edge].x <= 0 or trans_mod._corr_flow[edge].x <=0:
                pos_flow.remove_edge(edge[0],edge[1])
        elif 'trans' in model_list and 'hyb' in model_list:
            if trans_mod._corr_flow[edge].x <= 0 or hyb_mod._corr_flow[edge].x <=0:
                pos_flow.remove_edge(edge[0],edge[1])
        elif 'lin' in model_list and 'hyb' in model_list and 'trans' in model_list:
            if lin_mod._corr_flow[edge].x <= 0 or hyb_mod._corr_flow[edge].x <=0 or trans_mod._corr_flow[edge].x <= 0:
                pos_flow.remove_edge(edge[0],edge[1]) 
    
    return pos_flow
        

def fused_graph_neg(graph, model_list, lin_mod, hyb_mod, trans_mod):
    neg_flow = graph.copy()
    
    for edge in neg_flow.edges:
        if model_list == ['lin']:
            if lin_mod._corr_flow[edge].x >= 0:
                neg_flow.remove_edge(edge[0],edge[1])
        elif model_list == ['hyb']:
            if hyb_mod._corr_flow[edge].x >= 0:
                neg_flow.remove_edge(edge[0],edge[1])
        elif model_list == ['trans']:
            if trans_mod._corr_flow[edge].x >= 0:
                neg_flow.remove_edge(edge[0],edge[1])
        elif 'lin' in model_list and 'hyb' in model_list:
            if lin_mod._corr_flow[edge].x >= 0 or hyb_mod._corr_flow[edge].x >=0:
                neg_flow.remove_edge(edge[0],edge[1])
        elif 'lin' in model_list and 'trans' in model_list:
            if lin_mod._corr_flow[edge].x >= 0 or trans_mod._corr_flow[edge].x >=0:
                neg_flow.remove_edge(edge[0],edge[1])
        elif 'trans' in model_list and 'hyb' in model_list:
            if trans_mod._corr_flow[edge].x >= 0 or hyb_mod._corr_flow[edge].x >=0:
                neg_flow.remove_edge(edge[0],edge[1])
        elif 'lin' in model_list and 'hyb' in model_list and 'trans' in model_list:
            if lin_mod._corr_flow[edge].x >= 0 or hyb_mod._corr_flow[edge].x >=0 or trans_mod._corr_flow[edge].x >= 0:
                neg_flow.remove_edge(edge[0],edge[1]) 
    
    return neg_flow



##########################################
########### VALID INEQUALITIES ###########
##########################################

def theorem_7_cycle_basis(model,where):
    if where == gp.GRB.Callback.MIPSOL:
        path = random.choice(cycle_paths)
        cycle_paths.remove(path)
        path1 = list(nx.utils.pairwise(path[0]))
        path2 = list(nx.utils.pairwise(path[1]))
        mp_length = graph.edge_subgraph(path1).size(weight = 'branch_CR')
        start = path[0][0]
        end = path[-1][1]
        
        sp_e = path2
        sp_length = graph.edge_subgraph(sp_e).size(weight = 'branch_CR')

        if sp_length>mp_length:
            sp_e,path1 = path1,sp_e
            sp_length,mp_length = mp_length,sp_length
            
        for i in range(len(path1)):
            if path1[i][0] > path1[i][1]:
                path1[i] = tuple(reversed(path1[i]))
        for i in range(len(sp_e)):
            if sp_e[i][0] > sp_e[i][1]:
                sp_e[i] = tuple(reversed(sp_e[i]))
                
        edge_status = nx.get_edge_attributes(graph,'branch_status')
        model.cbLazy(model._bus_angle[start] - model._bus_angle[end] <= sp_length + 
                 (mp_length -sp_length) * (len(sp_e) - gp.quicksum([mod._expansion[edge] for edge in sp_e if edge_status[edge] == 0])))
        model.cbLazy(model._bus_angle[end] - model._bus_angle[start] <= sp_length + 
                 (mp_length -sp_length) * (len(sp_e) - gp.quicksum([mod._expansion[edge] for edge in sp_e if edge_status[edge] == 0])))   
        
        
def theorem_7_heuristic(model,where):
    if where == gp.GRB.Callback.MIPSOL:
        
        #if n[0] <= 0:
        #    return
        if len(path_pairs) == 0:
            return
        path = random.choice(path_pairs)
        path_pairs.remove(path)
        path1 = list(nx.utils.pairwise(path[0]))
        path2 = list(nx.utils.pairwise(path[1]))
        mp_length = graph.edge_subgraph(path1).size(weight = 'branch_CR')
        start = path[0][0]
        end = path[-1][1]
        
        sp_e = path2
        sp_length = graph.edge_subgraph(sp_e).size(weight = 'branch_CR')

        if sp_length>mp_length:
            sp_e,path1 = path1,sp_e
            sp_length,mp_length = mp_length,sp_length
            
        for i in range(len(path1)):
            if path1[i][0] > path1[i][1]:
                path1[i] = tuple(reversed(path1[i]))
        for i in range(len(sp_e)):
            if sp_e[i][0] > sp_e[i][1]:
                sp_e[i] = tuple(reversed(sp_e[i]))
                
        #n[0] = n[0] - 1
        edge_status = nx.get_edge_attributes(graph,'branch_status')
        model.cbLazy(model._bus_angle[start] - model._bus_angle[end] <= sp_length + 
                 (mp_length -sp_length) * (len(sp_e) - gp.quicksum([model._expansion[edge] for edge in sp_e if edge_status[edge] == 0])))
        model.cbLazy(model._bus_angle[end] - model._bus_angle[start] <= sp_length + 
                 (mp_length -sp_length) * (len(sp_e) - gp.quicksum([model._expansion[edge] for edge in sp_e if edge_status[edge] == 0])))   
                             
##############################################
########### END VALID INEQUALITIES ###########
##############################################



##############################################
############ INPUT SETTINGS ##################
##############################################


if len(sys.argv) > 1:
    model_flag = int(sys.argv[1])
    filenames = [sys.argv[2]]
    output_filename = sys.argv[3]
    SEED = int(sys.argv[4])
    pySEED = int(sys.argv[5])
else:
    model_flag = 8
    filename = "pglib-opf-master/pglib_opf_case2383_anor.m"
    filenames = [filename]
    output_filename = ""
    SEED = 0
    pySEED = 0

#ieee 300
if '300' in filenames[0]:
    LP_Length = [1945.15]
    cap = 1

#goc 500
elif '500' in filenames[0]:
    LP_Length = [20692]
    cap = 1.1

#goc 793
elif '793' in filenames[0]:
    LP_Length = [1147.08]
    cap = 1

#############################################
########### EXECUTION BEGINS ################
#############################################


#n = [50]


random.seed(pySEED)
#random.seed(748812)
#cap = .9
demand = .3

#filenames = ["pglib-opf-master/pglib_opf_case500_goc_tep.m"]

#LP_Length = 3015
cap = 1

gen = 1

#for goc 500 tep, used demand = .1, cap = .7
#for 2383, demand = .33, cap = .4



for filename in filenames:
    
    bus_data = np.array([])
    gen_data = np.array([])
    branch_data = np.array([])
    bus_data,gen_data,branch_data = mp.load_data(filename)
    
    graph = nx.Graph()
    mp.encode_graph(graph, bus_data, gen_data, branch_data,demand,cap,gen,1)
    
    
    graph_lr = graph.copy()
    
    
    ##########################################
    ########### Declare TEP Model ############
    ##########################################
    
    mod = gp.Model()
    mod.modelSense = gp.GRB.MINIMIZE
    mod.Params.LogFile = 'anor_test_0.txt'
    mod.Params.MIPGap = .001
    #mod.Params.OutputFlag = 0
    
    
    #Gather line status and cost properties for full graph
    M = 2*.6*max({key: 1/value for  (key,value) in nx.get_edge_attributes(graph,'branch_b').items()}.values())
    edge_status = nx.get_edge_attributes(graph,'branch_status')
    edge_b = nx.get_edge_attributes(graph, 'branch_b')
    expand_lines = [(i,j) for (i,j) in graph.edges if edge_status[i,j] ==0]
    recond_lines = [(i,j) for (i,j) in graph.edges if edge_status[i,j] ==1]
    expand_cost = nx.get_edge_attributes(graph,'branch_cand_cost')
    
    #Add binary decision variables    
    mod._expansion = mod.addVars(expand_lines, vtype = gp.GRB.BINARY, name = 'expansion_status', obj = expand_cost)
    
    #Add continuous variables
    
    ######## Flow Variables #######
    mod._corr_flow = mod.addVars(graph.edges, name='corr_flow', ub = gp.GRB.INFINITY, lb = -gp.GRB.INFINITY)
    

    ######## Generation Variables ########
    mod._gen = mod.addVars(graph.nodes, obj = nx.get_node_attributes(graph, 'gen_cost'), name = 'gen', ub = nx.get_node_attributes(graph, 'gen_Pmax'),
                           lb= nx.get_node_attributes(graph, 'gen_Pmin'))
    
    ######## Angle Variables ########
    mod._bus_angle = mod.addVars(graph.nodes, name = 'bus_angle', ub = 30, lb = -30)
    
    
    ####### Load All Constraints #######
    
    #Load flow constraints
    
    mod.addConstrs((-1/edge_b[i,j] * (mod._bus_angle[i] - mod._bus_angle[j] ) - mod._corr_flow[i,j] + 
                    (1 - mod._expansion[i,j])*M >= 0 for (i,j) in expand_lines),
                   name = "F_eq_pos")
    mod.addConstrs((-1/edge_b[i,j] * (mod._bus_angle[i] - mod._bus_angle[j] ) - mod._corr_flow[i,j] - 
                    (1 - mod._expansion[i,j])*M <= 0 for (i,j) in expand_lines),
                   name = "F_eq_neg")
    mod.addConstrs((-1/edge_b[i,j] * (mod._bus_angle[i] - mod._bus_angle[j] ) ==  mod._corr_flow[i,j] for (i,j) in recond_lines),
                   name = "F_eq")
      
    #Load capacity constraints
    
    mod.addConstrs((mod._corr_flow[i,j] <= graph.edges[i,j]['branch_cap'] * mod._expansion[i,j] for (i,j) in expand_lines), name="F_cap_pos")
    mod.addConstrs((mod._corr_flow[i,j] >= -graph.edges[i,j]['branch_cap'] * mod._expansion[i,j] for (i,j) in expand_lines), name="F_cap_neg")
    
    #Load balance constraints
    
    mod.addConstrs(gp.quicksum([mod._corr_flow[j,i] for j in graph.neighbors(i) if j < i]) - 
                           gp.quicksum([mod._corr_flow[i,j] for j in graph.neighbors(i) if j > i]) + 
                           (mod._gen.get(i) or 0) == nx.get_node_attributes(graph, 'bus_pd')[i] for i in graph.nodes)
    
    #mod.update()
    #mod.optimize()
    #sys.exit()


    ##########################################
    ############# End TEP Model ##############
    ##########################################


    ##########################################
    ########### Declare Linear Model #########
    ##########################################
    
    linear_mod = gp.Model()
    linear_mod.modelSense = gp.GRB.MINIMIZE
    linear_mod.Params.LogFile = 'anor_test_0.txt'
    linear_mod.Params.MIPGap = .001
    #linear_mod.Params.OutputFlag = 0
    
    
    #Gather line status and cost properties for full graph
    M = 2*.6*max({key: 1/value for  (key,value) in nx.get_edge_attributes(graph,'branch_b').items()}.values())
    edge_status = nx.get_edge_attributes(graph,'branch_status')
    edge_b = nx.get_edge_attributes(graph, 'branch_b')
    expand_lines = [(i,j) for (i,j) in graph.edges if edge_status[i,j] ==0]
    recond_lines = [(i,j) for (i,j) in graph.edges if edge_status[i,j] ==1]
    expand_cost = nx.get_edge_attributes(graph,'branch_cand_cost')
    
    #Add binary decision variables    
    linear_mod._expansion = linear_mod.addVars(expand_lines, ub = 1, name = 'expansion_status', obj = expand_cost)
    
    #Add continuous variables
    
    ######## Flow Variables #######
    linear_mod._corr_flow = linear_mod.addVars(graph.edges, name='corr_flow', ub = gp.GRB.INFINITY, lb = -gp.GRB.INFINITY)
    

    ######## Generation Variables ########
    linear_mod._gen = linear_mod.addVars(graph.nodes, obj = nx.get_node_attributes(graph, 'gen_cost'), name = 'gen', ub = nx.get_node_attributes(graph, 'gen_Pmax'),
                           lb= nx.get_node_attributes(graph, 'gen_Pmin'))
    
    ######## Angle Variables ########
    linear_mod._bus_angle = linear_mod.addVars(graph.nodes, name = 'bus_angle', ub = 30, lb = -30)
    
    
    ####### Load All Constraints #######
    
    #Load flow constraints
    
    linear_mod.addConstrs((-1/edge_b[i,j] * (linear_mod._bus_angle[i] - linear_mod._bus_angle[j] ) - linear_mod._corr_flow[i,j] + 
                    (1 - linear_mod._expansion[i,j])*M >= 0 for (i,j) in expand_lines),
                   name = "F_eq_pos")
    linear_mod.addConstrs((-1/edge_b[i,j] * (linear_mod._bus_angle[i] - linear_mod._bus_angle[j] ) - linear_mod._corr_flow[i,j] - 
                    (1 - linear_mod._expansion[i,j])*M <= 0 for (i,j) in expand_lines),
                   name = "F_eq_neg")
    linear_mod.addConstrs((-1/edge_b[i,j] * (linear_mod._bus_angle[i] - linear_mod._bus_angle[j] ) ==  linear_mod._corr_flow[i,j] for (i,j) in recond_lines),
                   name = "F_eq")
      
    #Load capacity constraints
    
    linear_mod.addConstrs((linear_mod._corr_flow[i,j] <= graph.edges[i,j]['branch_cap'] * linear_mod._expansion[i,j] for (i,j) in expand_lines), name="F_cap_pos")
    linear_mod.addConstrs((linear_mod._corr_flow[i,j] >= -graph.edges[i,j]['branch_cap'] * linear_mod._expansion[i,j] for (i,j) in expand_lines), name="F_cap_neg")
    
    #Load balance constraints
    
    linear_mod.addConstrs(gp.quicksum([linear_mod._corr_flow[j,i] for j in graph.neighbors(i) if j < i]) - 
                           gp.quicksum([linear_mod._corr_flow[i,j] for j in graph.neighbors(i) if j > i]) + 
                           (linear_mod._gen.get(i) or 0) == nx.get_node_attributes(graph, 'bus_pd')[i] for i in graph.nodes)
    
    

    ##########################################
    ############# End Linear Model ###########
    ##########################################
    
    ##########################################
    ########### Declare Transport Model ######
    ##########################################
    
    transport_mod = gp.Model()
    transport_mod.modelSense = gp.GRB.MINIMIZE
    transport_mod.Params.LogFile = 'anor_test_0.txt'
    transport_mod.Params.MIPGap = .001
    #transport_mod.Params.OutputFlag = 0
    
    
    #Gather line status and cost properties for full graph
    M = 2*.6*max({key: 1/value for  (key,value) in nx.get_edge_attributes(graph,'branch_b').items()}.values())
    edge_status = nx.get_edge_attributes(graph,'branch_status')
    edge_b = nx.get_edge_attributes(graph, 'branch_b')
    expand_lines = [(i,j) for (i,j) in graph.edges if edge_status[i,j] ==0]
    recond_lines = [(i,j) for (i,j) in graph.edges if edge_status[i,j] ==1]
    expand_cost = nx.get_edge_attributes(graph,'branch_cand_cost')
    
    #Add binary decision variables    
    transport_mod._expansion = transport_mod.addVars(expand_lines, ub = 1, name = 'expansion_status', obj = expand_cost)
    
    #Add continuous variables
    
    ######## Flow Variables #######
    transport_mod._corr_flow = transport_mod.addVars(graph.edges, name='corr_flow', ub = gp.GRB.INFINITY, lb = -gp.GRB.INFINITY)
    

    ######## Generation Variables ########
    transport_mod._gen = transport_mod.addVars(graph.nodes, obj = nx.get_node_attributes(graph, 'gen_cost'), name = 'gen', ub = nx.get_node_attributes(graph, 'gen_Pmax'),
                           lb= nx.get_node_attributes(graph, 'gen_Pmin'))
    
    ######## Angle Variables ########
    transport_mod._bus_angle = transport_mod.addVars(graph.nodes, name = 'bus_angle', ub = 30, lb = -30)
    
    
    ####### Load All Constraints #######

    #Load capacity constraints
    
    transport_mod.addConstrs((transport_mod._corr_flow[i,j] <= graph.edges[i,j]['branch_cap'] * transport_mod._expansion[i,j] for (i,j) in expand_lines), name="F_cap_pos")
    transport_mod.addConstrs((transport_mod._corr_flow[i,j] >= -graph.edges[i,j]['branch_cap'] * transport_mod._expansion[i,j] for (i,j) in expand_lines), name="F_cap_neg")
    
    #Load balance constraints
    
    transport_mod.addConstrs(gp.quicksum([transport_mod._corr_flow[j,i] for j in graph.neighbors(i) if j < i]) - 
                           gp.quicksum([transport_mod._corr_flow[i,j] for j in graph.neighbors(i) if j > i]) + 
                           (transport_mod._gen.get(i) or 0) == nx.get_node_attributes(graph, 'bus_pd')[i] for i in graph.nodes)
    

    ##########################################
    ############# End Transport Model ########
    ##########################################
    
    ##########################################
    ########### Declare Hybrid Model #########
    ##########################################
    
    hybrid_mod = gp.Model()
    hybrid_mod.modelSense = gp.GRB.MINIMIZE
    hybrid_mod.Params.LogFile = 'anor_test_0.txt'
    hybrid_mod.Params.MIPGap = .001
    #hybrid_mod.Params.OutputFlag = 0
    
    
    #Gather line status and cost properties for full graph
    M = 2*.6*max({key: 1/value for  (key,value) in nx.get_edge_attributes(graph,'branch_b').items()}.values())
    edge_status = nx.get_edge_attributes(graph,'branch_status')
    edge_b = nx.get_edge_attributes(graph, 'branch_b')
    expand_lines = [(i,j) for (i,j) in graph.edges if edge_status[i,j] ==0]
    recond_lines = [(i,j) for (i,j) in graph.edges if edge_status[i,j] ==1]
    expand_cost = nx.get_edge_attributes(graph,'branch_cand_cost')
    
    #Add binary decision variables    
    hybrid_mod._expansion = hybrid_mod.addVars(expand_lines, ub = 1, name = 'expansion_status', obj = expand_cost)
    
    #Add continuous variables
    
    ######## Flow Variables #######
    hybrid_mod._corr_flow = hybrid_mod.addVars(graph.edges, name='corr_flow', ub = gp.GRB.INFINITY, lb = -gp.GRB.INFINITY)
    

    ######## Generation Variables ########
    hybrid_mod._gen = hybrid_mod.addVars(graph.nodes, obj = nx.get_node_attributes(graph, 'gen_cost'), name = 'gen', ub = nx.get_node_attributes(graph, 'gen_Pmax'),
                           lb= nx.get_node_attributes(graph, 'gen_Pmin'))
    
    ######## Angle Variables ########
    hybrid_mod._bus_angle = hybrid_mod.addVars(graph.nodes, name = 'bus_angle', ub = 30, lb = -30)
    
    
    ####### Load All Constraints #######
    
    #Load flow constraints

    hybrid_mod.addConstrs((-1/edge_b[i,j] * (hybrid_mod._bus_angle[i] - hybrid_mod._bus_angle[j] ) ==  hybrid_mod._corr_flow[i,j] for (i,j) in recond_lines),
                   name = "F_eq")
      
    #Load capacity constraints
    
    hybrid_mod.addConstrs((hybrid_mod._corr_flow[i,j] <= graph.edges[i,j]['branch_cap'] * hybrid_mod._expansion[i,j] for (i,j) in expand_lines), name="F_cap_pos")
    hybrid_mod.addConstrs((hybrid_mod._corr_flow[i,j] >= -graph.edges[i,j]['branch_cap'] * hybrid_mod._expansion[i,j] for (i,j) in expand_lines), name="F_cap_neg")
    
    #Load balance constraints
    
    hybrid_mod.addConstrs(gp.quicksum([hybrid_mod._corr_flow[j,i] for j in graph.neighbors(i) if j < i]) - 
                           gp.quicksum([hybrid_mod._corr_flow[i,j] for j in graph.neighbors(i) if j > i]) + 
                           (hybrid_mod._gen.get(i) or 0) == nx.get_node_attributes(graph, 'bus_pd')[i] for i in graph.nodes)
    

    ##########################################
    ############# End Hybrid Model ###########
    ##########################################
    
    transport_mod.optimize()
    linear_mod.optimize()
    hybrid_mod.optimize()
    
    
    mod_list = ['lin','hyb']
    pos_flow = fused_graph_pos(graph, mod_list, linear_mod, hybrid_mod, transport_mod)
    neg_flow = fused_graph_neg(graph, mod_list, linear_mod, hybrid_mod, transport_mod)
    
    path_pairs = []
    
    
    for start in pos_flow.nodes():
        for end in pos_flow.nodes():
            if end>start:
                if nx.has_path(pos_flow,  start, end):
                    sp = list(islice(nx.shortest_simple_paths(pos_flow, start, end, weight = 'branch_CR'),200))
                    sp = nx.utils.pairwise(sp)
                    path_pairs += sp
    for start in neg_flow.nodes():
        for end in neg_flow.nodes():
            if end>start:
                if nx.has_path(neg_flow,  start, end):
                    sp = list(islice(nx.shortest_simple_paths(neg_flow, start, end, weight = 'branch_CR'),200))
                    sp = nx.utils.pairwise(sp)
                    path_pairs += sp
    
    for path in path_pairs:
        path_pairs.remove(path)
        path1 = list(nx.utils.pairwise(path[0]))
        path2 = list(nx.utils.pairwise(path[1]))
        mp_length = graph.edge_subgraph(path1).size(weight = 'branch_CR')
        start = path[0][0]
        end = path[-1][1]
        
        sp_e = path2
        sp_length = graph.edge_subgraph(sp_e).size(weight = 'branch_CR')

        if sp_length>mp_length:
            sp_e,path1 = path1,sp_e
            sp_length,mp_length = mp_length,sp_length
            
        for i in range(len(path1)):
            if path1[i][0] > path1[i][1]:
                path1[i] = tuple(reversed(path1[i]))
        for i in range(len(sp_e)):
            if sp_e[i][0] > sp_e[i][1]:
                sp_e[i] = tuple(reversed(sp_e[i]))
                
        edge_status = nx.get_edge_attributes(graph,'branch_status')
        mod.addConstr(mod._bus_angle[start] - mod._bus_angle[end] <= sp_length + 
                 (mp_length -sp_length) * (len(sp_e) - gp.quicksum([mod._expansion[edge] for edge in sp_e if edge_status[edge] == 0])))
        mod.addConstr(mod._bus_angle[end] - mod._bus_angle[start] <= sp_length + 
                 (mp_length -sp_length) * (len(sp_e) - gp.quicksum([mod._expansion[edge] for edge in sp_e if edge_status[edge] == 0])))   
    
    mod.update()
    r = mod.relax()
    r.update()
    r.optimize()
    print(linear_mod.getObjective().getValue())
    print(r.getObjective().getValue())
    sys.exit()
    
    
    
    
    
    
    
    
    
    
    #model_flag = 2
    '''if model_flag >= 2:
        linear_mod.optimize()
        transport_mod.optimize()
        hybrid_mod.optimize()
        mod.Params.Presolve = 0
        mod.Params.lazyConstraints = 1
        mod.optimize()
        sys.exit()'''
        
    if model_flag == 2:    
        mod_list = ['lin']
        pos_flow = fused_graph_pos(graph, mod_list, linear_mod, hybrid_mod, transport_mod)
        neg_flow = fused_graph_neg(graph, mod_list, linear_mod, hybrid_mod, transport_mod)
        
        path_pairs = []
        
        
        for start in pos_flow.nodes():
            for end in pos_flow.nodes():
                if end>start:
                    if nx.has_path(pos_flow,  start, end):
                        sp = list(islice(nx.shortest_simple_paths(pos_flow, start, end, weight = 'branch_CR'),200))
                        sp = nx.utils.pairwise(sp)
                        path_pairs += sp
        for start in neg_flow.nodes():
            for end in neg_flow.nodes():
                if end>start:
                    if nx.has_path(neg_flow,  start, end):
                        sp = list(islice(nx.shortest_simple_paths(neg_flow, start, end, weight = 'branch_CR'),200))
                        sp = nx.utils.pairwise(sp)
                        path_pairs += sp
        
        timing = open('records/' + filename[18:-2] + '_'.join(mod_list) + '.txt', 'a')
        t0 = time.time()
        mod.Params.lazyConstraints = 1
        mod.Params.seed = SEED
        mod.Params.MIPGap = .0001
        mod.optimize(theorem_7_heuristic)
        timing.write('LR VIs: ' + str(time.time()-t0) + ' obj:' + str(mod.objVal) + ' ' + 
                 str(filename) + ' ' + str(SEED) + ' ' + str(pySEED) + ' \n')
        timing.close()
        
    if model_flag == 3:    
        mod_list = ['hyb']
        pos_flow = fused_graph_pos(graph, mod_list, linear_mod, hybrid_mod, transport_mod)
        neg_flow = fused_graph_neg(graph, mod_list, linear_mod, hybrid_mod, transport_mod)
        
        path_pairs = []
        
        
        for start in pos_flow.nodes():
            for end in pos_flow.nodes():
                if end>start:
                    if nx.has_path(pos_flow,  start, end):
                        sp = list(islice(nx.shortest_simple_paths(pos_flow, start, end, weight = 'branch_CR'),200))
                        sp = nx.utils.pairwise(sp)
                        path_pairs += sp
        for start in neg_flow.nodes():
            for end in neg_flow.nodes():
                if end>start:
                    if nx.has_path(neg_flow,  start, end):
                        sp = list(islice(nx.shortest_simple_paths(neg_flow, start, end, weight = 'branch_CR'),200))
                        sp = nx.utils.pairwise(sp)
                        path_pairs += sp
        
        timing = open('records/' + filename[18:-2] + '_'.join(mod_list) + '.txt', 'a')
        t0 = time.time()
        mod.Params.lazyConstraints = 1
        mod.Params.seed = SEED
        mod.Params.MIPGap = .0001
        mod.optimize(theorem_7_heuristic)
        timing.write('HR VIs: ' + str(time.time()-t0) + ' obj:' + str(mod.objVal) + ' ' + 
                 str(filename) + ' ' + str(SEED) + ' ' + str(pySEED) + ' \n')
        timing.close()
            
    
    if model_flag == 4:    
        mod_list = ['trans']
        pos_flow = fused_graph_pos(graph, mod_list, linear_mod, hybrid_mod, transport_mod)
        neg_flow = fused_graph_neg(graph, mod_list, linear_mod, hybrid_mod, transport_mod)
        
        path_pairs = []
        
        
        for start in pos_flow.nodes():
            for end in pos_flow.nodes():
                if end>start:
                    if nx.has_path(pos_flow,  start, end):
                        sp = list(islice(nx.shortest_simple_paths(pos_flow, start, end, weight = 'branch_CR'),200))
                        sp = nx.utils.pairwise(sp)
                        path_pairs += sp
        for start in neg_flow.nodes():
            for end in neg_flow.nodes():
                if end>start:
                    if nx.has_path(neg_flow,  start, end):
                        sp = list(islice(nx.shortest_simple_paths(neg_flow, start, end, weight = 'branch_CR'),200))
                        sp = nx.utils.pairwise(sp)
                        path_pairs += sp
        
        timing = open('records/' + filename[18:-2] + '_'.join(mod_list) + '.txt', 'a')
        t0 = time.time()
        mod.Params.lazyConstraints = 1
        mod.Params.seed = SEED
        mod.Params.MIPGap = .0001
        mod.optimize(theorem_7_heuristic)
        timing.write('TR VIs: ' + str(time.time()-t0) + ' obj:' + str(mod.objVal) + ' ' + 
                 str(filename) + ' ' + str(SEED) + ' ' + str(pySEED) + ' \n')
        timing.close()  
        
    if model_flag == 5:    
        mod_list = ['lin','hyb']
        pos_flow = fused_graph_pos(graph, mod_list, linear_mod, hybrid_mod, transport_mod)
        neg_flow = fused_graph_neg(graph, mod_list, linear_mod, hybrid_mod, transport_mod)
        
        path_pairs = []
        
        
        for start in pos_flow.nodes():
            for end in pos_flow.nodes():
                if end>start:
                    if nx.has_path(pos_flow,  start, end):
                        sp = list(islice(nx.shortest_simple_paths(pos_flow, start, end, weight = 'branch_CR'),200))
                        sp = nx.utils.pairwise(sp)
                        print(sp)
                        path_pairs += sp
        for start in neg_flow.nodes():
            for end in neg_flow.nodes():
                if end>start:
                    if nx.has_path(neg_flow,  start, end):
                        sp = list(islice(nx.shortest_simple_paths(neg_flow, start, end, weight = 'branch_CR'),200))
                        sp = nx.utils.pairwise(sp)
                        path_pairs += sp
        
        timing = open('records/' + filename[18:-2] + '_'.join(mod_list) + '.txt', 'a')
        t0 = time.time()
        mod.Params.lazyConstraints = 1
        mod.Params.seed = SEED
        mod.Params.MIPGap = .0001
        mod.optimize(theorem_7_heuristic)
        timing.write('LR+HR VIs: ' + str(time.time()-t0) + ' obj:' + str(mod.objVal) + ' ' + 
                 str(filename) + ' ' + str(SEED) + ' ' + str(pySEED) + ' \n')
        timing.close()
        
    if model_flag == 6:    
        mod_list = ['lin','trans']
        pos_flow = fused_graph_pos(graph, mod_list, linear_mod, hybrid_mod, transport_mod)
        neg_flow = fused_graph_neg(graph, mod_list, linear_mod, hybrid_mod, transport_mod)
        
        path_pairs = []
        
        
        for start in pos_flow.nodes():
            for end in pos_flow.nodes():
                if end>start:
                    if nx.has_path(pos_flow,  start, end):
                        sp = list(islice(nx.shortest_simple_paths(pos_flow, start, end, weight = 'branch_CR'),200))
                        sp = nx.utils.pairwise(sp)
                        path_pairs += sp
        for start in neg_flow.nodes():
            for end in neg_flow.nodes():
                if end>start:
                    if nx.has_path(neg_flow,  start, end):
                        sp = list(islice(nx.shortest_simple_paths(neg_flow, start, end, weight = 'branch_CR'),200))
                        sp = nx.utils.pairwise(sp)
                        path_pairs += sp
        
        timing = open('records/' + filename[18:-2] + '_'.join(mod_list) + '.txt', 'a')
        t0 = time.time()
        mod.Params.lazyConstraints = 1
        mod.Params.seed = SEED
        mod.Params.MIPGap = .0001
        mod.optimize(theorem_7_heuristic)
        timing.write('LR+tR VIs: ' + str(time.time()-t0) + ' obj:' + str(mod.objVal) + ' ' + 
                 str(filename) + ' ' + str(SEED) + ' ' + str(pySEED) + ' \n')
        timing.close()
        
    if model_flag == 7:    
        mod_list = ['trans','hyb']
        pos_flow = fused_graph_pos(graph, mod_list, linear_mod, hybrid_mod, transport_mod)
        neg_flow = fused_graph_neg(graph, mod_list, linear_mod, hybrid_mod, transport_mod)
        
        path_pairs = []
        
        
        for start in pos_flow.nodes():
            for end in pos_flow.nodes():
                if end>start:
                    if nx.has_path(pos_flow,  start, end):
                        sp = list(islice(nx.shortest_simple_paths(pos_flow, start, end, weight = 'branch_CR'),200))
                        sp = nx.utils.pairwise(sp)
                        path_pairs += sp
        for start in neg_flow.nodes():
            for end in neg_flow.nodes():
                if end>start:
                    if nx.has_path(neg_flow,  start, end):
                        sp = list(islice(nx.shortest_simple_paths(neg_flow, start, end, weight = 'branch_CR'),200))
                        sp = nx.utils.pairwise(sp)
                        path_pairs += sp
        
        timing = open('records/' + filename[18:-2] + '_'.join(mod_list) + '.txt', 'a')
        t0 = time.time()
        mod.Params.lazyConstraints = 1
        mod.Params.seed = SEED
        mod.Params.MIPGap = .0001
        mod.optimize(theorem_7_heuristic)
        timing.write('TR+HR VIs: ' + str(time.time()-t0) + ' obj:' + str(mod.objVal) + ' ' + 
                 str(filename) + ' ' + str(SEED) + ' ' + str(pySEED) + ' \n')
        timing.close()
        
    if model_flag == 8:    
        mod_list = ['lin','hyb','trans']
        pos_flow = fused_graph_pos(graph, mod_list, linear_mod, hybrid_mod, transport_mod)
        neg_flow = fused_graph_neg(graph, mod_list, linear_mod, hybrid_mod, transport_mod)
        
        path_pairs = []
        
        
        for start in pos_flow.nodes():
            for end in pos_flow.nodes():
                if end>start:
                    if nx.has_path(pos_flow,  start, end):
                        sp = list(islice(nx.shortest_simple_paths(pos_flow, start, end, weight = 'branch_CR'),200))
                        sp = nx.utils.pairwise(sp)
                        path_pairs += sp
        for start in neg_flow.nodes():
            for end in neg_flow.nodes():
                if end>start:
                    if nx.has_path(neg_flow,  start, end):
                        sp = list(islice(nx.shortest_simple_paths(neg_flow, start, end, weight = 'branch_CR'),200))
                        sp = nx.utils.pairwise(sp)
                        path_pairs += sp
        
        timing = open('records/' + filename[18:-2] + '_'.join(mod_list) + '.txt', 'a')
        t0 = time.time()
        mod.Params.lazyConstraints = 1
        mod.Params.seed = SEED
        mod.Params.MIPGap = .0001
        mod.optimize(theorem_7_heuristic)
        timing.write('LR+HR+TR VIs: ' + str(time.time()-t0) + ' obj:' + str(mod.objVal) + ' ' + 
                 str(filename) + ' ' + str(SEED) + ' ' + str(pySEED) + ' \n')
        timing.close()
    
    elif model_flag == 0:
        timing = open('records/' + filename[18:-2] + 'orig.txt', 'a')
        t0 = time.time()
        mod.Params.lazyConstraints = 1
        mod.Params.presolve = 0
        mod.Params.seed = SEED
        mod.Params.MIPGap = .0001
        mod.optimize()
        timing.write('orig solve: ' + str(time.time()-t0) + ' obj:' + str(mod.objVal) + ' ' + 
                 str(filename) + ' ' + str(SEED) + ' ' + str(pySEED) + ' \n')
        mod.write(filename[18:-2] + '.sol')
        timing.close()
    elif model_flag == 1:
        timing = open('records/' + filename[18:-2]+ 'cycle_vis.txt', 'a')
        t0 = time.time()
        cycle_paths = []
        cycleBasis = nx.cycle_basis(graph)
        for cycle in cycleBasis:
            for node in cycle[1:-1]:
                path1 = cycle[:cycle.index(node)+1]
                path2 = cycle[cycle.index(node)+1:]
                path2.append(cycle[0])
                cycle_paths.append((path1,path2))
        
        mod.Params.lazyConstraints = 1
        mod.Params.seed = SEED
        mod.Params.MIPGap = .0001
        mod.optimize(theorem_7_cycle_basis)
        timing.write('cycle vis: ' + str(time.time()-t0) + ' obj:' + str(mod.objVal) + ' ' + 
                     str(filename) + ' ' + str(SEED) + ' ' + str(pySEED) + ' \n')
        timing.close()

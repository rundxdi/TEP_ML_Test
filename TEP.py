import mpinput as mp
import networkx as nx
import gurobipy as gp
import numpy as np
from miplearn import Instance
import random
from pyomo.environ import *


class TEPInstance(Instance):

    def __init__(self, filename):
        self.demand_scale = .01
        self.capacity_scale = 1
        self.gen_scale = 1
        self.cost_scale = 1
        self.bus_data = np.array([])
        self.gen_data = np.array([])
        self.branch_data = np.array([])
        self.graph = nx.Graph()
        self.filename = filename
        self.found_violated_lazy_constraints = []
        self.lp_value=0

    def cycle_base(self):
        cycle_paths = []
        cycleBasis = nx.cycle_basis(self.graph)
        for cycle in cycleBasis:
            for node in cycle[1:-1]:
                path1 = cycle[:cycle.index(node) + 1]
                path2 = cycle[cycle.index(node):]
                path2.append(cycle[0])
                cycle_paths.append((path1, path2))

        cycle_cuts = []

        while len(cycle_paths) > 0:
            path = random.choice(cycle_paths)
            cycle_paths.remove(path)
            path1 = list(nx.utils.pairwise(path[0]))
            path2 = list(nx.utils.pairwise(path[1]))
            mp_length = self.graph.edge_subgraph(path1).size(weight='branch_CR')
            start = path[0][0]
            end = path[-1][1]

            sp_e = path2
            sp_length = self.graph.edge_subgraph(sp_e).size(weight='branch_CR')

            if sp_length > mp_length:
                sp_e, path1 = path1, sp_e
                sp_length, mp_length = mp_length, sp_length

            for i in range(len(path1)):
                if path1[i][0] > path1[i][1]:
                    path1[i] = tuple(reversed(path1[i]))
            for i in range(len(sp_e)):
                if sp_e[i][0] > sp_e[i][1]:
                    sp_e[i] = tuple(reversed(sp_e[i]))

            cycle_cuts.append((start, end, sp_length, mp_length, sp_e))

        return cycle_cuts

    '''
    def to_model(self):
        self.bus_data, self.gen_data, self.branch_data = mp.load_data(self.filename)
        mp.encode_graph(self.graph, self.bus_data, self.gen_data, self.branch_data, self.demand_scale, self.capacity_scale, self.gen_scale, self.cost_scale)
        mod = gp.Model()

        # Gather line status and cost properties for full graph
        M = 2 * .6 * max({key: 1 / value for (key, value) in nx.get_edge_attributes(self.graph, 'branch_b').items()}.values())
        edge_status = nx.get_edge_attributes(self.graph, 'branch_status')
        edge_b = nx.get_edge_attributes(self.graph, 'branch_b')
        expand_lines = [(i, j) for (i, j) in self.graph.edges if edge_status[i, j] == 0]
        recond_lines = [(i, j) for (i, j) in self.graph.edges if edge_status[i, j] == 1]
        expand_cost = nx.get_edge_attributes(self.graph, 'branch_cand_cost')

        # Add binary decision variables
        mod._expansion = mod.addVars(expand_lines, vtype=gp.GRB.BINARY, name='expansion_status', obj=expand_cost)

        # Add continuous variables

        ######## Flow Variables #######
        mod._corr_flow = mod.addVars(self.graph.edges, name='corr_flow', ub=gp.GRB.INFINITY, lb=-gp.GRB.INFINITY)

        ######## Generation Variables ########
        mod._gen = mod.addVars(self.graph.nodes, obj=nx.get_node_attributes(self.graph, 'gen_cost'), name='gen', ub=nx.get_node_attributes(self.graph, 'gen_Pmax'),
                               lb=nx.get_node_attributes(self.graph, 'gen_Pmin'))

        ######## Angle Variables ########
        mod._bus_angle = mod.addVars(self.graph.nodes, name='bus_angle', ub=30, lb=-30)

        ####### Load All Constraints #######

        # Load flow constraints

        mod.addConstrs((-1 / edge_b[i, j] * (mod._bus_angle[i] - mod._bus_angle[j]) - mod._corr_flow[i, j] +
                        (1 - mod._expansion[i, j]) * M >= 0 for (i, j) in expand_lines),
                       name="F_eq_pos")
        mod.addConstrs((-1 / edge_b[i, j] * (mod._bus_angle[i] - mod._bus_angle[j]) - mod._corr_flow[i, j] -
                        (1 - mod._expansion[i, j]) * M <= 0 for (i, j) in expand_lines),
                       name="F_eq_neg")
        mod.addConstrs((-1 / edge_b[i, j] * (mod._bus_angle[i] - mod._bus_angle[j]) == mod._corr_flow[i, j] for (i, j) in recond_lines),
                       name="F_eq")

        # Load capacity constraints

        mod.addConstrs((mod._corr_flow[i, j] <= self.graph.edges[i, j]['branch_cap'] * mod._expansion[i, j] for (i, j) in expand_lines), name="F_cap_pos")
        mod.addConstrs((mod._corr_flow[i, j] >= -self.graph.edges[i, j]['branch_cap'] * mod._expansion[i, j] for (i, j) in expand_lines), name="F_cap_neg")

        # Load balance constraints

        mod.addConstrs((gp.quicksum([mod._corr_flow[j, i] for j in self.graph.neighbors(i) if j < i]) -
                        gp.quicksum([mod._corr_flow[i, j] for j in self.graph.neighbors(i) if j > i]) +
                        (mod._gen.get(i) or 0) == nx.get_node_attributes(self.graph, 'bus_pd')[i] for i in self.graph.nodes), name="Flow_Balance")

        self.cycle_cuts = self.cycle_base()

        return mod
    '''



    def to_model(self):
        self.bus_data, self.gen_data, self.branch_data = mp.load_data(self.filename)
        mp.encode_graph(self.graph, self.bus_data, self.gen_data, self.branch_data, self.demand_scale, self.capacity_scale, self.gen_scale, self.cost_scale)
        mod = ConcreteModel()
        mod.graph = self.graph
        # Gather line status and cost properties for full graph
        mod.M = 2 * .6 * max({key: 1 / value for (key, value) in nx.get_edge_attributes(self.graph, 'branch_b').items()}.values())
        mod.edge_status = nx.get_edge_attributes(self.graph, 'branch_status')
        mod.edge_b = nx.get_edge_attributes(self.graph, 'branch_b')
        mod.expand_lines = [(i, j) for (i, j) in self.graph.edges if mod.edge_status[i, j] == 0]
        mod.recond_lines = [(i, j) for (i, j) in self.graph.edges if mod.edge_status[i, j] == 1]
        mod.expand_cost = nx.get_edge_attributes(self.graph, 'branch_cand_cost')
        mod.found_violated_lazy_constraints = ConstraintList()

        # Add binary decision variables
        mod.expansion = Var(mod.expand_lines, within=Binary)

        # Add continuous variables
        def flow_bounds(model, i, j):
            return (-model.graph.edges[i, j]['branch_cap'], model.graph.edges[i, j]['branch_cap'])

        mod.corr_flow = Var(self.graph.edges, within=Reals, bounds=flow_bounds)

        def gen_bounds(model, node):
            return (model.graph.nodes[node]['gen_Pmin'], model.graph.nodes[node]['gen_Pmax'])

        mod.gen = Var(self.graph.nodes, within=NonNegativeReals, bounds=gen_bounds)

        mod.bus_angle = Var(self.graph.nodes, bounds=(-30, 30))

        # Load flow constraints

        def load_cons_1(model,i,j):
            return (-1 / model.edge_b[i, j]) * (model.bus_angle[i] - model.bus_angle[j]) - model.corr_flow[i, j] + (1 - model.expansion[i, j]) * model.M >= 0

        mod.load_cons_pos = Constraint(mod.expand_lines, rule = load_cons_1)

        def load_cons_2(model,i,j):
            return (-1 / model.edge_b[i, j]) * (model.bus_angle[i] - model.bus_angle[j]) - model.corr_flow[i, j] - (1 - model.expansion[i, j]) * model.M <= 0

        mod.load_cons_neg = Constraint(mod.expand_lines, rule = load_cons_2)

        def load_cons_3(model,i,j):
            return (-1/model.edge_b[i,j])*(model.bus_angle[i] - mod.bus_angle[j]) == mod.corr_flow[i,j]

        mod.load_cons_eq = Constraint(mod.recond_lines, rule=load_cons_3)

        # Load Capacity Constraints

        def cap_cons_pos(model, i, j):
            return model.corr_flow[i, j] <= model.graph.edges[i, j]['branch_cap'] * model.expansion[i, j]
        def cap_cons_neg(model, i, j):
            return model.corr_flow[i, j] >= -model.graph.edges[i, j]['branch_cap'] * model.expansion[i, j]

        mod.cap_cons_pos = Constraint(mod.expand_lines, rule = cap_cons_pos)
        mod.cap_cons_neg = Constraint(mod.expand_lines, rule = cap_cons_neg)

        # Load balance constraints

        def balance_rule(model, i):
            return sum([model.corr_flow[j, i] for j in model.graph.neighbors(i) if j < i]) - sum([model.corr_flow[i, j] for j in model.graph.neighbors(i) if j > i]) + model.gen[i] == nx.get_node_attributes(model.graph, 'bus_pd')[i]

        mod.load_balance = Constraint(mod.graph.nodes, rule=balance_rule)

        # Objective Function

        def obj_rule(model):
            return sum(model.expand_cost[i,j]*model.expansion[i,j] for (i,j) in model.expand_lines) + sum(model.graph.nodes[i]['gen_cost']*model.gen[i] for i in model.graph.nodes)

        mod.obj = Objective(rule=obj_rule)

        self.cycle_cuts = self.cycle_base()


        return mod

    def get_instance_features(self):
        return np.array([1])

    def get_variable_features(self, var_name, idx):
        return np.array([1])

    def find_violated_lazy_constraints(self, model):
        edge_status = nx.get_edge_attributes(self.graph, 'branch_status')
        violations = []
        for c in self.cycle_cuts:
            start, end, sp_length, mp_length, sp_e = c
            lhs = model.bus_angle[start].x - model.bus_angle[end].x
            rhs = sp_length + (mp_length - sp_length) * (len(sp_e) - gp.quicksum([model.expansion[edge].x for edge in sp_e if edge_status[edge] == 0]))
            if 1:
                print(c)
                violations.append(c)
                model.found_violated_lazy_constraints.append(c)
        return violations

    def build_lazy_constraint(self, model, idx):
        edge_status = nx.get_edge_attributes(self.graph, 'branch_status')
        start, end, sp_length, mp_length, sp_e = self.cycle_cuts[idx]
        lhs = model.bus_angle[start] - model.bus_angle[end]
        rhs = sp_length + (mp_length - sp_length) * (len(sp_e) - gp.quicksum([model.expansion[edge] for edge in sp_e if edge_status[edge] == 0]))
        #self.found_violated_lazy_constraints.add(abs(lhs) <= rhs)
        return model.found_violated_lazy_constraints.add(abs(lhs) <= rhs)

    def find_violated_user_cuts(self, model):
        return self.find_violated_lazy_constraints(model)

    def build_user_cut(self, model, violation):
        return self.build_lazy_constraint(model, violation)

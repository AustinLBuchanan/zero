import math
import networkx as nx
import gurobipy as gp
from gurobipy import GRB
from coarsen import hop_coarsen, distance_to_vertex_set

# In the spirit of "ILP-based local search for graph partitioning"
#    by A Henzinger, A Noe, C Schulz - Journal of Experimental Algorithmics, 2020
#    https://scholar.google.com/scholar?cluster=9173497710013156715&hl=en&as_sdt=0,37
def mip_local_search(G, initial_plan, h=1, max_iterations=10, preserve_splits=True, time_limit=600, verbose=False):
    
    if verbose:
        print(f"Applying MIP-based local search to improve the MIP warm start (hop parameter: {h})...")
    new_plan = initial_plan.copy()
    new_obj = float('infinity') # arbitrary "big" number
    grb_time = 0
    
    # LOCAL SEARCH
    if verbose:
        print("iter \t\t obj \t\t time")
    for iteration in range(1,max_iterations+1):
        
        old_plan = new_plan.copy()
        old_obj = new_obj
        ( GH, coarsen_partition, coarsened_old_plan ) = hop_coarsen(G, old_plan, h, preserve_splits )
        GH._k = G._k; GH._L = G._L; GH._U = G._U
        
        (coarsened_new_plan, new_obj, runtime) = best_neighboring_plan(GH, coarsened_old_plan, h=h, preserve_splits=preserve_splits, time_limit=time_limit, verbose=verbose)
        
        new_plan = list()
        for coarsened_district in coarsened_new_plan:
            district = list()
            for i in coarsened_district:
                district += coarsen_partition[i]
            new_plan.append(district)
        
        grb_time += runtime
        if verbose:
            print(iteration,'\t','{0:.8f}'.format(new_obj),'\t','{0:.2f}'.format(runtime))
        
        # because of floating point numerical issues, objectives may not be monotone, especially near the end
        if new_obj >= old_obj - 1e-8 or new_plan == old_plan:
            break
    
    return new_plan


# G = input graph (with attributes L, U, k)
# initial_plan = starting point for local search (list of lists)
# h = # hops permitted
#
def best_neighboring_plan(G, initial_plan, h=1, minority=None, preserve_splits=False, time_limit=600, verbose=False):

    k = len(initial_plan)
    
    # create a (sparse) directed "assignment graph" AG 
    # that has vertices V = V(G) + {0, 1, ..., k-1}  (which may overlap)
    # and edge (i,j) if node i can be assigned to district j in {0, 1, ..., k-1}
    nodes = list( set(G.nodes) | set(range(k)) )
    AG = nx.DiGraph()
    AG.add_nodes_from( nodes )
    
    # add edges to AG
    for j in range(k):
        district = initial_plan[j]
        dist = distance_to_vertex_set(G, district, cutoff=h)
        AG.add_edges_from( [ (i,j) for i in dist.keys() ] )
        
    # initialize MIP
    with gp.Env(empty=True) as env:
        env.setParam('OutputFlag', verbose)
        env.start()
        with gp.Model(env=env) as m:
            m.Params.TimeLimit = time_limit
                
            x = m.addVars(AG.edges, vtype=GRB.BINARY)
        
            # assignment constraints
            m.addConstrs( sum(  x[i,j] for j in AG.neighbors(i) ) == 1 for i in G.nodes )
        
            # population balance
            m.addConstrs( sum( G.nodes[i]['TOTPOP'] * x[i,j] for i in AG.predecessors(j) ) >= G._L for j in range(k) )
            m.addConstrs( sum( G.nodes[i]['TOTPOP'] * x[i,j] for i in AG.predecessors(j) ) <= G._U for j in range(k) )
        
            # fix nodes based on preserve_splits status
            if preserve_splits:
                counties = { G.nodes[i]['GEOID20'][0:5] for i in G.nodes }
                support = { c : list() for c in counties }
                
                for j in range(k):
                    for i in initial_plan[j]:
                        c = G.nodes[i]['GEOID20'][0:5]
                        if j not in support[c]:
                            support[c].append(j)
                
                nonsupport = { c : set(range(k)) - set(support[c]) for c in counties }
                
                for i in G.nodes:
                    c = G.nodes[i]['GEOID20'][0:5]
                    for j in nonsupport[c]:
                        if j in AG.neighbors(i):
                            x[i,j].UB = 0
            
            m.update()
            
            # fix nodes that have only one choice
            num_nodes_fixed = 0
            for i in G.nodes:
                choices = [ j for j in AG.neighbors(i) if x[i,j].UB > 0.5 ]
                if len(choices) == 1:
                    j = choices[0]
                    x[i,j].LB = 1
                    num_nodes_fixed += 1
                    
            m.update()
            
            if verbose:
                print("this many nodes fixed:",num_nodes_fixed)
                print("out of a total:",G.number_of_nodes())
                print("remaining nodes:",G.number_of_nodes()-num_nodes_fixed)

            DG = nx.DiGraph(G)
            
            # minimize the number of people who change districts
            for i,j in AG.edges:
                x[i,j].obj = G.nodes[i]['TOTPOP']
            for j in range(k):
                for i in initial_plan[j]:
                    x[i,j].obj = 0
            
            ###########################################
            # ADD CONSTRAINT: NUM CUT EDGES <= PREV CUT EDGES
            ###########################################
            
            # Interpretation: y[u,v,j]=1  <=> directed edge (u,v) is cut because u->j but v!->j 
            # Implementation: we create y[u,v,j] variable 
            #          <=> 
            # (u,v) is directed edge and x[u,j] was created (i.e., if j in AG.neighbors(u) )
            #
            triples = list()
            
            for u,v in DG.edges:
                triples += [ (u,v,j) for j in AG.neighbors(u) ]
            
            if verbose:
                print("this many triples (y cut edge vars) =",len(triples))
            y = m.addVars(triples, vtype=GRB.BINARY)
            
            # add cut edge constraints
            for u,v in DG.edges:
                for j in AG.neighbors(u): # x[u,j] exists
                    if j in AG.neighbors(v): # x[v,j] exists
                        m.addConstr( x[u,j] - x[v,j] <= y[u,v,j] )
                    else: # x[v,j] does not exist (assumed to be zero)
                        m.addConstr( x[u,j] == y[u,v,j] )
                        
            # fix y vars when safe
            for u,v in DG.edges:
                for j in AG.neighbors(u):
                    if j in AG.neighbors(v): # x[v,j] exists
                        # if x[u,j]=0 or x[v,j]=1 then y[u,v,j]=0
                        if x[u,j].UB < 0.5 or x[v,j].LB > 0.5:
                            y[u,v,j].UB = 0
                        # if x[u,j]=1 and x[v,j]=0 then y[u,v,j]=1
                        if x[u,j].LB > 0.5 and x[v,j].UB < 0.5:
                            y[u,v,j].LB = 1
                    else: # x[v,j] does not exist
                        if x[u,j].UB < 0.5:
                            y[u,v,j].UB = 0
                        if x[u,j].LB > 0.5:
                            y[u,v,j].LB = 1
                        
            # the following undirected cut edge vars are not strictly necessary, but seem to help computationally
            is_cut = m.addVars( G.edges, vtype=GRB.BINARY )
            m.addConstrs( is_cut[u,v] == sum( y[u,v,j] for j in AG.neighbors(u) ) for u,v in G.edges )
            m.addConstrs( is_cut[u,v] == sum( y[v,u,j] for j in AG.neighbors(v) ) for u,v in G.edges )
        
            # Add cut edges constraint
            a = { i : j for j in range(k) for i in initial_plan[j] }
            prev_cut_edges = sum( G.edges[i,j]['weight'] for i,j in G.edges if a[i] != a[j] )
            m.addConstr( sum( G.edges[e]['weight'] * is_cut[e] for e in G.edges ) <= prev_cut_edges )
            
            # add contiguity constraints via callback
            m.Params.LazyConstraints = 1
            m._x = x
            m._k = k
            m._DG = DG
            m._AG = AG
            m._callback = AG_contiguity_callback
            
            # solve MIP 
            m.Params.IntFeasTol = 1e-7
            m.Params.FeasibilityTol = 1e-7
            m.optimize(m._callback)
            
            if m.solCount == 0:
                return (initial_plan, float('inf'), m.runtime) 
                                    
            plan = [ [ i for i in AG.predecessors(j) if m._x[i,j].x > 0.5 ] for j in range(k) ]
            return (plan, m.objVal, m.runtime)


def find_fischetti_separator(DG, component, b):
    neighbors_component = { i : False for i in DG.nodes }
    for i in nx.node_boundary(DG, component, None):
        neighbors_component[i] = True
    
    visited = { i : False for i in DG.nodes }
    child = [ b ]
    visited[b] = True
    
    while child:
        parent = child
        child = list()
        for i in parent:
            if not neighbors_component[i]:
                for j in DG.neighbors(i):
                    if not visited[j]:
                        child.append(j)
                        visited[j] = True
    
    C = [ i for i in DG.nodes if neighbors_component[i] and visited[i] ]
    return C


def AG_contiguity_callback(m, where):

    # check if LP relaxation at this branch-and-bound node has an integer solution
    if where != GRB.Callback.MIPSOL: 
        return
        
    # retrieve the LP solution
    xval = m.cbGetSolution(m._x)
    DG = m._DG
    AG = m._AG

    # check if each district is connected
    for j in range(m._k):

        # which nodes are selected?
        S = list( { i for i in AG.predecessors(j) if xval[i,j] > 0.5 } )
        if len(S)==1 or nx.is_strongly_connected( DG.subgraph(S) ):
            continue

        # for smaller connected components, find some vertex 'a' and add cut.
        max_comp = None
        
        for comp in sorted(nx.strongly_connected_components( DG.subgraph(S) ), key=len, reverse=True):
            if max_comp == None:
                max_comp = comp
                b = nx.utils.arbitrary_element( comp )
                continue

            # pick some vertex from this component
            a = nx.utils.arbitrary_element( comp )

            # add a,b-separator inequality
            C = find_fischetti_separator(DG, comp, b)
            m.cbLazy( m._x[a,j] + m._x[b,j] <= 1 + sum( m._x[c,j] for c in C if c in AG.predecessors(j)) )
            
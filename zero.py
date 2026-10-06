from argparse import ArgumentParser
from utils import read_graph_from_json, read_plan_from_csv, write_plan_to_csv
from cleanup import mip_local_search
from math import floor, ceil

def main():
    parser = ArgumentParser()
    parser.add_argument("--graph", required=True)
    parser.add_argument("--inplan", required=True)
    parser.add_argument("--outplan", required=True)
    parser.add_argument("--deviation", default=1, type=float)
    parser.add_argument("--timelimit", default=60, type=float)
    parser.add_argument("--iterations", default=20, type=int)
    parser.add_argument("--hops", default=2, type=int)
    parser.add_argument("--nopreserve", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    
    args = parser.parse_args()
    print(args)

    # read input graph
    G = read_graph_from_json(args.graph)
    print(f"Input graph has {G.number_of_nodes()} vertices and {G.number_of_edges()} edges.")
    geoid_to_node = { G.nodes[i]['GEOID20'] : i for i in G.nodes }

    # read input plan
    geoid_to_district = read_plan_from_csv(args.inplan)
    G._k = 1 + max( geoid_to_district.values() )
    plan = [ list() for _ in range(G._k) ]
    for g, d in geoid_to_district.items():
        v = geoid_to_node[g]
        plan[d].append(v)

    def pops(plan): return [ sum(G.nodes[i]['TOTPOP'] for i in district) for district in plan ]
    ideal = sum(pops(plan)) / G._k
    def deviation(plan): return max( abs(ideal-pop) for pop in pops(plan) )
    def cut_edges(plan): 
        a = { i : j for j in range(len(plan)) for i in plan[j] }
        return sum( 1 for i,j in G.edges if a[i]!=a[j] )

    # get (more?) compact plan; reduce deviation by half; repeat as needed/able
    for iteration in range(1+args.iterations):
        if iteration > 0:
            plan = mip_local_search(G, plan, h=args.hops, max_iterations=3, preserve_splits=not args.nopreserve, time_limit=args.timelimit, verbose=args.verbose)
            if deviation(plan) > next_dev:
                print("***Got stuck. Exiting now.***")
                print("(Maybe try adjusting the options: preserve_splits, hops, time_limit?)")
                print("(Or, set verbose=True to diagnose where the MIP goes wrong?)")
                break
        dev = deviation(plan)
        print(f"Plan {iteration} has {cut_edges(plan)} cut edges and deviation of +/- {dev:.2f}.")
        if dev <= args.deviation + 1e-6:
            print("-> Plan meets target deviation! Exiting.")
            break
        if iteration < args.iterations:
            G._L = ceil(  min( ideal-args.deviation, ideal-dev/2 ) )
            G._U = floor( max( ideal+args.deviation, ideal+dev/2 ) )
            next_dev=max( G._U-ideal, ideal-G._L )
            print(f"---Imposing L={G._L} and U={G._U}, which would yield deviation <= {next_dev:.2f}.")

    print(f"Writing plan with {cut_edges(plan)} cut edges and deviation of +/- {dev}.")
    geoid_to_district = { G.nodes[i]['GEOID20'] : j for j in range(len(plan)) for i in plan[j] }
    write_plan_to_csv(geoid_to_district, args.outplan)
    
if __name__ == "__main__":
    main()
    
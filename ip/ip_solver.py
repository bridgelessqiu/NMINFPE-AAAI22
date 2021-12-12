import algo
import numpy as np
import scipy as sp
import networkx as nx
import sys

if __name__ == "__main__":
    # --------------------------- #
    #     Command line inputs     #
    # --------------------------- #
    network_name = str(sys.argv[1])

    exp_type = str(sys.argv[2]) # random, uniform

    if exp_type == "uniform":
        k = int(sys.argv[3]) # k-uniform threshold
    
    # ---------------------- #
    #      Path to data      #
    # ---------------------- #
    # The path to the network file
    edgelist_path = "../networks/real/" + network_name + "/" + network_name + ".edges"

    print("Network: {}".format(network_name))
    
    if exp_type == "random":
        threshold_path = "../networks/real/" + network_name + "/" + network_name + "_random_thresh.txt"
        result_file = "results/" + exp_type + "_threshold/" + network_name + ".txt"
    elif exp_type == "uniform":
        threshold_path = "../networks/real/" + network_name + "/" + network_name + "_" + str(k) + "_uniform_thresh.txt"
        result_file = "results/" + exp_type + "_threshold/" + network_name + "_uniform_" + str(k) + ".txt"
    else:
        print("unknown type")
    
    # ------------------------------- #
    #        Network preporcessing    #
    # ------------------------------- #
    # Read in the network
    G = nx.read_edgelist(edgelist_path)
    n = G.number_of_nodes() 
    # the ordering of the node in the matrix (THIS IS VERY VERY IMPORTNAT)
    node_order = [str(i) for i in range(n)]
 
    # The adjacency matrix
    A = nx.to_scipy_sparse_matrix(G, nodelist = node_order)
    A.setdiag(1) # Set all diagonal entries to 1, very important since we consider the closed neighborhood

    # double check the ordering
    for u in range(n):
        if(G.degree(str(u)) + 1 != A[u].count_nonzero()):
            print("wrong")
            print(G.degree(str(u)), A[u].count_nonzero())
    
    # -------------------------- #
    #     Extract thresholds     #
    # -------------------------- #
    threshold_file_obj = open(threshold_path, 'r')
    
    tau = np.array([0] * n)
    for line in threshold_file_obj:
        u = int(line.split(' ')[0])
        t = int(line.split(' ')[1])
        tau[u] = t

    threshold_file_obj.close()

    # ------------------------ #
    #       Run the solver     #
    # ------------------------ #
    x = algo.ip_nmin_fpe(A, tau)
    nmin_fpe = [u for u, s in enumerate(x) if s != 0]
    
    f = open(result_file, 'w')
    f.write(str(len(nmin_fpe)) + '\n')
    f.close()

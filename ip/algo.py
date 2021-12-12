import numpy as np 
import scipy as sp 
from scipy.sparse import diags
import gurobipy as gp
from gurobipy import GRB

def ip_nmin_fpe(A, tau):
    """
    Decription
    ----------
    Solve the NMIN-FPE as an integer program. Formulation see paper.

    Input
    -----
    A: scipy sparse matrix
       The adjacency matrix of the underlying graph G. MUST HAVE ALL ONES ON DIAGOANL

    tau: numpy array
       The threshold vector

    Output
    ------
    C: numpy binary array
       The fixed point configuration
    """

    A.setdiag(1) # VERY IMPORTANT! CLOSED NEIGHBORHOOD
    
    n = np.shape(A)[0] # number of nodes

    delta = n + 1 # constant that is greater than the maximum degree of 

    one = np.ones(n) # vector of all ones

    zero = np.zeros(n)

    # -------------------------------- #
    #      The model and variables     #
    # -------------------------------- #
    M = gp.Model("NMIN-FPE")

    # Variables: {0, 1}
    x = M.addMVar(shape = n, vtype = GRB.BINARY, name = "x")

    # ------------------ #
    #      Objective     #
    # ------------------ #
    M.setObjective(one @ x, GRB.MINIMIZE)

    # ------------------- #
    #     Constraints     #
    # ------------------- #
    # Constraint 1: \tau_v * x_v <= \sum_{each closed neighbor u of v} x_u
    T = diags(tau, shape=(n, n))

    M.addConstr(T @ x <= A @ x, name = "c1")

    # Constraint 2: one * x >= 1
    M.addConstr(one @ x >= 1, name = "c2")

    # Constraint 3: delta * x_v + \tau_v > \sum_{each closed neighbor u of v} x_u
    M.addConstr(delta * x + tau >= A @ x + one, name = "c3")

    # -------------------- #
    #      Time limit      #
    # -------------------- #
    #M.Params.timeLimit = 100.0

    # -------------------- #
    #     Optimization     #
    # -------------------- #
    M.optimize()

    # ---------------- #
    #      Output      #
    # ---------------- #
    print("The Hamming weight of a nontrivial minimum fixed point: {}".format(M.objVal))
    return(x.X)


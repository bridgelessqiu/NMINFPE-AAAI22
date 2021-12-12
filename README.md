# NMIN-FPE
Paper: *Finding Nontrivial Minimum Fixed Points in Discrete Dynamical Systems:Complexity, Special Case Algorithms and Heuristics*
**Full version** (contains complete proofs of theorems): [https://github.com/bridgelessqiu/NMIN-FPE/blob/main/NMIN_FPE_Full_version.pdf]

### Directory layout

	├──  ip/
	|   |_ algo.py
	|   |_ ip_solver.py
	|   |_ results/
	|         |_ random_threshold/
	|         |_ uniform_threshold/
	|
	|_ heuristics/
	|      |_ syds.h
	|      |_ heuristic.h
	|      |_ heuristic.cpp
	|      |_ other_function.h
	|      |_ other_function.cpp
	|      |_ main.cpp
	|      |_ results/
	|            |_ random_thresh/
	|            |_ uniform_thresh/
	|  
	|_ networks/
	      |_ real
		  |_ google_plus/ (a representative network)
		  |_ arena/ (a representative network)


### Directory Overview 

  **ip/**: consists of the implementations of ILP formulation of NMIN-FPE

  **heuristics/**: consists of implementations of the Greedy family

  **network/**: consists of the testing case (random, uniform) for the representative networks google_plus and arena


### Instruction on running the heuristics    

Required software: C++11

Steps:
1. `cd heuristics/`
2. To run the heuristic on google+ network (20,000+ vertices) under random threshold:
   - GreedyThresh: 
	`./greedy google_plus 1 random`

   - GreedyNP: 
	`./greedy google_plus 2 random`

   - GreedySub: 
	`./greedy google_plus 3 random`

   - GreedyFull: 
	`./greedy google_plus 4 random`

3.  To run the heuristic on arena network (10,000+ vertices) under random threshold:
   - GreedyThresh:
   	`./greedy arena 1 random`

   - GreedyNP:
   	`./greedy arena 2 random`

   - GreedySub:
   	`./greedy arena 3 random`

   -GreedyFull:
   	`./greedy arena 4 random`

4. To run the heuristic on google+ network (20,000+ vertices) under uniform threshold with \tau = 8:
   - GreedyThresh: 
	`./greedy google_plus 1 uniform 8`

   - GreedyNP:
	`./greedy google_plus 2 uniform 8`

   - GreedySub:
	`./greedy google_plus 3 uniform 8`

   - GreedyFull
	`./greedy google_plus 4 uniform 8`

5. To run the heuristic on arena network (10,000+ vertices) under uniform threshold with \tau = 8:
   - GreedyThresh:
   	`./greedy arena 1 uniform 8`

   - GreedyNP:
   	`./greedy arena 2 uniform 8`

   - GreedySub:
   	`./greedy arena 3 uniform 8`

   - GreedyFull:
   	`./greedy arena 4 uniform 8`

If the ./greedy is corrputed, please use the following commnad to manually complie the code:
	`g++ -std=c++11 -O3 heuristic.cpp other_function.cpp main.cpp -o greedy`


### Instruction on running the ILP solver 

Required software: (1) python 3.6/7/8; (2) gurobi 9.1.1

Steps:
1. `cd ip/`
3. `Load gurobi`
3. `export PYTHONPATH=$GUROBI_HOME/lib/python3.6_utf32` (The command is not needed if you use python3.7)
4. To run the solve on google+ network under random threshold:
   	`python3 ip_solver.py google_plus random`
5. To run the solver on google+ network under uniform threshold with \tau = 8:
   	`python3 ip_solver.py google_plus uniform 8`

# Complex6GNFPlacement
A Complex-based 6G Network Functions Placement based on multiplex networks and evolutionary game theory.
The “Complex6GNF” Algorithm is based on a complex-based approach aimed at introducing a novel 6G NF Placement taking into account multiplex representation and evolutionary game theory. 
We provide an algorithm where the following conditions are: 

-It provides a placement of 6G-RU, 6G-DU, 6G-CU, and 6G-UPF functions such that each 6G-RU can be connected to a 6G-UPF with a latency lower than a threshold value τ.
-It minimizes the interference between 6G-RUs not connected to the same 6G-DU.  
-It minimizes the number of active nodes, i.e. nodes hosting virtual network functions. As 6G-RUs are physical components and their placement is given, we can achieve this objective by acting on the placement of 6G-DUs, 6G-CUs, and 6G-UPFs

Given a large-scale topology of hierarchical multistage metro-aggregation represented as graph G(V,E) with V the set of network nodes and E the set of optical fiber links, the algorithm “Complex6GNF” requires also an interference matrix I, delay matrix L, and path delay matrix D. Our model assumes that the 6GRUs are distributed over the selected area.
The transmitted power of each 6G-RU such that the coverage radius of each 6G-RU is set at half the distance from the closest 6G-RU.

In this way it is possible building the weighted multiplex network representation and running the algorithm to get the 6G NF placement. 

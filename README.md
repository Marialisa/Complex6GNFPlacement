# Complex-based 6G Network Function Placement
## Description 

Towards 6G, a key challenge lies in the placement of virtual network functions on physical resources. This becomes complex due to the dynamic nature of mobile environments, making the design a major point of research. We propose a framework that sees this challenge as a complex and dynamic collective process, presenting a novel perspective which encompasses transport network and wireless segment aspects. The framework is built around an analytical modeling and algorithmic tools that rely on complex systems' paradigm as multiplex networks and evolutionary game theory. The multiplex network enables capturing the layered and heterogeneous nature of the environment. Evolutionary game theory models the dynamical behavior of the system as a collective social process, where each decision on functions influences the overall outcome. Our model allows us to achieve a placement scheme that optimizes 6G functions deployment and minimizes the number of active computational nodes. Compared to traditional transport network centric approach, it effectively reduces interference, ensuring the network's effective operation and performance. Results show the efficacy of the strategy, enabling the dynamic distribution of functions as the outcome of a social dilemma, and highlight the potential applicability of this approach to tackle the network function placement problem in 6G networks.

The *“Complex6GNF” Algorithm* is based on a complex-based approach aimed at introducing the novel 6G NF Placement taking into account multiplex representation and evolutionary game theory:

-It provides a placement of 6G-RU, 6G-DU, 6G-CU, and 6G-UPF functions such that each 6G-RU can be connected to a 6G-UPF with a latency lower than a threshold value τ.
-It minimizes the interference between 6G-RUs not connected to the same 6G-DU.  
-It minimizes the number of active nodes, i.e. nodes hosting virtual network functions. As 6G-RUs are physical components and their placement is given, we can achieve this objective by acting on the placement of 6G-DUs, 6G-CUs, and 6G-UPFs

Given a large-scale topology of hierarchical multistage metro-aggregation represented as graph G(V,E) with V the set of network nodes and E the set of optical fiber links, the algorithm “Complex6GNF” requires also an interference matrix I, delay matrix L, and path delay matrix D as input dataset. 
Our model assumes that the 6GRUs are distributed over the selected area.
The transmitted power of each 6G-RU such that the coverage radius of each 6G-RU is set at half the distance from the closest 6G-RU.
In this way it is possible building the weighted multiplex network representation and running the algorithm to get the 6G NF placement. 

## Usage 

### Installation 
#### Install R:

If you haven't already installed R on your system, you can download it from the R Project website https://www.r-project.org/ 

#### Install RStudio:

RStudio is an integrated development environment (IDE) for R. You can download it from the RStudio website https://posit.co/

#### Clone or Download the Repository:

Clone the repository to your local machine using Git:

git clone https://github.com/Marialisa/Complex6GNFPlacement.git

Alternatively, you can download the repository as a ZIP file and extract it to your desired location.

#### Open the Project in RStudio:

Launch RStudio and navigate to the directory where you cloned or extracted the repository. Open the .Rproj file to load the project in RStudio.

#### Run the Project:

- Once you've installed the dependencies, you can run the project by executing the main script *C6GNFP.R* or following any specific instructions provided in the project's README or documentation.
- The script *C6GNFP.R* is referred to a network with a population size N=500.
- Before you run the project you need to check the input data for the case study shared in github in the .rar *input data csv (N=500)*, and moreover you need to fix the *your_local_path* as indicated in the script, to run the algorithm with the input data. 

### Citation

If you use this code in your research, please cite the following papers:

#### Scatà Marialisa, Aurelio La Corte, Andrea Marotta, Fabio Graziosi, Dajana Cassioli. "A 6G Function Placement Framework Based on Complex Networks and Evolutionary Dynamics." 2nd International Conference on 6G Networking (6GNet), 2023.

#### Bibtex Citation :
@inproceedings{scata20236g,
  title={A 6G Function Placement Framework Based on Complex Networks and Evolutionary Dynamics},
  author={Scat{\`a}, M and La Corte, A and Marotta, A and Graziosi, F and Cassioli, D},
  booktitle={2023 2nd International Conference on 6G Networking (6GNet)},
  pages={1--4},
  year={2023},
  organization={IEEE}
}

#### Scatà Marialisa, Aurelio La Corte, Andrea Marotta, Fabio Graziosi, Dajana Cassioli. "A Complex Network and Evolutionary Game Theory Framework for 6G Function Placement." IEEE Open Journal of the Communications Society, 2024.(*currently in revision*)



## Extensions
In addition to the original algorithm, this repository, in future reserves the options to include other extensions developed to enhance its functionality or address specific use cases. These extensions will be provided within the repository. Each extension may will have its own documentation and usage instructions.

## License
The code in this repository is licensed under the MIT License. Please see the LICENSE file for more information.

## Contributions
Contributions to this repository, including bug fixes, improvements, and new extensions, are welcome. 

# cross-impact
Matlab implementation of a cross-impact method for generating Bayes-networks based on expert judgments

The methodology is detailed in this paper https://doi.org/10.1002/ffo2.165

The three files ending in network.m are runnable scripts that correspond to the examples in the paper. The energy_network.m file is for the simplified example found in section 3 of the paper, the printing_network.m file is for the case study in section 4, and the tree_netwrok.m file is for the tree-structured network in appendix B.

The code was developed only for this research paper and might not be trivially applicable to other contexts. The comments in the code are in Finnish, and thus not that easily understood. If you want to do something with this on your own, I recommend using the energy_network.m as basis. The network skeleton and the inputs are found at the start of the file before the first for-loop. By default the script does not really output anything interesting. After running it, you will have the joint probability distribution for your uncertainty factors in the jointpd variable. The conditional probability distributions for individual uncertainty factors in the network can be found in the nodes struct. By default it is commented out, but the command genie_parser(nodes,states,"tiedoston_nimi.xdsl"); can be used to automatically generate a Bayesian network file that can be opened using the GeNIe software, that can be downloaded from https://www.bayesfusion.com/genie/

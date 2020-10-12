# STDF-network
Compte et al. (2003) network with a short-term plasticity mechanism (both depression and facilitation)

The main function is called Network-STDAll.m



BuidNetwork320.mat contains two different tables:
  - ExcInh which is a 1x320 vector such that for position i,
      i = 1 if the neuron corresponding to that position is an inhibitory neuron, and 
      i = 0 if the neuron is excitatory.
  - P which is a 320x320 table such that for position (i,j),
      (i,j) = 1 if neurons i and j are connected,  
      (i,j) = 0 otherwise.

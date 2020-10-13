# STDF-network
Compte et al. (2003) network with a short-term plasticity mechanism (both depression and facilitation)


* The main function is called Network-STDAll.m

* Inside the main function, there are some parameters that can be switched according to the preferences of the user to run depression, facilitation, both and also include bicuculline factors or inject some pulse of current. They are:

Dnumber=1.0;                   <-- Depression factor. It should be a value lying in [0,1]. 1 means no depression. 
Fnumber=0.0;                   <-- Facilitation factor. It should be a value lying in [0,1]. 0 means no facilitation. 
BIC=0;                         <-- Bicuculline value. 0 means no Bicuculline
Sm='000';                      <-- Only for saving purposes. Write de number of Depression, of facilitation or bicuculline
iduration=0;                   <-- Duration of the stimulus
AmplitudCurrent=0;             <-- Amplitud of the stimulus
numInjNeurons=0;               <-- Number of neurons we want to inject

caseDF=0;                       <-- Only for saving purposes. State it as 1 for depression, 0 for facilitation, 10 for Depression+Facilitation
DF='F';                         <-- Only for saving purposes. Write:
                                      D for depression
                                      F for facilitation
                                      D-F for depression-facilitation
                                      BIC for Bicuculline
                                      Iapp for applied current


* Rest of functions are called from the main one such that:

- NetworkField.m
  Contains the differential equations that model de network
- Prelease.m
  Update the probability of release of each neuron at each time-step
  
- rk45Network.m
  Is the integrator to solve the differential equations that model the network.



* BuidNetwork320.mat contains two different tables:
  - ExcInh which is a 1x320 vector such that for position i,
      i = 1 if the neuron corresponding to that position is an inhibitory neuron, and 
      i = 0 if the neuron is excitatory.
  - P which is a 320x320 table such that for position (i,j),
      (i,j) = 1 if neurons i and j are connected,  
      (i,j) = 0 otherwise.
      
This table is called from the main function to consider a specific network configuration. 


The outputs of the main function are

(1) Recursive files which names start by solRK_ containing the solution of all the differential equations per each neuron
(2) A file which name starts by voltage_ containing the time course membrane potential of each neuron in the network
(3) A file which name starts by oscillation320_ containing the onset time of spikes for each of the neuron in the network
(4) 3 different files which names start by gAMPA_ , gGABA_ , gNMDA_ containing the time course of the AMPA, GABA and NMDA conductances impinding into the different neurons.


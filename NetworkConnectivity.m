function [P,EI]=NetworkConnectivity(nNeuron,ProbInh,nsynapses,xDist,sigma_exc,sigma_inh)

P=zeros(nNeuron,nNeuron);
EI=zeros(1,nNeuron);


% 1. Decide which neurons are Inhibitoris

n_inh=0;
while n_inh<round(nNeuron*ProbInh)
%     % Using normal distribution
%     EIrand=normrnd(0,1);
%     pos=fix(EIrand*nNeuron);
%     if (EI(pos)==0)
%         EI(pos)=1;
%         n_inh=n_inh+1;
%     end
    pos=random('unid',nNeuron);
    if EI(pos)==0  
        EI(pos)=1;
        n_inh=n_inh+1;
    end
end


% 2. Connections between nsynapses neurons

hN=xDist/nNeuron;           % distance between two cosecutive neurons
for preS=1:nNeuron          % fixed a neuron, we want to connect it to other
    nsyn=0;
    while nsyn<nsynapses    % we are looking for the postsynaptic neurons
        postS=random('unid',nNeuron);
        if (postS~=preS)&&(P(preS,postS)==0)
            dist_PreSPostS=abs(preS-postS)*hN;
            xi=normrnd(0,1);
            % excitatory or inhibitory synapse
            if EI(preS)==1 % the synapse is inhibitory
                parGauss=-(dist_PreSPostS*dist_PreSPostS)/(2*sigma_inh*sigma_inh); 
                gauss=exp(parGauss);
            else % the synapse is excitatory
                parGauss=-(dist_PreSPostS*dist_PreSPostS)/(2*sigma_exc*sigma_exc); 
                gauss=exp(parGauss);
            end
            % connect or not connect
            if abs(gauss)>xi    % connect
                P(preS,postS)=1;
                nsyn=nsyn+1;
            end   
        end
    end
end










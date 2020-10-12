function [pRel,pRel_STF] = rk45Network (RHS, t0, x0, tf, N , I, EI, nvar, neq, Prob, pRel, pRel_STF, Medir, potasio, npreinj, leak1, leak2, gsdrnd, sumes,Vthre, p0_AMPA,p0_NMDA,p0_GABA,factor_relAMPA,factor_relNMDA,factor_relGABA,p0_stfAMPA,p0_stfNMDA,p0_stfGABA,factor_stfAMPA,factor_stfNMDA,factor_stfGABA,fD_AMPA,fD_NMDA,fD_GABA,fF_AMPA, fF_NMDA,fF_GABA,nNeurons,Sm,tON,tOFF,DF,BIC)

longData=1000;
neqn=length(x0);
ti=zeros(longData+1);
% pRelTime=zeros(1,longData+1);
% pRel_STFTime=zeros(1,longData+1);
ti(1) = t0;
wi(1:neqn, 1) = x0';
pRelTime=zeros(longData+1,3*nNeurons);
pRelTime(1,1:end)=pRel;

pRel_STFTime=zeros(longData+1,3*nNeurons);
pRel_STFTime(1,1:end)=pRel_STF;


%%% CREATE OUTPUT FILES %%%

i = 2;
saveFile=2;
numFile=1;
if(N <= 0)
    disp( 'N must be positive and different to 0' );
    return;
else
    h = (tf - t0)/N;
end
while(t0 < tf)
   
    % RK45-Fiel integrator
    pre = x0;
    k1 = h * feval(RHS, t0, x0, I, EI, nvar, neq, Prob, pRel, pRel_STF, Medir, potasio, npreinj, leak1, leak2, gsdrnd, sumes,tON,tOFF,BIC);
    k2 = h * feval(RHS, t0 + h/2, x0 + k1/2, I, EI, nvar, neq, Prob, pRel, pRel_STF, Medir, potasio, npreinj, leak1, leak2, gsdrnd, sumes,tON,tOFF,BIC);
    k3 = h * feval(RHS, t0 + h/2, x0 + k2/2, I, EI, nvar, neq, Prob, pRel, pRel_STF, Medir, potasio, npreinj, leak1, leak2, gsdrnd, sumes,tON,tOFF,BIC);
    k4 = h * feval(RHS, t0 + h, x0 + k3, I, EI, nvar, neq, Prob, pRel, pRel_STF, Medir, potasio, npreinj, leak1, leak2, gsdrnd, sumes,tON,tOFF,BIC);
    

    x0 = x0 + (k1 + 2*k2 + 2*k3 + k4)/6;
    t0 = t0 + h;

    ti(i) = t0;
    wi(1:neqn, i) = x0';
	i = i + 1;
    
      
    %pRel changes
    
    for j=1:nNeurons
        [pRel,pRel_STF] = Prelease(j,neq,nNeurons,EI,pRel,pRel_STF, x0, pre, Vthre, p0_AMPA,p0_NMDA,p0_GABA,factor_relAMPA,factor_relNMDA,factor_relGABA,p0_stfAMPA,p0_stfNMDA,p0_stfGABA,factor_stfAMPA,factor_stfNMDA,factor_stfGABA,fD_AMPA,fD_NMDA,fD_GABA,fF_AMPA, fF_NMDA,fF_GABA);
    end
    
    % Save data
    pRelTime(i-1,1:end)=pRel;
    pRel_STFTime(i-1,1:end) = pRel_STF;
    
    
    saveFile=saveFile+1;
    if saveFile>length(ti)
        disp('save file');
        strFile=['solRK_',DF,Sm,'_',int2str(numFile),'.mat'];
        save(strFile,'ti','wi','pRelTime','pRel_STFTime');
        saveFile=0;
        numFile=numFile+1;
        clear wi; clear ti; clear pRelTime; clear pRel_STFTime;
%         wi(1:neqn, 1) = x0';
        ti=zeros(longData);
        pRel_STFTime=zeros(longData,3*nNeurons);
        pRelTime=zeros(longData,3*nNeurons);
        i=1;
    end

  
    
end

% fclose(fileSolution);
    
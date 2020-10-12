%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  MAIN PROGRAM  ALL FUNCTIONS TOGETHER %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc

tstart = tic;

disp('       STD - Network       ');

% NETWORK %

% Network Connectivity
% ProbInh=0.2;
% nsynapses=20;
% [P, ExcInh]=NetworkConnectivity(nNeurons,ProbInh,nsynapses,xDist,sigma_exc,sigma_inh);

% Read the builded network
disp('Reading P');
BN=load('BuildNetwork320.mat');
P=BN.P;
disp('Reading EI');
ExcInh=BN.ExcInh;

nNeurons=length(ExcInh);

% DEPRESSION, FACILITATION and BICUCULLINE FACTORS
Dnumber=1.00;                  % 1 means no depression
Fnumber=1.0;                   % 1 means no facilitation
BIC=0;                         % 0 means no Bicuculline
Sm='100';                      % Write de number of Depression, of facilitation or of bicuculline
iduration=30;                  % Duration of the stimulus
AmplitudCurrent=1;             % Amplitud of the stimulus

DF='Iapp1_tOn3670_dt30_D';     % WRITE:
                                 % D for depression
                                 % F for facilitation
                                 % BIC for Bicuculline
                                 % Iapp1_tOn3670_dt30_D  for applied current

                                 
% RUNNING TIME
tf=10000;        % Time of simulation in miliseconds (ms)
t0=0;   
titer=floor(tf/0.05) + 1;
numFile=floor((titer-1)/1000)-1;

% EXTERNAL APPLIED CURRENT     % We want to apply external current
tON=3670;    % ms
tOFF=tON + iduration;   % ms
choiceNeuron=92;                % Injected neuron   
Iinj=zeros(1,nNeurons);         % It has to be inicializated to inject current
for numInj=0:3
    Iinj(1,choiceNeuron+numInj)=AmplitudCurrent;         % Neurons that we want to inject.
end


% Iinj1=zeros(1,nNeurons);
% Iinj2=zeros(1,nNeurons);




%%% FIXED PARAMETERS %%%
xDist=5000;		% Measured distance in micrometers
sigma_exc=250;	% Exitatory variance in micrometers
sigma_inh=125;	% Inhibitory variance
%neuprint=20;	% Random neurons to be printed with all values
%variables=47;   % Number of principle variables - InOutput file

%nNeurons=1280;
%nNeurons=320;
% nprint=10;      % Number of steps between outputs (20 = 1ms)
neq=20;         % Number of equations (variables)
Vthre=-50.;		
potasio=0;		% For potasium in field Vk=-100+a and Vk=-90+a
nmdachg=1;		% Porcentual change for reducing gNMDA (1 no change)
gabachg=1;
Iahpchg=1;
TOL=0.05;       % Runge-Kuta step
% moreStimul=1;   % Number of Pre-synaptic neurons to stimulate
% nMedir=1;       % Number of Post-synaptic neurons to be measure
% current=200;    % Frequency of the injected current (every 200ms=5Hz)
% corriente=3.;   % Magnitude of the applied current Iapp
% tms=1;          % Time the current is injected  1ms = 20steps of simulation
% presynaptic=1;  % Number of presynaptic to print different from stimulate
% postsynaptic=1; % Number of postsynaptic to print
% n1=0;
n2=0;
% n3=0;
% n4=0;
% nIter=0;
% nSec=1;
% volt=0;
% voltcont=1;
% fs=0;
% seconds=0;
% spiketrain=0;	% If 1 then the injection will be with spike trains
% freqtrain=20;	% Frequency of the spikes during the train (50ms = 20Hz)
% strain=5;		% Number of spikes during the train
fD_AMPA=Dnumber;	% Synaptic Depression Factor (0<fD<1)
fD_NMDA=Dnumber;
fD_GABA=Dnumber;
fF_AMPA=Fnumber;		% Synaptic Facilitation Factor (0<fF<1)
fF_NMDA=Fnumber;
fF_GABA=Fnumber;

tau_relAMPA=400;    % Rate-Time for release probability (depression)
tau_relNMDA=400;
tau_relGABA=400;
tau_stfAMPA=50;		% Rate-Time for release probability (facilitation)
tau_stfNMDA=50;
tau_stfGABA=50;
% fre=0;				% Excitatory Frequency
% fri=0;
% arch_tmp=0;			% If 1 the Voltages.out and pRel.out will be printed
% tIapp=12000;    % Time where the injected current starts
% tIapp=tIapp+fix((normrnd(0,1)*100));  % time variable where the injection is applied
%                                       % round((normrnd(0,1)*100)) per
%                                       % arrodonir

nvar=neq*nNeurons;  %Number of variables to integrate


%%% VECTORS NEEDED %%%
ex=zeros(1,nvar+1);  % Initial condition vector
%vv=zeros(1,nvar+1);
%dydx=zeros(1,nvar+1);
%pre=zeros(1,nvar+1);
%presyn=zeros(1,nNeurons);


leak1=zeros(1,nNeurons);
leak2=zeros(1,nNeurons);
gsdrnd=zeros(1,nNeurons);

%Ptmp=zeros(1,nNeurons*nNeurons);
%eitmp=zeros(1,nNeurons*4);
%SynExc=zeros(1,nNeurons);  %Number of excitatory sinaptic connexions per neuron
%SynInh=zeros(1,nNeurons);
Nmeasure=zeros(1,nNeurons);

%spikevs=zeros(1,nNeurons);
%spikevd=zeros(1,nNeurons);
%spikeintN=zeros(1,nNeurons);

pRel=zeros(1,3*nNeurons);
pRel_STF=zeros(1,3*nNeurons);


%%% TAU VALUES %%%
% Release Probability: Assuming the same for AMPA, NMDA and GABA
% pRel just affects to excitation neurons

p0_AMPA=1;    %Initial probability of release
p0_NMDA=1;
p0_GABA=1;
p0_stfAMPA=1;
p0_stfNMDA=1;
p0_stfGABA=1;

factor_relAMPA=exp(-TOL/tau_relAMPA);	% Depression
factor_relNMDA=exp(-TOL/tau_relNMDA);
factor_relGABA=exp(-TOL/tau_relGABA);

factor_stfAMPA=exp(-TOL/tau_stfAMPA);   % Facilitation
factor_stfNMDA=exp(-TOL/tau_stfNMDA);
factor_stfGABA=exp(-TOL/tau_stfGABA);

pk_AMPA=1;
pk_NMDA=1;
pk_GABA=1;

pk_stfAMPA=1;
pk_stfNMDA=1;
pk_stfGABA=1;




disp('Initializations and SD');
% leakage
for j=1:nNeurons
    pRel(j)=1;                %For depression 
    pRel(j+nNeurons)=1;
    pRel(j+2*nNeurons)=1;
    
    pRel_STF(j)=1;             %For facilitation
	pRel_STF(j+nNeurons)=1;
	pRel_STF(j+2*nNeurons)=1;
    
    SD=SDNumber();
    if ExcInh(j)==0 % The neuron is excitatory
		leak1(j)=0.0667+(0.0067*SD); % gL of neuron j
    else
		leak1(j)=0.1025+(0.0025*SD);
    end
    
    SD=SDNumber();
    if ExcInh(j)==0
		leak2(j)=-60.95+(0.3*SD); % vL of neuron j
    else
		leak2(j)=-63.8+(0.15*SD);
    end
    
    SD=SDNumber();
    if ExcInh(j)==0
        gsdrnd(j)=0.1*SD; % SD of gsd (i.e. gsd = 1.75+gsdrnd)
    else
        gsdrnd(j)=0;
    end
end

% ex=x0 (the Initial condition)
for j=1:nNeurons
    if ExcInh(j)==0 %Pyramidal Exc.
        ex((j-1)*neq+2)=-60+(normrnd(0,1)*5);
		ex((j-1)*neq+3)=-60+(normrnd(0,1)*5);
    else  % Interneurons Inh.
        ex((j-1)*neq+14)=-60-(normrnd(0,1)*5);
    end
end
    
% nw=1;               % The first output (nw=0) is printed apart
% iny=current;        % Counter of injection on external frequency
% train=freqtrain;    % Conunter for the spike train, the frequency is indicated
	
cont=1;
cont2=1;
renglon=1;
num_sp=1;


N=(tf-t0)/TOL;

disp('Rk45');
[pRel,pRel_STF]=rk45Network ('NetworkField', t0, ex, tf, N , Iinj, ExcInh, nvar, neq, P, pRel, pRel_STF, Nmeasure, potasio, n2, leak1, leak2, gsdrnd, 1,Vthre, p0_AMPA,p0_NMDA,p0_GABA,factor_relAMPA,factor_relNMDA,factor_relGABA,p0_stfAMPA,p0_stfNMDA,p0_stfGABA,factor_stfAMPA,factor_stfNMDA,factor_stfGABA,fD_AMPA,fD_NMDA,fD_GABA,fF_AMPA, fF_NMDA,fF_GABA,nNeurons,Sm,tON,tOFF,DF,BIC);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   TABLE VOLTAGE   %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except Sm titer numFile ExcInh DF
file=['solRK_',DF,Sm,'_'];
variable = ['voltage_',DF,Sm];

conn=load('BuildNetwork320.mat');
EI=conn.ExcInh;

nNeuron=320;

neq=20;
var=neq*nNeuron;
Time=zeros(titer,1);
Voltage=zeros(titer,nNeuron);


k=1;
for i=1:numFile
%    disp('Reading File %d',i);
   str=[file,int2str(i),'.mat'];
   
   CurrentFile=load(str);
   T=CurrentFile.ti;
   Sol=CurrentFile.wi;
   Sol=transpose(Sol);
   
   % Voltage
   for j=1:length(T)
       Time(k,1)=T(j,1);
       for s=1:nNeuron
           if EI(s)==0
               Voltage(k,s)=Sol(j,2+(s-1)*neq);
           else
               Voltage(k,s)=Sol(j,14+(s-1)*neq);
           end
 
       end
       k=k+1;
   end
   
   clear T; clear Sol; 

end
save(variable,'Time','Voltage');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   TREAT DATA   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except Sm titer numFile ExcInh DF
strCurrent=['voltage_',DF,Sm,'.mat'];
CurrentFile=load(strCurrent);
strFD=['f',DF,'=',Sm];
strOscillation=['oscillation320-',DF,Sm];


T=CurrentFile.Time;
rkSol=CurrentFile.Voltage;

neq=20;
nNeuron=320;
Vthre=-30;
k=1;

for j=1:nNeuron
    for ts=2:length(T)-1
        if ((rkSol(ts,j)>Vthre) && rkSol(ts-1,j)<rkSol(ts,j) && rkSol(ts,j)>rkSol(ts+1,j))
            TIMEspike(1,k)=T(ts,1);
            TIMEspike(2,k)=j;
            k=k+1;
        end
    end
    
end

h=figure(1);

plot(TIMEspike(1,:),TIMEspike(2,:),'.');
set(gca,'FontSize',14);
title(strFD,'FontSize',16);
ylabel('Neurons','FontSize',16);
xlabel('time (ms)','FontSize',16);
ylim([1 320])
xlim([0 10000])
saveas(h,strOscillation);
saveas(h,strOscillation,'jpg');




save(strOscillation,'TIMEspike');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   TABLE CONDUCTANCES   %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except Sm titer numFile DF
numD=1;

file=['solRK_',DF,Sm,'_'];

conn=load('BuildNetwork320.mat');
EI=conn.ExcInh;

nNeuron=320;
neq=20;
var=neq*nNeuron;
time=zeros(titer,1);
gAMPA=zeros(titer,nNeuron);
gGABA=zeros(titer,nNeuron);
gNMDA=zeros(titer,nNeuron);


k=1;
for i=1:numFile
    %    disp('Reading File %d',i);
    str=[file,int2str(i),'.mat'];
    
    %         cd(way);
    CurrentFile=load(str);
    %         cd('/Users/CatiVich/Documents/MATLAB/Matlab programes/XarxaSTD')
    T=CurrentFile.ti;
    Sol=CurrentFile.wi;
    Sol=transpose(Sol);
    Prel=CurrentFile.pRelTime;
    Prel_STF=CurrentFile.pRel_STFTime;
    
    % Voltage
    
    for j=1:length(T)
        time(k,1)=T(j,1);
        
        for s=1:nNeuron
            if EI(1,s)==0 % Excitadora
                gAMPA(k,s)=(5.4/10000)*Sol(j,10+(s-1)*neq)*Prel(j,s)*Prel_STF(j,s);               
                gNMDA(k,s)=(0.9/10000)*Sol(j,11+(s-1)*neq)*Prel(j,s+nNeuron)*Prel_STF(j,s+nNeuron);       
                gGABA(k,s)=(4.15/10000)*Sol(j,13+(s-1)*neq)*Prel(j,s+2*nNeuron)*Prel_STF(j,s+2*nNeuron);
            else   % Inhibidora
                gAMPA(k,s)=(2.25/10000)*Sol(j,17+(s-1)*neq)*Prel(j,s)*Prel_STF(j,s);               
                gNMDA(k,s)=(0.5/10000)*Sol(j,18+(s-1)*neq)*Prel(j,s+nNeuron)*Prel_STF(j,s+nNeuron);       
                gGABA(k,s)=(0.165/10000)*Sol(j,20+(s-1)*neq)*Prel(j,s+2*nNeuron)*Prel_STF(j,s+2*nNeuron);
            end                                   
        end
        k=k+1;
    end
    clear T; clear Sol;
    
end

varAMPA = ['gAMPA_',DF,Sm,'_',int2str(nNeuron),'N'];
varGABA = ['gGABA_',DF,Sm,'_',int2str(nNeuron),'N'];
varNMDA = ['gNMDA_',DF,Sm,'_',int2str(nNeuron),'N'];


save(varAMPA,'time','gAMPA');
save(varNMDA,'time','gNMDA');
save(varGABA,'time','gGABA');


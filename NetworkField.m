function dx=NetworkField(t,x, I, EI, nvar, neq, Prob, pRel, pRel_STF, Medir, potasio, npreinj, leak1, leak2, gsdrnd, sumes,tON,tOFF,BIC)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compte et al's pyramidal cell model %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% t= time
% x= vector field
% param = (I, pre, EI, nvar, neq, Prob, pRel, pRel_STF, nw, temp,
%          Medir, potasio, npreinj, nmda_chg, gaba_chg, Iahp_chg, 
%          arch_currents, leak1, leak2, gsdrnd, ninh_medir, sumes, 
%          aleatSingleExc, aleatSingleInh)
%

dx=zeros(1,nvar+1);

% Initilizations
sIl=0;
sIna=0;
sIk=0;
sIa=0;
sIks=0;
sIkna=0;
sIca=0;
sIkca=0;
sINap=0;
sIar=0;
sEGABA=0;
sENMDA=0;
sEAMPA=0;
sCoupling1=0;
sCoupling2=0;
nExc=0;
	
sIIl=0;
sIIna=0;
sIIk=0;
sIGABA=0;
sINMDA=0;
sIAMPA=0;
nInh=0;


% Parameters
Cm=1.0;
% gl, gna, gk
ga=1; gks=0.576; gkna=1.33; gkca=0.57; gca=0.43; gnap=0.0686; gar=0.0257;
vca=120; % vl, vna, vk, a

    % Somas and dendrite areas
As=0.015;
Ad=0.035;
AintN=0.02;
%ninh2=ninh_medir;
neuron2=npreinj;
a=potasio;

    % Synaptic dynamics parameters
aAMPA=3.48;
aNMDA=0.5;
aGABA=1.0;
ax=3.48;
taux=2;
tauAMPA=2.0;
tauNMDA=100.0;
tauGABA= 10.0;

jNaux=nvar/neq;

currpre=zeros(1,19);
currmed=zeros(1,19);

sAMPA=zeros(1,jNaux);
sNMDA=zeros(1,jNaux);
sGABA=zeros(1,jNaux);

synAMPA=zeros(1,jNaux);
synNMDA=zeros(1,jNaux);
synGABA=zeros(1,jNaux);

for jN=0:jNaux-1
    index=jN*neq;
   
    sAMPA(jN+1)=x(index+10);
    sNMDA(jN+1)=x(index+11);
    sGABA(jN+1)=x(index+13);
    
    synAMPA(jN+1)=x(index+17);
    synNMDA(jN+1)=x(index+18);
    synGABA(jN+1)=x(index+20);
end


% Calculating the equations for every neuron.
for jN=0:jNaux-1
    index=jN*neq;
    
    % ionic channels kinetics variables
    vs=x(index+2); % soma's voltage
    vd=x(index+3); % dendrite's voltage
    h=x(index+4);
    n=x(index+5);
    ha=x(index+6);
    mks=x(index+7);
    Na=x(index+8);
    Ca=x(index+9);
    xNMDAs=x(index+12);
     
    % interneurons variables
    v=x(index+14); % voltage of the current neuron
    h_intN=x(index+15);
    n_intN=x(index+16);
    xNMDA=x(index+19);
    
    % other parameters
    
    
    if t>=tON && t<tOFF
        Iapp=I(jN+1);
    else
        Iapp=0;
    end
    
    ant_excinh=EI(jN+1);
    excinh=EI(jN+1);
    leakage1=leak1(jN+1);
    leakage2=leak2(jN+1);
    
    gsd=(1.75+gsdrnd(jN+1))*0.1; %0.1 perquè quadrin les unitats
    
    neuron=Medir(jN+1);
    
    nexc=0;
    ninh=0;
    
    
    fact_AMPA=0;
    fact_NMDA=0.;
	fact_GABA=0.;
	intfact_AMPA=0.;
	intfact_NMDA=0.;
	intfact_GABA=0.;
    
    if(excinh==0)   % We have not synaptic device. We just are look for the 
                    % currents into the soma and the dendrite
        
        %%% currents in SOMA (equations for vs) %%%
        
        % leak
        gl=leakage1;
        vl=leakage2;
        Il=gl*(vs-vl);
        
        % sodium
        if neuron==1
            Ina=0.0;
        else
            vNa=55;
            gna=50;
            am=0.1*(vs+33)/(1-exp(-(vs+33)/10));
            bm=4*exp(-(vs+53.7)/12);
            minf=am/(am+bm);
            m=minf;
            m3=m*m*m;
            Ina=gna*m3*h*(vs-vNa);
            ah=0.07*exp(-(vs+50)/10);
            bh=1/(1+exp(-(vs+20)/10));
        end
        phi=4.0;
        
        % delayed rectifier
        vk=-100+a;
        gk=10.5;
        n2=n*n;
        Ik=gk*(n2*n2)*(vs-vk);
        an=0.01*(vs+34)/(1-exp(-(vs+34)/10));
        bn=0.125*exp(-(vs+44)/25);
        
        % Fast A-type K channel
        haInf=1/(1+exp((vs+80)/6));
        maInf=1/(1+exp(-(vs+50)/20));
        Ia=ga*maInf*maInf*maInf*ha*(vs-vk);
        tauHa=15;
        phiHa=1;
        
        % non-inactivating K channel
        Iks=gks*mks*(vs-vk);
        mksinf=1/(1+exp(-(vs+34)/6.5));
        tks=8/(exp(-(vs+55)/30)+exp((vs+55)/30));
        phiks=1;
        
        % Na dependent K channel
        if neuron==1
            Ikna=0.0;
        else
            wNa=0.37/(1+(38.7/Na)^3.5); 
            Ikna=gkna*wNa*(vs-vk);
            ana=0.01*10; % *10 perquè les unitats quadrin a mM 
            NaEq=9.5;
            NaEq3=NaEq*NaEq*NaEq;
            Rpump=0.018;
        end
        
        
        %%% currents in DENDRITE (equations for vD) %%%
        
        % Calcium channel
        mCainf=1/(1+exp(-(vd+20)/9));
       	Ica=gca*mCainf*mCainf*(vd-vca);
        
        % The Ca dependent K channel
        Kd=30;
        Ikca=gkca*Ca/(Ca+Kd)*(vd-vk);
        aCa=0.005*10;  % *10 perquè les unitats quadrin a mM 
        tauCa=150;
        
        % Persistens sodium channel
        if neuron==1
            INap=0.0;
        else
            mNapinf=1/(1+exp(-(vd+55.7)/7.7));
            INap=gnap*mNapinf*mNapinf*mNapinf*(vd-vNa);
            
        end
        
        % Inward rectifying K channel
        hArinf=1/(1+exp((vd+75)/4.));
        Iar=gar*hArinf*(vd-vk);

        % Synapsis
        VsynAMPA=0.0;
        VsynNMDA=0.0;
        VsynGABA=-70;
        
            % Adding the values of all posible synapsis
        for j=1:jNaux
           %index2=(j-1)*neq;
            pRels_AMPA=pRel(j);
            pRels_NMDA=pRel(j+jNaux);
            pRels_GABA=pRel(j+2*jNaux);
            pRels_stfAMPA=pRel_STF(j);
			pRels_stfNMDA=pRel_STF(j+jNaux);
			pRels_stfGABA=pRel_STF(j+2*jNaux);
            
            gAMPA=5.4/10000;
            gNMDA=0.9/10000;
            gGABA=(4.15/10000)*(1-BIC);
            
            fact_AMPA=gAMPA*(sAMPA(j)*Prob(jN+1,j))*pRels_AMPA*pRels_stfAMPA+fact_AMPA;
			fact_NMDA=gNMDA*(sNMDA(j)*Prob(jN+1,j))*pRels_NMDA*pRels_stfNMDA+fact_NMDA;
            fact_GABA=gGABA*(sGABA(j)*Prob(jN+1,j))*pRels_GABA*pRels_stfGABA+fact_GABA;
            
        end
        
            % Adding the current neuron
        fVpre_s=(1/(1+exp(-(vs-20)/2)));
        
        if ninh<=0
            ninh=1;
        end
        if nexc<=0
            nexc=1;
        end
        
        
        %%% Summed currents per blocks %%%
        Tot_exc=Il+Ina+Ik+Ia+Iks+Ikna;      % Total Excitation current
		Tot_inh=Ica+Ikca+INap+Iar;          % Total Inhibition current
		Isyns=(fact_GABA)*(vs-VsynGABA);    % Total Synaptic current at soma
		Isynd1=(fact_NMDA)*(vd-VsynNMDA);
        Isynd2=(fact_AMPA)*(vd-VsynAMPA);
        Isynd=Isynd1+Isynd2;                % Total Synaptic current at dendrite
        
        
        
        %%% equations for pyramical neurons %%%
        dx(index+2)=(-(Il+Ina+Ik+Ia+Iks+Ikna)-Isyns/As-gsd*(vs-vd)/As+Iapp/As)/Cm;
        dx(index+3)=(-(Ica+Ikca+INap+Iar)-Isynd/Ad-gsd*(vd-vs)/Ad)/Cm;
        dx(index+4)=phi*(ah*(1-h)-bh*h);
        dx(index+5)=phi*(an*(1-n)-bn*n);
        dx(index+6)=phiHa*(haInf-ha)/tauHa;
        dx(index+7)=phiks*(mksinf-mks)/tks;
		dx(index+8)=-ana*(As*Ina+Ad*INap)-Rpump*((Na*Na*Na)/((Na*Na*Na)+3375)-NaEq3/(NaEq3+3375));
        dx(index+9)=-aCa*(Ad*Ica)-(Ca/tauCa);

        dx(index+10)=aAMPA*fVpre_s-sAMPA(jN+1)/tauAMPA;
        dx(index+11)=aNMDA*xNMDAs*(1-sNMDA(jN+1))-sNMDA(jN+1)/tauNMDA;
        dx(index+12)=ax*fVpre_s-xNMDAs/taux;
        dx(index+13)=aGABA*fVpre_s-sGABA(jN+1)/tauGABA; %0
       
        for j=14:20
            dx(index+j)=0;   % Interneurons derivations are equal to 0
        end
        
        % Resting excitatory currents in the soma
        if sumes==1  
            sIl=sIl-Il;
            sIna=sIna-Ina;
            sIk=sIk-Ik;
            sIa=sIa-Ia;
            sIks=sIks-Iks;
            sIkna=sIkna-Ikna;
            sIca=sIca-Ica;
            sIkca=sIkca-Ikca;
            sINap=sINap-INap;
            sIar=sIar-Iar;
            sEGABA=sEGABA-(fact_GABA)*(vd-VsynGABA)/As;
            sENMDA=sENMDA-(fact_NMDA)*(vd-VsynNMDA)/Ad;
            sEAMPA=sEAMPA-(fact_AMPA)*(vd-VsynAMPA)/Ad;
            sCoupling1=sCoupling1-gsd*(vs-vd)/As;
            sCoupling2=sCoupling2-gsd*(vd-vs)/Ad;
            nExc=nExc+1; % We add 1 excitation into the counter.
        end
        
        
    else % We have synaptic device. Then, we treat the interneurons

        % Leakage
        gl=leakage1;
        vl=leakage2;
        Il=gl*(v-vl);
        
        % Sodium
        gna=35.0;
        vNa=55.0;
        am=0.5*(v+35)/(1-exp(-(v+35)/10));
        bm=20*exp(-(v+60)/18);
        minf=am/(am+bm);
        m=minf;
        m3=m*m*m;
        Ina=gna*m3*h_intN*(v-vNa);
        ah=0.35*exp(-(v+58)/20);
        bh=5/(1+exp(-(v+28)/10));
        phi_intN=1;
        
        % Delayed rectifier potassium
        gk=9.0;
        vk=-90.0+a;
        n2=n_intN*n_intN;
        Ik=gk*(n2*n2)*(v-vk);
        an=0.05*(v+34)/(1-exp(-(v+34)/10));
        bn=0.625*exp(-(v+44)/80);
        
        % Synapsis
        VsynAMPA=0.0;
        VsynNMDA=0.0;
		VsynGABA=-70.0;
        
            % Adding the values of all possible synapses
        for j=1:jNaux
            %index2=(j-1)*neq;
			pRels_AMPA=pRel(j);
			pRels_NMDA=pRel(j+jNaux);
			pRels_GABA=pRel(j+2*jNaux);

			pRels_stfAMPA=pRel_STF(j);
			pRels_stfNMDA=pRel_STF(j+jNaux);
			pRels_stfGABA=pRel_STF(j+jNaux);
				
			gAMPA=2.25/10000;     %original value 2.25
			gNMDA=0.5/10000;      %original (0.5/10000)
			gGABA=0.165/10000;    %original 0.165/10000 for modified (Wang proportion) 1.73
            
            intfact_AMPA=gAMPA*(sAMPA(j)*Prob(jN+1,j))*pRels_AMPA*pRels_stfAMPA+intfact_AMPA;
			intfact_NMDA=gNMDA*(sNMDA(j)*Prob(jN+1,j))*pRels_NMDA*pRels_stfNMDA+intfact_NMDA;
            intfact_GABA=gGABA*(sGABA(j)*Prob(jN+1,j))*pRels_GABA*pRels_stfGABA+intfact_GABA;
        end
        fVpre=(1/(1+exp(-(v-20)/2)));
		
		Tot_intN=Il+Ina+Ik; %total interneuron current
			
		IsynintN1=(intfact_GABA)*(v-VsynGABA);
		IsynintN2=(intfact_NMDA)*(v-VsynNMDA);
		IsynintN3=(intfact_AMPA)*(v-VsynAMPA);
        IsynintN=IsynintN3+IsynintN2+IsynintN1;
        
        %%% Equations for interneurons %%%
        dx(index+14)=(-(Il+Ina+Ik)-IsynintN/AintN+Iapp/AintN)/Cm;
		dx(index+15)=phi_intN*(ah*(1-h_intN)-bh*h_intN);
        dx(index+16)=phi_intN*(an*(1-n_intN)-bn*n_intN);

        dx(index+17)=aAMPA*fVpre-synAMPA(jN+1)/tauAMPA;   %0
        dx(index+18)=aNMDA*xNMDA*(1-synNMDA(jN+1))-synNMDA(jN+1)/tauNMDA; %0
        dx(index+19)=ax*fVpre-xNMDA/taux; %0
	    dx(index+20)=aGABA*fVpre-synGABA(jN+1)/tauGABA; %0
        
        for j=2:13
            dx(index+j)=0;
        end
        
        if sumes==1         % Resting inhibitory currents in the soma
            sIIl=sIIl-Il;
            sIIna=sIIna-Ina;
            sIIk=sIIk-Ik;
            sIGABA=sIGABA-(intfact_GABA)*(v-VsynGABA);
            sINMDA=sINMDA-(intfact_NMDA)*(v-VsynNMDA);
            sIAMPA=sIAMPA-(intfact_AMPA)*(v-VsynAMPA);
            nInh=nInh+1; % We add 1 inhibition into the counter.
        end

    end
       
    if EI(jN+1)==0  % Sum of Currents of all excitatory neurons
		currpre(1)=Ikna+currpre(1);
		currpre(2)=Na+currpre(2);
		currpre(3)=wNa+currpre(3);
		currpre(4)=dx(index+8)+currpre(4);
		currpre(5)=-ana*(As*Ina+Ad*INap)+currpre(5);
		currpre(6)=-Rpump*((Na*Na*Na)/((Na*Na*Na)+3375)-NaEq3/(NaEq3+3375))+currpre(6);
		currpre(7)=Ina+currpre(7);
		currpre(8)=INap+currpre(8);
        
        if (jN+1)==neuron2 % Current of a presynaptic (injected) neuron
            currmed(1)=Ikna;
			currmed(2)=Na;
			currmed(3)=wNa;
			currmed(4)=dx(index+8);
			currmed(5)=-ana*(As*Ina+Ad*INap);
			currmed(6)=-Rpump*((Na*Na*Na)/((Na*Na*Na)+3375)-NaEq3/(NaEq3+3375));
			currmed(7)=Ina;
			currmed(8)=INap;
        end
    end
       
end

dx(1)=0;





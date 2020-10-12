function [pRel,pRel_STF] = Prelease(j,neq,nNeurons,ExcInh,pRel,pRel_STF, ex, pre, Vthre, p0_AMPA,p0_NMDA,p0_GABA,factor_relAMPA,factor_relNMDA,factor_relGABA,p0_stfAMPA,p0_stfNMDA,p0_stfGABA,factor_stfAMPA,factor_stfNMDA,factor_stfGABA,fD_AMPA,fD_NMDA,fD_GABA,fF_AMPA, fF_NMDA,fF_GABA)

pRel(j)= p0_AMPA*(1-factor_relAMPA)+pRel(j)*factor_relAMPA;
pRel(j+nNeurons)=p0_NMDA*(1-factor_relNMDA)+pRel(j+nNeurons)*factor_relNMDA;
pRel(j+2*nNeurons)=p0_GABA*(1-factor_relGABA)+pRel(j+2*nNeurons)*factor_relGABA;

pRel_STF(j)=p0_stfAMPA*(1-factor_stfAMPA)+pRel_STF(j)*factor_stfAMPA;
pRel_STF(j+nNeurons)=p0_stfNMDA*(1-factor_stfNMDA)+pRel_STF(j+nNeurons)*factor_stfNMDA;
pRel_STF(j+2*nNeurons)=p0_stfGABA*(1-factor_stfGABA)+pRel_STF(j+2*nNeurons)*factor_stfGABA;

l=(j-1)*neq;
if ExcInh(j)==0
    if (pre(l+2)-Vthre <0) && (ex(l+2)-Vthre>0)
        pRel(j)=pRel(j)*fD_AMPA;
        pRel(j+nNeurons)=pRel(j+nNeurons)*fD_NMDA;
		pRel_STF(j)=pRel_STF(j)+(1-pRel_STF(j))*fF_AMPA;
		pRel_STF(j+nNeurons)=pRel_STF(j+nNeurons)+(1-pRel_STF(j+nNeurons))*fF_NMDA;
    else
        pRel(j)=pRel(j);
        pRel_STF(j)=pRel_STF(j);
    end
else
    if((pre(l+14)-Vthre)<0 && (ex(l+14)-Vthre)>0)
        pRel(j+2*nNeurons)=pRel(j+2*nNeurons)*fD_GABA;
		pRel_STF(j+2*nNeurons)=pRel_STF(j+2*nNeurons)+(1-pRel_STF(j+2*nNeurons))*fF_GABA;
    else
        pRel(j)=pRel(j);
        pRel_STF(j)=pRel_STF(j);
    end
end

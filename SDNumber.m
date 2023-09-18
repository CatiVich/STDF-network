function SD = SDNumber()
SD=normrnd(0,1);
if (SD>1)
    SD=SD-floor(SD);
end
if (SD<-1)
    SD=SD-ceil(SD);
end
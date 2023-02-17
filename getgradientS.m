function GS = getgradientS(X,Y,mk,Delta,type)

covMX = rankCovIID(X);
covMY = rankCovIID(Y);
p=size(covMY,1);
M_k=reshape(mk,p,p);
if type==1
    TX=2/pi*asin(covMX);
    GS=(Delta*covMY*M_k')'.*(pi/2*cos(pi/2*TX))+(M_k'*covMY*Delta)'.*(pi/2*cos(pi/2*TX));
    GS=GS/2-M_k.*(pi/2*cos(pi/2*TX));
else
    TY=2/pi*asin(covMY);
    GS=(Delta*covMX*M_k')'.*(pi/2*cos(pi/2*TY))+(M_k'*covMX*Delta)'.*(pi/2*cos(pi/2*TY));
    GS=GS/2+M_k.*(pi/2*cos(pi/2*TY));
end

end

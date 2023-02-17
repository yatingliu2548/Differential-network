function eStd3 = varmatrix(covMX,covMY,sigmagx,sigmagy,Delta,iter,n1,n2)
% this function is use to get the estimator of variance matrix of gradient lasso

p=size(covMX,1);
TX=2/pi*asin(covMX);
TY=2/pi*asin(covMY);
if iter>1
    GS=(Delta*covMY).*(pi/2*cos(pi/2*TX))+(covMY*Delta).*(pi/2*cos(pi/2*TX));
    GS=GS/2-(pi/2*cos(pi/2*TX));
    GSY=(Delta*covMX).*(pi/2*cos(pi/2*TY))+(covMX*Delta).*(pi/2*cos(pi/2*TY));
    GSY=GSY/2-(pi/2*cos(pi/2*TY));
else
    GS=(pi/2*cos(pi/2*TX));
    GSY=(pi/2*cos(pi/2*TY));
end
diagsigmagx=diag(sigmagx);
matrixsigmagx=reshape(diagsigmagx,p,p);
diagsigmagy=diag(sigmagy);
matrixsigmagy=reshape(diagsigmagy,p,p);
eStd3=4/n1*GS.*GS.*matrixsigmagx+4/n2*GSY'.*GSY.*matrixsigmagy;
end
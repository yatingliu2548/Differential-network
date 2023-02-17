
function var = jacknifeEstimatorVar(X,Y,mk,Delta,type)

[nx, p] = size(X);
[ny, p] = size(Y);
if type==1
    S0=gethatU(X,Y,mk,Delta);
elseif type==2 
    S0=getstatisticsS(X,Y,mk,Delta);
elseif type==3 
    S0=gethatdiff(X,Y,0,0);
else
    S0=gethatdiff(X,Y,1,2);
end

Sadim=zeros(nx,1);
for a=1:nx
    X2=X(setdiff(1:nx,a),:);
    if type==1
        Sa=gethatU(X2,Y,mk,Delta);
    elseif type==2
        Sa=getstatisticsS(X2,Y,mk,Delta);
    elseif type==3
        Sa=gethatdiff(X2,Y,0,0);
    else
        Sa=gethatdiff(X2,Y,1,2);
    end

    Sadim(a)=nx*S0-(nx-1)*Sa;
end
Sx=1/nx*sum(Sadim);
Sbdim=zeros(ny,1);
for b=1:ny
    Y2=X(setdiff(1:ny,b),:);
    if type==1
        Sb=gethatU(X,Y2,mk,Delta);
    elseif type==2
        Sb=getstatisticsS(X,Y2,mk,Delta);
    elseif type==3
        Sb=gethatdiff(X,Y2,0,0);    
    else 
        Sb=gethatdiff(X,Y2,1,2);
    end
    
    Sbdim(b)=ny*S0-(ny-1)*Sb;
end
Sy=1/ny*sum(Sbdim);
%calculate var
Sxvec=zeros(nx,1);
Sxvec=Sx;
Syvec=zeros(ny,1);
Syvec=Sy;
var=sum((Sadim-Sxvec).*(Sadim-Sxvec))/nx/(nx-1)+sum((Sbdim-Syvec).*(Sbdim-Syvec))/ny/(ny-1);


end

function [hatDelta,hatM,hatGamma] = oracleEstimator(S,covMX,covMY)

p=size(covMX,1);
[row,col]=find(S);
D=[row,col];
ss=size(D,1);
D2 = repmat(D,ss,1);
D1 = kron(D,ones(ss,1));
D3=[D1,D2];

%
Gamma=zeros(p^2,p^2);
hatM=zeros(p^2,p^2);
for rep=1:size(D3,1)
    j=D3(rep,1);
    k=D3(rep,2);
    l=D3(rep,3);
    m=D3(rep,4);
    r=p*(j-1)+k;
    s=p*(l-1)+m;
    Gamma(r,s)=1/2*(covMX(j,l)*covMY(k,m)+covMY(j,l)*covMX(k,m));
end
[G_row,G_col]=find(Gamma);
GammaSS=Gamma(unique(G_row),unique(G_col));
hatMSS=inv(GammaSS);
hatGamma=Gamma;
hatM(unique(G_row),unique(G_col))=hatMSS;
%diff
k=sub2ind(size(S),row,col);
covMXvec=reshape(covMX,[],1);
covMYvec=reshape(covMY,[],1);
diff=covMXvec(k)-covMYvec(k);
hatDeltaSS=GammaSS\diff;
%
hatDelta=zeros(p,p);
hatDelta(k)=hatDeltaSS;

end
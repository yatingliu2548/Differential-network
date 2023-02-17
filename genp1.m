function result = genp1(p,n,sigma)

A= ones(p, p);
l=find(triu(A,1));
m=randsample(l,n);
theta=zeros(p,p);
theta(m)=normrnd(0,1,n,1)*sigma;
result=theta+theta.';

end


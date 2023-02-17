%verify two ways variance:
% case i
M_k=reshape(hatmk,p,p);
A=1/2*((hatDelta*hatcovMY*M_k')').*(pi/2*cos(pi/2*hatTX));% sigmay \kron sigmax
vecA=A(:);
var1=vecA'*hatSigmagX*vecA;
% case ii
B=(hatmk*(reshape(hatDelta,[],1)')).*(kron(hatcovMY,pi/2*cos(pi/2*hatTX)));
vecB=reshape(B,[],1);
one=ones(size(B,1),1);
onem=one*one.';
var2=vecB'*(kron(onem,hatSigmagX))*vecB;
[var1,var2]

%estimator for case i
A2=hatDelta*hatcovMY*M_k';
A2=A2';
coss=pi/2*cos(pi/2*hatTX);
hatTX(:)'*diag(coss(:))*A2(:)

% original estimator
trace((pi/2*cos(pi/2*hatTX).*(hatTX))*hatDelta*hatcovMY*M_k')
% estimator for case ii
one2=ones(size(hatcovMY,1),1);
onemm=one2*one2';
Q=kron(onemm,hatTX).*kron(hatcovMY,pi/2*cos(pi/2*hatTX));
hatmk'*Q*hatDelta(:)
%
H=(hatmk*hatDelta(:)').*kron(hatcovMY,pi/2*cos(pi/2*hatTX));
B=kron(onemm,hatTX);
H(:)'*B(:)

% running sample covariance for estimator
numRep=200;
n=100;
p1=10;

value1=ones(n,1);
value2=ones(n,1);
value3=ones(n,p*p);
for i=2:n
    signx=sign(repmat(X(i,:),1)-X(1,:))';
    onex=kron(one2,signx);
    %case ii value
    value2(i)=H(:)'*kron(onex,onex);
    %case i value
    value1(i)=A2(:)'*diag(coss(:))*kron(signx,signx);
    value3(i,:)=kron(signx,signx);
end
    
value3=value3';

H(:)'*kron(onem,cov(value3))*H(:)
A3=A2(:)'*diag(coss(:));
A3(:)'*cov(value3)*A3(:)
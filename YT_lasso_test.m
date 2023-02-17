rng(0);
tol = 1e-3;
p=20;
%form sample from huiliyuan and tony cai's method
Omega1 = genp(p,0.2,0.5); %p=50
Omega2 = Omega1+genp1(p,100,0.5);
tdelta = Omega2-Omega1;
SigmaX=Omega1\eye(size(Omega1));%inv(Omega1);
SigmaY=inv(Omega2);
n=200;
mu=zeros(p,1);
X=mvnrnd(mu,SigmaX,n);
Y=mvnrnd(mu,SigmaY,n);
nlambda=10;

lambdaMinRatio=0.04;
trueTX=2/pi*asin(SigmaX);
trueTY=2/pi*asin(SigmaY);
%estimate sigma
hatcovMX = rankCovIID(X);
hatcovMY = rankCovIID(Y);
hatTX=asin(hatcovMX) * 2 / pi;
hatTY=asin(hatcovMY) * 2 / pi;
e=hatcovMY-hatcovMX;
lambdaMax = 2 * max(abs(e(:)));
lambdaMin = lambdaMinRatio * lambdaMax;
lambda = exp(linspace(log(lambdaMax), log(lambdaMin), nlambda));

rho=1;
shrink= 1.5; 
lambda =lambda - shrink*0.001;

Ematrix=ones(p,p);

iternum=1000;
tol_D=10^(-8);

%score is BIC or AIC

TP1seq=cell(length(lambda),1);
TP2seq=cell(length(lambda),1);
TP3seq=cell(length(lambda),1);
score1seq=ones(length(lambda),1);
score2seq=ones(length(lambda),1);
score3seq=ones(length(lambda),1);
distD1seq=ones(length(lambda),1);
distD2seq=ones(length(lambda),1);
distD3seq=ones(length(lambda),1);
distdelta1seq=ones(length(lambda),1);
distdelta2seq=ones(length(lambda),1);
distdelta3seq=ones(length(lambda),1);
hatdelta1seq=cell(length(lambda),1);
hatdelta2seq=cell(length(lambda),1);
hatdelta3seq=cell(length(lambda),1);
% error type in BIC, error = [max_abs_err, l1_err, matrix_l1_err, spectral_err, frobenius_err, nuclear_err]
for la=1:length(lambda)
    D_isE=1;
    %using sample covariance
    [Dseq1,iter1,hatdelta1seq{la},TP1seq{la},true_MatrixD1,loss1,score1seq(la),distD1seq(la),distdelta1seq(la)]=lasso_kendall(X,Y,iternum,lambda(la),rho,tol_D,D_isE,tdelta,1,"BIC",2,0);
    %using kendall's tau
    [Dseq3,iter3,hatdelta3seq{la},TP3seq{la},true_MatrixD3,loss3,score3seq(la),distD3seq(la),distdelta3seq(la)]=lasso_kendall(X,Y,iternum,lambda(la),rho,tol_D,D_isE,tdelta,1,"BIC",1,0);
    
    D_isE=2;
    [Dseq2,iter2,hatdelta2seq{la},TP2seq{la},true_MatrixD2,loss2,score2seq(la),distD2seq(la),distdelta2seq(la)]=lasso_kendall(X,Y,iternum,lambda(la),rho,tol_D,D_isE,tdelta,1,"BIC",1,0);
    
end
[min_value1, min_index1] = min(score1seq);
[min_value2, min_index2] = min(score2seq);
lambdachose=lambda(min_index2);
hatM_lasso=zeros(p*p,p*p);
for i=1:(p*p)
    [hatM_lasso_matrix,iterseq]=L1_dts(hatcovMX,hatcovMY,rho,lambdachose,Ematrix,i);
    hatM_lasso(:,i)=reshape(hatM_lasso_matrix,[],1);
end
TP1seq{min_index1}
TP2seq{min_index2}
%TD,TP,TN
%want TP,TD close to 1
hatdelta1=hatdelta1seq{min_index1};
hatdelta2=hatdelta2seq{min_index2};
% estimated Hessian M
S1=(abs(hatdelta1)>1e-8);
S2=(abs(hatdelta2)>1e-8);
%[row1,col1]=find(S1);
%k1=sub2ind(size(S1),row1,col1);
[~,hatM1,~] = oracleEstimator(S1,hatcovMX,hatcovMY);
[~,hatM2,~] = oracleEstimator(S2,hatcovMX,hatcovMY);
trueS=(abs(tdelta)>1e-8);
[~,trueM,~] = oracleEstimator(trueS,SigmaX,SigmaY);
[~,trueM2,~] = oracleEstimator(trueS,hatcovMX,hatcovMY);

% tic;

% Insert your code here

% elapsed_time = toc;
% fprintf('The elapsed time is %f seconds\n', elapsed_time);

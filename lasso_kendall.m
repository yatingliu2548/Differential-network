function [Dseq,iter2,hatdelta2,TD,true_MatrixD,loss,score,distD,distdelta]=lasso_kendall(X,Y,iternum,lambda,rho,tol_D,D_isE,trueD,error_type,criterion,covtype,lasso_type)
%input
% D_isE=1 means penalty matrix "Upsilon_Delta" is all one matrix
%True D is true delta matrix
%error_type =c([max_abs_err, l1_err, matrix_l1_err, spectral_err, frobenius_err, nuclear_err];)
%           if error_type=1 means we use max_abs_error in BIC
% criterion=c("BIC","AIC","nBIC")
% lasso_type=0 means we are estimating the delta, >0 means we are estimating m_k
%covtype=1 means using kendall's tau to estimate Sigma, =2 means using sample covariance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%output
%loss is BIC norm loss
%score is BIC
%Dseq is sequence that stores the iteration of penalty matrix
%iter2 is sequence that stores the iteration of delta
%true_MatrixD is the penalty matrix(variance) using the true Delta
%distdelta is the frobeniuos form of estimated delta and true delta.
%TD is the set that stores the output [TD,TP,TN]
[nx,p]=size(X);
[ny,p]=size(Y);


if covtype==1
    hatcovMY = rankCovIID(Y);
    hatcovMX = rankCovIID(X);
else
    hatcovMY=cov(Y);
    hatcovMX=cov(X);
end

hatSigmagX=hatvargx(X,hatcovMX);
hatSigmagY=hatvargx(Y,hatcovMY);    
Dseq=[];
oldD=ones(p,p);
ED=ones(p,p);
true_MatrixD=varmatrix(hatcovMX,hatcovMY,hatSigmagX,hatSigmagY,trueD,2,nx,ny);
if D_isE==1
    %running the lasso with identity penalty matrix
    [hatdelta2,iter2]=L1_dts(hatcovMX,hatcovMY,rho,lambda,oldD,0);
    newD=oldD;
else
    iter2=[];
   
    hatdelta2=ED;
    for iter = 1:iternum
        newD=varmatrix(hatcovMX,hatcovMY,hatSigmagX,hatSigmagY,hatdelta2,iter,nx,ny);
        distD=distance(oldD,newD);
        Dseq=[Dseq,distD];
        % running the lasso with fixed penalty matrix
        [hatdelta2,iterseq]=L1_dts(hatcovMX,hatcovMY,rho,lambda,newD,lasso_type);
        oldD=newD;
        iter2=[iter2,iterseq];
        if (distD < tol_D)
            break;
        end
    end
end

%compute the BIC
[loss,score]=lasso_loss(hatdelta2,hatcovMY,hatcovMX, error_type,criterion,nx,ny);
%compute the TP
id1=(hatdelta2~=0);
id2=(trueD~=0); %trueD is truedelta
%TD
nom1=sum(sum(id1==id2 & id1==1));
denom1=sum(sum(id1));
%TP
nom2=sum(sum(id1==id2 & id1==1));
denom2=sum(sum(id2));
%TN
nom3=sum(sum(id1==id2 & id1==0));
denom3=sum(sum(id2==0));

distD=distance(true_MatrixD,newD);
distdelta=distance(trueD,hatdelta2);
TD=[nom1/denom1,nom2/denom2,nom3/denom3];
end


numRep=200;
n=100;

eStd = zeros(numRep, 1);
eStd2 = zeros(numRep, 1);
eStd21 = zeros(numRep, 1);
eStd3 = zeros(numRep, 1);
eStd4 = zeros(numRep, 1);
hatU = zeros(numRep, 1);
hatS0 = zeros(numRep, 1);
S0 = zeros(numRep, 1);
p1=10;
hatDeltavec = zeros(numRep, p1^2);
trueDeltavec = zeros(numRep, p1^2);
Ematrix=ones(10,10);
rho=0.4;
for rep=1:numRep
    %using rocket paper's genData function
    [X, precMX,covMX] = YT_genData(1, rho,n);
    [Y, precMY,covMY] = YT_genData(1, rho+0.1,n);
    trueDelta=precMY-precMX;
    tdelta=trueDelta;
    trueDeltavec(rep,:)=reshape(trueDelta,[],1)';
    p=size(trueDelta,1);
    onep=ones(p);
    trueTX=2/pi*asin(covMX);
    trueTY=2/pi*asin(covMY);
    Gamma=1/2*(kron(covMY,covMX)+kron(covMX,covMY));
    M=inv(Gamma);
    hatcovMX = rankCovIID(X);
    hatcovMY = rankCovIID(Y);
    hatTX=asin(hatcovMX) * 2 / pi;
    hatTY=asin(hatcovMY) * 2 / pi;  
    %get estimate delta and M
    D_isE=1; %using identity penalty matrix
    tol_D=10^(-8);
    iternum=1000;
    %just use lambda=0.127
    [Dseq4,iter4,hatDelta,TD4,true_MatrixD4,loss4,score4,distD4,distdelta4]=lasso_kendall(X,Y,iternum,0.127,1,tol_D,D_isE,tdelta,1,"BIC",1);
    S=(abs(trueDelta)>1e-8);
    hatS=(abs(hatDelta)>1e-8);
    [~,hatM2,tildeGamma] = oracleEstimator(hatS,hatcovMX,hatcovMY);
    hatM=zeros(p*p,p*p);
    for i=1:(p*p)
        [hatM_lasso_matrix,iterseq]=L1_dts(hatcovMX,hatcovMY,1,0.127,Ematrix,i);
        hatM(:,i)=reshape(hatM_lasso_matrix,[],1);
    end
    %[hatDelta,hatM,tildeGamma2] = oracleEstimator(S,covMX,covMY); use to get oracle estimator
    hatDeltavec(rep,:)=reshape(hatDelta,[],1)';
    
    %estimate sigma gx
    hatSigmagX=hatvargx(X,hatcovMX);
    hatSigmagY=hatvargx(Y,hatcovMY);
    
    hatmk=hatM(1,:)';
    mk=M(:,1);
    %hatU(rep)=gethatU(X,Y,hatmk,hatDelta);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %use case ii method to estimated variance
    [hatH_1_1_1,hatH_1_1_2,hatH_1_2_1,hatH_1_2_2,hatH_2_1,hatH_2_2,hatH_3_1,hatH_3_2] = getusefulmatrix(hatmk,hatDelta,hatcovMX,hatcovMY,hatTX,hatTY);
    % get S(m_k,delta,sigma_x,sigma_y)
    hatS0(rep)=getstatisticsS(X,Y,hatmk,hatDelta);
    S0(rep)=getstatisticsS(X,Y,mk,trueDelta);
    % estimated var hx and hy
    [hatvarhY,hatvarhX] = varkernelfunction(hatSigmagX,hatSigmagY,hatH_1_1_1,hatH_1_1_2,hatH_1_2_1,hatH_1_2_2,hatH_3_1,hatH_3_2);
    eStd(rep)=4/n*hatvarhX+4/n*hatvarhY;
    %using case i method to estimated variance with plug in the hatsigma,hatdelta,hatmk
    %get gradient S
    GSX=getgradientS(X,Y,hatmk,hatDelta,1);
    vecGSX=reshape(GSX,[],1);
    GSY=getgradientS(X,Y,hatmk,hatDelta,2);
    vecGSY=reshape(GSY,[],1);
    hatvarsX=vecGSX'*hatSigmagX*vecGSX;
    hatvarsY=vecGSY'*hatSigmagY*vecGSY;
    eStd2(rep)=4/n*hatvarsX+4/n*hatvarsY;
    %using case i method to estimated variance with truedelta truemk
    trueGSX=getgradientS(X,Y,mk,trueDelta,1);
    truevecGSX=reshape(trueGSX,[],1);
    trueGSY=getgradientS(X,Y,mk,trueDelta,2);
    truevecGSY=reshape(trueGSY,[],1);
    truevarsX=truevecGSX'*hatSigmagX*truevecGSX;
    truevarsY=truevecGSY'*hatSigmagY*truevecGSY;
    eStd3(rep)=4/n*truevarsX+4/n*truevarsY;
    
end

subplot(1,6,1)
qqplot((hatDeltavec(:,1)-trueDeltavec(:,1))./eStd)
subplot(1,6,2)
hist(eStd)
subplot(1,6,3)
hist(eStd2)
subplot(1,6,4)
hist(eStd3)
subplot(1,6,5)
qqplot((hatDeltavec(:,1)-trueDeltavec(:,1))./eStd2)
subplot(1,6,6)
hist(S0)

[mean(eStd), mean(eStd2), mean(eStd3), var(S0),var(hatS0)]

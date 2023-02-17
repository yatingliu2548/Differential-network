function [X, precM,covM] = YT_genData(modelType, rho,n)

% large grid example

% defaultStream = RandStream.getGlobalStream();
% reset(defaultStream, 1);

% outdir = '/home/mkolar/scratch-midway/copulaInf/simulation13/data';
% [~, ~, ~] = mkdir(outdir);

% set params
%n = 300; p = 100;
% numRep = 500; 

% covM = eye(p);
% covM(1,2) = rho; covM(2,1) = rho;
p=10;

covM = rho.^abs(repmat(1:p, p, 1) - repmat((1:p)', 1, p));
sqrtCov = sqrtm(covM);
precM = inv(covM);

if modelType == 2

    pArr = [10, 5, 5];

    a = [1, 2, 4];
    precM = [];
    covM = [];


    for k=1:3
        tmpPrecM = a(k) * eye(pArr(k)+2);
        for j=1:pArr(k)
            tmpPrecM(j, j+1) = (rho+0.1) * a(k);

            tmpPrecM(j+1, j) = (rho+0.1) * a(k);
        
            tmpPrecM(j, j+2) = 0.4  * a(k);
            tmpPrecM(j+2, j) = 0.4  * a(k);        
        end
        tmpPrecM = tmpPrecM(1:pArr(k),1:pArr(k));
    
        precM = blkdiag(precM, tmpPrecM);
        covM = blkdiag(covM, inv(tmpPrecM));
    end

    dCovM = diag( 1 ./ sqrt(diag(covM)) );
    corrM = dCovM * covM * dCovM;
    covM=corrM;
    precM = inv(corrM);
    sqrtCov = sqrtm(corrM);
    p = size(corrM, 1);
end



% precM = cl_sim13_genModel();
% sqrtCov = sqrtm(inv(precM));

% for modelType=1:6
%     for rep=1:numRep
%         X = generateData(modelType, sqrtCov, n, p);
%         save(sprintf('%s/data_%d.mat', outdir, (modelType-1)*numRep + rep), 'precM', 'X');
%     end
% end

X = generateData(modelType, sqrtCov, n, p);
end


function X = generateData(modelType, sqrtCov, n, p)

if modelType == 1
    X = normrnd(0, 1, [n, p]) * sqrtCov;
    X = marginalTranform(X);    %can be removed
elseif modelType == 2
    X = normrnd(0, 1, [n, p]) * sqrtCov;
    X = marginalTranform(X);    
elseif modelType == 3
    X = normrnd(0, 1, [n, p]) * sqrtCov;
    Y = chi2rnd(5, [n, 1]);
    X = repmat(sqrt(5 ./ Y), 1, p) .* X;    
elseif modelType == 4
    X = normrnd(0, 1, [n, p]) * sqrtCov;
    Y = chi2rnd(5, [n, 1]);
    X = repmat(sqrt(5 ./ Y), 1, p) .* X;
    X = marginalTranform(X);    
elseif modelType == 5
    X = normrnd(0, 1, [n, p]);
    X = repmat(sqrt(1./sum(X.*X, 2)), 1, p) .* X;
    X = X * sqrtCov;
    sc= abs(trnd(1.5,n,1));
    X = repmat(sc, 1, p) .* X;
elseif modelType == 6
    X = normrnd(0, 1, [n, p]);
    X = repmat(sqrt(1./sum(X.*X, 2)), 1, p) .* X;
    X = X * sqrtCov;
    sc= abs(trnd(1.5,n,1));
    X = repmat(sc, 1, p) .* X;    
    X = marginalTranform(X);    
end

    
end

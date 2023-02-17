function S = getstatisticsS(X,Y,mk,Delta)

hatcovMX = rankCovIID(X);
hatcovMY = rankCovIID(Y);
   
hatGamma=1/2*(kron(hatcovMX,hatcovMY)+kron(hatcovMY,hatcovMX));
S=mk'*(hatGamma*reshape(Delta,[],1)-reshape(hatcovMX-hatcovMY,[],1));
end

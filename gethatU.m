function hatU = gethatU(X,Y,hatmk,hatDelta)
p=size(X,2);
onep=ones(p);
hatcovMX = rankCovIID(X);
hatcovMY = rankCovIID(Y);
hatTX=asin(hatcovMX) * 2 / pi;
hatTY=asin(hatcovMY) * 2 / pi;

[hatH_1_1_1,hatH_1_1_2,hatH_1_2_1,hatH_1_2_2,hatH_2_1,hatH_2_2,hatH_3_1,hatH_3_2] = getusefulmatrix(hatmk,hatDelta,hatcovMX,hatcovMY,hatTX,hatTY);
hatU=trace(hatH_1_1_1'*kron(onep,hatTX))+trace(hatH_1_1_2'*kron(hatTX,onep))+trace(hatH_1_2_1'*kron(onep,hatTY))+trace(hatH_1_2_2'*kron(hatTY,onep))+hatH_3_1'*reshape(hatTX,[],1)-hatH_3_2'*reshape(hatTY,[],1);

end
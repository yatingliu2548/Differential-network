function hatdiff = gethatdiff(X,Y,a,b)
p=size(X,2);

hatcovMX = rankCovIID(X);
hatcovMY = rankCovIID(Y);
if a>=1
    hatcovMX=hatcovMX(a,b);
    hatcovMY=hatcovMY(a,b);
end
hatdiff=hatcovMX-hatcovMY;
end
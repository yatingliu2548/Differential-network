function vargx = hatvargx(X,covM)

[n, p] = size(X);

mu=reshape(covM,[],1);
tempcov=zeros(p^2);
for a=1:n
    t = sign(repmat(X(a, :), (n-1), 1) - X(setdiff(1:n,a), :)); 
    temp=reshape(t'*t,[],1)/(n-1)-mu;
    tempcov = tempcov+temp*temp';
end
vargx = tempcov / n;


end


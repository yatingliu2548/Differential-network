function result= genp(p,sparsity,sigma)


Delta= zeros(p, p);
for i=1:(p-1)
    for j=(i+1):p
        if unifrnd(0,1)<sparsity
            Delta(i,j)=sign(normrnd(0,1))*sigma;
        end
    end
end

Delta = Delta+Delta.';
[U,meig] = eig(Delta);
[meig, idxx] = sort(diag(meig), 'descend');% Dx= is vector
U = U(:, idxx);
result = Delta-(meig(p)-1)*diag(ones(p,1));

end


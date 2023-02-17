function [loss,Score]=lasso_loss(D,SigmaY,SigmaX, type,criterion,nx,ny)
%D is hatD
S1=SigmaX;
S0=SigmaY;
err = (S1 * D * S0 + S0 * D * S1) / 2 - S1 + S0;
err_abs = abs(err);
max_abs_err = max(err_abs(:));
l1_err = sum(err_abs(:));
matrix_l1_err = max(sum(err_abs, 2));
p=size(D,1);
[~, S, ~] = svd(err);
singular_values = sort(diag(S), 'descend');
spectral_err = singular_values(1);
frobenius_err = norm(err, 'fro');
nuclear_err = sum(singular_values);
error = [max_abs_err, l1_err, matrix_l1_err, spectral_err, frobenius_err, nuclear_err];
loss=error(type);
%lowertri=D(find(tril(D)));
[s1,s2]=size(D(D~=0));
numnonzero=s1*s2/2+p;
if criterion=="BIC"
    Score=(nx+ny)*loss+log(nx+ny)*numnonzero;
elseif criterion=="AIC"
    Score=(nx+ny)*loss+2*numnonzero;
else
    Score=(nx+ny)*loss+2*log(nx+ny)*numnonzero/(nx+ny);
end

end

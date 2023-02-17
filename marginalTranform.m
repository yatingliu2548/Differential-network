function [Y] = marginalTranform(X)

% f_1(x)=x, 
% f_2(x)=sign(x)|x|^{1/2}
% f_3(x)=x^3
% f_4=\Phi(x) 
% f_5(x)=\exp(x)

f = { ...
      @(x) x, ...
      @(x) sign(x) .* sqrt(abs(x)), ...
      @(x) x.^3, ...
      @(x) normcdf(x), ...
      @(x) exp(x) ...
    };

[n, p] = size(X);

Y = zeros(n, p);

for i=1:p
    fInd = mod((i-1), 5) + 1;
    Y(:, i) = f{fInd}(X(:,i));
end


end
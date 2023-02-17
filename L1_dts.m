function [Delta3,itererror] = L1_dts(SigmaX,SigmaY,rho,lambda,Gamma,type)
% Gamma is matrix of p*p instead of diagonal matrix
%Gamma is penalty matrix
%type>0 means we are doing lasso for m_k type=k

[n1,p1]=size(SigmaX); %n1=p1=p
[n2,p2]=size(SigmaY);
Delta0=inv(SigmaY+diag(ones(n1,1)))-inv(SigmaX+diag(ones(n1,1)));
Lambda0=zeros(n1,p1);
tol = 1e-5;
[n,p]= size(Delta0);
[Ux,Dx] = eig(SigmaX);

[Dx, idxx] = sort(diag(Dx), 'descend');% Dx= is vector
Ux = Ux(:, idxx);

[Uy,Dy] = eig(SigmaY);
[Dy, idxy] = sort(diag(Dy), 'descend');% Dx= is vector
Uy = Uy(:, idxy);


D1 = zeros(p,p);
D2 = zeros(p,p);
for i=1:p
    for j=1:p
        D1(i,j)=1/(Dy(i)*Dx(j)+4*rho); % for Delta 1
        D2(i,j)= 1/(Dy(j)*Dx(i)+4*rho); % for Delta 2
    end
end
%initialization
 Delta1 = Delta0;
 Delta2 = Delta0;
 Delta3 = Delta0;
 Lambda1 = Lambda0;
 Lambda2 = Lambda0;
 Lambda3 = Lambda0;
 % given delta, lambda
 tUy = Uy.';
 tUx = Ux.';
 itererror=[];
 for k =1:1000
     if type>0
         indi=mod(type-1,p)+1;
         indj=ceil(type/p);
         Ematrix=zeros(p,p);
         Ematrix(indi,indj)=1;
         C1=Ematrix + 2*rho * Delta2 + 2*rho * Delta3 + 2*Lambda1 - 2*Lambda3;
     else
         C1 = SigmaX - SigmaY + 2*rho * Delta2 + 2*rho * Delta3 + 2*Lambda1 - 2*Lambda3;
     end
     z1 = Ux*(D1.*(tUx * C1 * tUy))*Uy; % 
     
     if type>0
         indi=mod(type-1,p)+1;
         indj=ceil(type/p);
         Ematrix=zeros(p,p);
         Ematrix(indi,indj)=1;
         C2=Ematrix + 2*rho * z1 + 2*rho * Delta3 + 2*Lambda3 - 2*Lambda2;
     else
         C2 = SigmaX - SigmaY + 2*rho * z1 + 2*rho * Delta3 + 2*Lambda3 - 2*Lambda2;
     end
     
     z2 = Uy*(D2.*(tUy * C2 * tUx))*Ux; % 

     A = 1/(2*rho)*(rho*(z1 + z2) + Lambda2 - Lambda1);
     
     idx1= A>lambda;
     idx2= A<-lambda;
     idx3= A<=lambda & A>= -lambda;
     idx4= Gamma==0;
     A(idx1)=A(idx1)-lambda*sqrt(Gamma(idx1));
     A(idx2)=A(idx2)+lambda*sqrt(Gamma(idx2));
     A(idx3)=0;
     A(idx4)=0;

     z3=A;

     % Calculate dis1, dis2, and dis3
     a = sum(Delta1(:).^2); 
     b = sum(z1(:).^2);
     c = sum((Delta1(:) - z1(:)).^2);
     dis1 = c / max(1, min(a, b));

     a = sum(Delta2(:).^2);
     b = sum(z2(:).^2);
     c = sum((Delta2(:) - z2(:)).^2);
     dis2 = c / max(1, min(a, b));
     a = sum(z1(:).^2);
     b = sum(z2(:).^2);
     c = sum((z1(:) - z2(:)).^2);
     dis3 = c / max(1, min(a, b));

     itererror=[itererror,dis3];
     % Check if dis1, dis2, and dis3 are less than tol
     if (dis1 < tol && dis2 < tol && dis3 < tol)
         break;
     end 

     % Update Delta, Lambda
     Delta1 = z1;
     Delta2 = z2;
     Delta3 = z3;
     Lambda1 = Lambda1 + rho * (Delta3 - Delta1);
     Lambda2 = Lambda2 + rho * (Delta2 - Delta3);
     Lambda3 = Lambda3 + rho * (Delta1 - Delta2);
end
end
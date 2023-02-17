function [varhy,varhx] = varkernelfunction(sigmagX,sigmagY,H_1_1_1,H_1_1_2,H_1_2_1,H_1_2_2,H_3_1,H_3_2)
x_3_1=H_3_1;
x_3_2=H_3_2;
x_1_1_1=reshape(H_1_1_1,[],1);
x_1_1_2=reshape(H_1_1_2,[],1);
x_1_2_1=reshape(H_1_2_1,[],1);
x_1_2_2=reshape(H_1_2_2,[],1);
psquare=size(H_1_1_1,1);
vecsigmagx=reshape(sigmagX,[],1);
vecsigmagy=reshape(sigmagY,[],1);
one=ones(psquare,1);
onem=one*one.';

varhx=(x_1_1_1')*(kron(onem,sigmagX))*x_1_1_1+(x_1_1_2')*(kron(sigmagX,onem))*x_1_1_2+2*(kron(kron(one,vecsigmagx),one)')*reshape(x_1_1_1*(x_1_1_2'),[],1)+x_3_1'*sigmagX*x_3_1+2*kron(one,vecsigmagx)'*reshape(x_1_1_1*(x_3_1'),[],1)+2*kron(vecsigmagx,one)'*reshape(x_3_1*(x_1_1_2'),[],1);
varhy=(x_1_2_1')*(kron(onem,sigmagY))*x_1_2_1+(x_1_2_2')*(kron(sigmagY,onem))*x_1_2_2+2*(kron(kron(one,vecsigmagy),one)')*reshape(x_1_2_1*(x_1_2_2'),[],1)+x_3_2'*sigmagY*x_3_2-2*kron(one,vecsigmagy)'*reshape(x_1_2_1*(x_3_2'),[],1)-2*kron(vecsigmagy,one)'*reshape(x_3_2*(x_1_2_2'),[],1);

end
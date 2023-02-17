function [H_1_1_1,H_1_1_2,H_1_2_1,H_1_2_2,H_2_1,H_2_2,H_3_1,H_3_2] = getusefulmatrix(mk,Delta,covMX,covMY,TX,TY)

vecDelta=reshape(Delta,[],1);
H_1_1_1=(mk*vecDelta')'.*(kron(covMY,pi/2*cos(pi/2*TX)));
H_1_1_2=(mk*vecDelta')'.*(kron(pi/2*cos(pi/2*TX),covMY));
H_1_2_1=(mk*vecDelta')'.*(kron(covMX,pi/2*cos(pi/2*TY)));
H_1_2_2=(mk*vecDelta')'.*(kron(pi/2*cos(pi/2*TY),covMX));
H_3_1=mk.*reshape(pi/2*cos(pi/2*TX),[],1);
H_3_2=mk.*reshape(pi/2*cos(pi/2*TY),[],1);

H_2_1=(mk*vecDelta')'.*kron(pi/2*cos(pi/2*TX),pi/2*cos(pi/2*TY));
H_2_2=(mk*vecDelta')'.*kron(pi/2*cos(pi/2*TY),pi/2*cos(pi/2*TX));


end

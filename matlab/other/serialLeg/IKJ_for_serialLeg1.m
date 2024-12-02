function [tauf,taub] = IKJ_for_serialLeg1(Llinks,Fx,Fz,anglef,angleb)

% temp1 = Llinks(1)*sin(anglef-angleb);
% temp2 = Llinks(2)*sin(anglef-angleb);
% 
% JI = [-sin(anglef)/temp1,-cos(anglef)/temp1;sin(angleb)/temp2,cos(angleb)/temp2];
J = [-Llinks(1)*cos(angleb),-Llinks(2)*cos(anglef);Llinks(1)*sin(angleb),Llinks(2)*sin(anglef)];
F = [Fx;Fz];
tau = J'*F;
tauf = tau(2);
taub = tau(1);
end
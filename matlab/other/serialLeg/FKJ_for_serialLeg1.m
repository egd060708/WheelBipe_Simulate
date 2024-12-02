function [Fx,Fz] = FKJ_for_serialLeg1(Llinks,tauf,taub,anglef,angleb)

J = [-Llinks(1)*cos(angleb),-Llinks(2)*cos(anglef);Llinks(1)*sin(angleb),Llinks(2)*sin(anglef)];
tau = [taub;tauf];
F = J*tau;
Fx = F(1,1);
Fz = F(2,1);
end
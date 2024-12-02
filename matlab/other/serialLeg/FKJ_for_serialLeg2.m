function [Fx,Fz] = FKJ_for_serialLeg2(Llinks,tauf,taub,anglef,angleb)
a1 = angleb;
a2 = anglef;
x3 = -Llinks(4)*sin(a2);
z3 = -Llinks(4)*cos(a2);
x2 = -Llinks(1)*sin(a1);
z2 = -Llinks(1)*cos(a1);

L6 = sqrt((x3-x2).^2+(z3-z2).^2);
a1p = acos((L6.^2+Llinks(1).^2-Llinks(2).^2)/(2*L6*Llinks(1)));
a4p = acos((L6.^2+Llinks(4).^2-Llinks(3).^2)/(2*L6*Llinks(4)));

a4 = a1+a1p+a4p-pi;

dL6da1 = Llinks(1)*Llinks(4)/L6*sin(a2-a1);
dL6da2 = -Llinks(1)*Llinks(4)/L6*sin(a2-a1);
da1pda1 = -((L6.^2-Llinks(1).^2+Llinks(2).^2)/(2*sin(a1p)*Llinks(1)*L6.^2))*dL6da1;
da4pda1 = -((L6.^2-Llinks(4).^2+Llinks(3).^2)/(2*sin(a1p)*Llinks(4)*L6.^2))*dL6da1;
da1pda2 = -((L6.^2-Llinks(1).^2+Llinks(2).^2)/(2*sin(a1p)*Llinks(1)*L6.^2))*dL6da2;
da4pda2 = -((L6.^2-Llinks(4).^2+Llinks(3).^2)/(2*sin(a1p)*Llinks(4)*L6.^2))*dL6da2;

da4da1 = 1 + da1pda1 + da4pda1;
da4da2 = 1 + da1pda2 + da4pda2;

dxda1 = -Llinks(1)*cos(a1)-(Llinks(4)+Llinks(5))*cos(a4)*da4da1;
dxda2 = -(Llinks(4)+Llinks(5))*cos(a4)*da4da2;
dzda1 = Llinks(1)*sin(a1)+(Llinks(4)+Llinks(5))*sin(a4)*da4da1;
dzda2 = (Llinks(4)+Llinks(5))*sin(a4)*da4da2;
J = [dxda1,dxda2;dzda1,dzda2];
tau = [taub,tauf];
F = (J')\tau;
Fx = F(1);
Fz = F(2);



end
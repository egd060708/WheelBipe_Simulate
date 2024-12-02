function [x,z] = FK_for_serialLeg2(Llinks,anglef,angleb)
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
x = -Llinks(1)*sin(a1)-(Llinks(4)+Llinks(5))*sin(a4);
z = -Llinks(1)*cos(a1)-(Llinks(4)+Llinks(5))*cos(a4);

end
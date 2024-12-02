function [x,z] = FK_for_serialLeg1(Llinks,anglef,angleb)
% 串联腿第一类正运动学建模
x = - Llinks(1)*sin(angleb) - Llinks(2)*sin(anglef);
z = - Llinks(1)*cos(angleb) - Llinks(2)*cos(anglef);


end
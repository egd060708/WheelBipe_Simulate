function [P,M,I] = leg_param_cal(Lleg,Llinks,Ilinks,Mlinks)
% 输入Llinks三段连杆数组，Lleg当前腿长
% Ilinks后两段连杆在其重心处的转动惯量，Mlinks后两段连杆重量
[alpha1,alpha2] = inverse_kinematics_solution(0,-Lleg,Llinks);
P1 = [0.5*(Llinks(1)+Llinks(2)*cos(alpha1)), 0.5*Llinks(2)*sin(alpha1)];
P2 = [0.5*(0.5*Llinks(1)+Llinks(2)*cos(alpha1)), 0.5*(Llinks(2)*sin(alpha1)+Lleg)];
Pm = (Mlinks(1)*P1(2)+Mlinks(2)*P2(2))/(Mlinks(1)+Mlinks(2));
P = Lleg - Pm;
I = 2*(Ilinks(1)+Mlinks(1)*(P1(1)^2+(P1(2)-Pm)^2)+Ilinks(2)+Mlinks(2)*(P2(1)^2+(P2(2)-Pm)^2));
M = 2*Mlinks(1)+2*Mlinks(2);
end
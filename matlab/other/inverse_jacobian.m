function [Fx,Fy]=inverse_jacobian(tau1,tau2,a1,a2)
%逆向雅可比求解
[L1,L2,L3,L4,L5]=deal(0.15,0.288,0.288,0.15,0.15);
A = (L1.^2)*((sin(a2)-sin(a1)).^2)+(L5+L1*(cos(a1)-cos(a2))).^2;
B = -2*(L1.^2)*(sin(a2)-sin(a1))*cos(a1)-2*(L5+L1*(cos(a1)-cos(a2)))*L1*sin(a1);
C = 2*(L1.^2)*(sin(a2)-sin(a1))*cos(a2)+2*(L5+L1*(cos(a1)-cos(a2)))*L1*sin(a2);
D = sqrt((L2.^2)/A-0.25);
E = (L1*(sin(a1)-sin(a2))*(L2.^2))/(2*(A.^2)*D);
F = ((L5+L1*(cos(a1)-cos(a2)))*(L2.^2))/(2*(A.^2)*D);

x_a1 = E*B-L1*cos(a1)*D-0.5*L1*sin(a1);
x_a2 = E*C+L1*cos(a2)*D-0.5*L1*sin(a2);
y_a1 = F*B+L1*sin(a1)*D-0.5*L1*cos(a1);
y_a2 = F*C-L1*sin(a2)*D-0.5*L1*cos(a2);

J = [x_a1,x_a2;...
    y_a1,y_a2;];
tau = [tau1,tau2];
F = (J')\(tau');%左乘转置矩阵的倒数
Fx = F(1,1);
Fy = F(2,1);

end
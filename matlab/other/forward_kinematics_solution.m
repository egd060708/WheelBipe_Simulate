function [x,y] = forward_kinematics_solution(Alpha1, Alpha2)
[L1,L2,L3,L4,L5]=deal(0.15,0.288,0.288,0.15,0.15);
X3 = L5/2+L1*cos(Alpha1);
Y3 = -L1*sin(Alpha1);
X4 = -L5/2+L4*cos(Alpha2);
Y4 = -L4*sin(Alpha2);
X43 = X3-X4;
Y43 = Y3-Y4;
X6 = (X3+X4)/2;
Y6 = (Y3+Y4)/2;

D43 = sqrt(X43^2+Y43^2);

x = sqrt((L2^2-D43^2/4)/D43^2)*(Y43)+X6;
y = sqrt((L3^2-D43^2/4)/D43^2)*(-X43)+Y6;

end
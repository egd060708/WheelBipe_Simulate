clc;
clear;

% syms L1 L2 theta1 theta2 real
% J = [-L1*cos(theta1),-L2*cos(theta2);L1*sin(theta1),L2*sin(theta2)];
% Ji = inv(J);
% pretty(Ji);

eqLength = 0.13:0.01:0.4;
num = size(eqLength,2);
thigh = zeros(1,num);
knee = zeros(1,num);
Llinks = [0.212,0.245];
for i=1:1:num
    [thigh(1,i),knee(1,i)] = IK_for_serialLeg(Llinks,eqLength(i),0);
end

figure;
plot(eqLength,thigh,eqLength,knee);
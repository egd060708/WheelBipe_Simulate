clc;
clear;

% syms L1 L2 theta1 theta2 real
% J = [-L1*cos(theta1),-L2*cos(theta2);L1*sin(theta1),L2*sin(theta2)];
% Ji = inv(J);
% pretty(Ji);

% eqLength = 0.13:0.01:0.4;
% num = size(eqLength,2);
% thigh = zeros(1,num);
% knee = zeros(1,num);
% Llinks = [0.212,0.245];
% for i=1:1:num
%     [thigh(1,i),knee(1,i)] = IK_for_serialLeg(Llinks,eqLength(i),0);
% end
% 
% figure;
% plot(eqLength,thigh,eqLength,knee);



% [af,ab] = inverse_kinematics_solution(0,-0.3);
% disp(af);
% disp(ab);
% [tf,tb] = inverse_jacobian(0,-10,af,ab);
% disp(tf);
% disp(tb);
% 
% [ab,af] = IK_for_serialLeg1(Llinks,0.3,0);
% disp(af);
% disp(ab);
% [tf,tb] = IKJ_for_serialLeg1(Llinks,0,-10,af,ab);
% disp(tf);
% disp(tb);

Llinks1 = [0.212,0.245];
Llinks2 = [0.212,0.07,0.212,0.07,0.245-0.07];

% [x,z] = FK_for_serialLeg1(Llinks1,-1,0.5);
% disp(x);
% disp(z);
% 
% [x,z] = FK_for_serialLeg2(Llinks2,-1,0.5);
% disp(x);
% disp(z);
% 
% [anglef,angleb] = IK_for_serialLeg1(Llinks1,0.3,0.5);
% disp(anglef);
% disp(angleb);
% 
% [anglef,angleb] = IK_for_serialLeg2(Llinks2,0.3,0.5);
% disp(anglef);
% disp(angleb);

% [af,ab] = IK_for_serialLeg1(Llinks1,0.3,0);
% disp(af);
% disp(ab);
% [tf,tb] = IKJ_for_serialLeg1(Llinks1,0,-10,af,ab);
% disp(tf);
% disp(tb);
% 
% [af,ab] = IK_for_serialLeg2(Llinks2,0.3,0);
% disp(af);
% disp(ab);
% [tf,tb] = IKJ_for_serialLeg2(Llinks2,0,-10,af,ab);
% disp(tf);
% disp(tb);

eqLength = 0.13:0.01:0.4;
num = size(eqLength,2);
af = zeros(1,num);
ab = zeros(1,num);
tf = zeros(1,num);
tb = zeros(1,num);
for i=1:1:num
    [af(1,i),ab(1,i)] = IK_for_serialLeg1(Llinks1,eqLength(i),1);
    [tf(1,i),tb(1,i)] = IKJ_for_serialLeg1(Llinks1,0,-10,af(1,i),ab(1,i));
end

figure;
plot(eqLength,af,eqLength,ab);
figure;
plot(eqLength,tf,eqLength,tb);

for i=1:1:num
    [af(1,i),ab(1,i)] = IK_for_serialLeg2(Llinks2,eqLength(i),1);
    [tf(1,i),tb(1,i)] = IKJ_for_serialLeg2(Llinks2,0,-10,af(1,i),ab(1,i));
end

figure;
plot(eqLength,af,eqLength,ab);
figure;
plot(eqLength,tf,eqLength,tb);
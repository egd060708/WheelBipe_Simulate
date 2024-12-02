%轮足五连杆工作空间求解

clear,clc;
% i = 1;
% ang_min = -30;%设置轮腿角度范围
% ang_max = 80;
% for a1 = ang_min/180*pi:0.01:ang_max/180*pi%离散化
%     for a2 = ang_min/180*pi:0.01:ang_max/180*pi
%         [x(i),y(i)] = forward_kinematics_solution(a1,pi-a2);%正运动学映射
%         i = i + 1;
%     end
% end
% 
% figure;     % 绘制工作空间
% scatter(x,y);

x = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16];
y = [4.46 8.76 11.57 16.65 22.57 24.62 29.18 32.91 37.19 41.75 45.93 50.36 54.42 59.08 63.15 68.10];
fit_feedback(x,y,1);

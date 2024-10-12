%% 系统建模
clc;
clear;

% 获取当前脚本的完整路径
currentFile = mfilename('fullpath');
% 使用fileparts分离出路径部分
[currentPath,~,~] = fileparts(currentFile);
load(strcat(currentPath,'/BalanceTurn_Model.mat'));
syms t real; %定义微分时间
syms Twl Twr Tjl Tjr; %定义输入力矩，[左右轮，左右关节]
syms Mw R Jw D real; %定义驱动轮相关物理量，[轮质量，半径，惯量，轮距]
syms Ml Mr Hll Hlr Jll Jlr Pll Plr real; %定义腿杆相关物理量，[质量，长度，惯量，重心离轮轴位置]
syms Mb Pb Jbpitch Jbyaw real; %定义机体相关物理量，[质量，重心离关节位置，pitch惯量，yaw惯量]
syms thetawl(t) thetawr(t) thetall(t) thetalr(t) thetab(t); %定义广义坐标系，[位移，旋转角，左腿倾角，右腿倾角，机体倾角]
syms g real; %定义重力加速度

% 导入物理参数
webots_param;
% 求解部分物理参数
Lleg = 0.3;
[Pleg,Mleg,Ileg] = leg_param_cal(Lleg,Llinks,Ilinks,Mlinks);
Ml_= Mleg;%kg
Mr_= Mleg;
Hll_= Lleg;
Hlr_= Lleg;
Jll_= Ileg;
Jlr_= Ileg;
Pll_= Pleg;
Plr_= Pleg;

% 被替代的参数
symbolic_param = [Mw R Jw D Ml Mr Hll Hlr Jll Jlr Pll Plr Mb Pb Jbpitch Jbyaw g];
% 用于替代的实参
real_param = [Mw_ R_ Jw_ D_ Ml_ Mr_ Hll_ Hlr_ Jll_ Jlr_ Pll_ Plr_ Mb_ Pb_ Jbpitch_ Jbyaw_ g_];

%参数实例化
M=subs(M,symbolic_param,real_param);
M=vpa(M,7);%规范小数位
G=subs(G,symbolic_param,real_param);
G=vpa(G,7);
J=subs(J,R,R_);
%传递函数
H=M \ (-G);
I=M \ J;
H=double(H);
I=double(I);

A_ss = [0,1,0,0,0,0,0,0,0,0;...
    H(1,1),0,H(1,2),0,H(1,3),0,H(1,4),0,H(1,5),0;...
    0,0,0,1,0,0,0,0,0,0;...
    H(2,1),0,H(2,2),0,H(2,3),0,H(2,4),0,H(2,5),0;...
    0,0,0,0,0,1,0,0,0,0;...
    H(3,1),0,H(3,2),0,H(3,3),0,H(3,4),0,H(3,5),0;...
    0,0,0,0,0,0,0,1,0,0;...
    H(4,1),0,H(4,2),0,H(4,3),0,H(4,4),0,H(4,5),0;...
    0,0,0,0,0,0,0,0,0,1;...
    H(5,1),0,H(5,2),0,H(5,3),0,H(5,4),0,H(5,5),0];
B_ss = [0,0,0,0;...
    I(1,1),I(1,2),I(1,3),I(1,4);...
    0,0,0,0;...
    I(2,1),I(2,2),I(2,3),I(2,4);...
    0,0,0,0;...
    I(3,1),I(3,2),I(3,3),I(3,4);...
    0,0,0,0;...
    I(4,1),I(4,2),I(4,3),I(4,4);...
    0,0,0,0;...
    I(5,1),I(5,2),I(5,3),I(5,4)];
B_ss = double(B_ss);
C_ss = eye(10);
D_ss = zeros(10,4);

%系统离散化
sys_c = ss(A_ss, B_ss, C_ss, D_ss);
sys_d = c2d(sys_c, Ts);

%保存连续状态空间矩阵
save('E:\Git_Project\gitea\wheelbipe_simulate\轮足\matlab\Wheelbipe_cState.mat','sys_c');
%保存离散化状态空间矩阵
save("E:\Git_Project\gitea\wheelbipe_simulate\轮足\matlab\Wheelbipe_dState.mat",'sys_d');
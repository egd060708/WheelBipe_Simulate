%% 系统建模
clc;
clear;

Ts = 0.002;%离散时间

% 获取当前脚本的完整路径
currentFile = mfilename('fullpath');
% 使用fileparts分离出路径部分
[currentPath,~,~] = fileparts(currentFile);
% load(strcat(currentPath,'/../model/BalanceTurn_Model.mat'));
load('D:/Git_Project/github/WheelBipe_Simulate/matlab/model/BalanceTurn_Model_Simple.mat')
syms t real; %定义微分时间
syms Twl Twr Tjl Tjr; %定义输入力矩，[左右轮，左右关节]
syms Mw R Jw D real; %定义驱动轮相关物理量，[轮质量，半径，惯量，轮距]
syms Ml Mr Hll Hlr Jll Jlr Pll Plr real; %定义腿杆相关物理量，[质量，长度，惯量，重心离轮轴位置]
syms Mb Pb Jbpitch Jbyaw real; %定义机体相关物理量，[质量，重心离关节位置，pitch惯量，yaw惯量]
% syms thetawl(t) thetawr(t) thetall(t) thetalr(t) thetab(t); %定义广义坐标系，[位移，旋转角，左腿倾角，右腿倾角，机体倾角]
syms x(t) fi(t) thetall(t) thetalr(t) thetab(t); %定义广义坐标系，[位移，旋转角，左腿倾角，右腿倾角，机体倾角]
syms g real; %定义重力加速度

% 导入物理参数
webots_param;
robot_type = 1;% 0为并联腿，1为串联腿
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
X=subs(M,symbolic_param,real_param);
X=vpa(X,7);%规范小数位
Y=subs(G,symbolic_param,real_param);
Y=vpa(Y,7);
Z=subs(J,[R D],[R_ D_]);
Z=vpa(Z,7);
%传递函数
H=X \ (-Y);
I=X \ Z;
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


disp(sys_d.A);
disp(sys_d.B);

%% 单层lqr
%clc;
%判断可控性
if (rank(ctrb(sys_d.A,sys_d.B))==10)
    lqr_Q=eye(10);
    lqr_R=eye(4);
    
    % 转向模型参数
    % lqr_Q(1,1) = 0.00001;
    % lqr_Q(2,2) = 10;
    % lqr_Q(3,3) = 0.00001;
    % lqr_Q(4,4) = 10;
    % lqr_Q(5,5) = 10;
    % lqr_Q(6,6) = 75;
    % lqr_Q(7,7) = 10;
    % lqr_Q(8,8) = 75;
    % lqr_Q(9,9) = 5000;
    % lqr_Q(10,10) = 1;
    % 
    % lqr_R(1,1) = 1;
    % lqr_R(2,2) = 1;
    % lqr_R(3,3) = 1;
    % lqr_R(4,4) = 1;

    % 简化转向模型的参数
    % lqr_Q(1,1) = 0.00001;
    % lqr_Q(2,2) = 50;
    % lqr_Q(3,3) = 0.00001;
    % lqr_Q(4,4) = 25;
    % lqr_Q(5,5) = 0.5;
    % lqr_Q(6,6) = 0.05;
    % lqr_Q(7,7) = 0.5;
    % lqr_Q(8,8) = 0.05;
    % lqr_Q(9,9) = 20000;
    % lqr_Q(10,10) = 0.5;
    % 
    % lqr_R(1,1) = 7;
    % lqr_R(2,2) = 7;
    % lqr_R(3,3) = 1;
    % lqr_R(4,4) = 1;

    %调试mpc参数
    % rho = 0.;
    % lqr_Q(1,1) = rho + 0.00001;
    % lqr_Q(2,2) = rho + 100;
    % lqr_Q(3,3) = rho + 0.00001;
    % lqr_Q(4,4) = rho + 25;
    % lqr_Q(5,5) = rho + 200;
    % lqr_Q(6,6) = rho + 1.;
    % lqr_Q(7,7) = rho + 200;
    % lqr_Q(8,8) = rho + 1.;
    % lqr_Q(9,9) = rho + 10000;
    % lqr_Q(10,10) = rho + 0.5;
    % 
    % lqr_R(1,1) = rho + 10;
    % lqr_R(2,2) = rho + 10;
    % lqr_R(3,3) = rho + 1;
    % lqr_R(4,4) = rho + 1;

    %实车
    lqr_Q(1,1) =   0.00000000000000001;
    lqr_Q(2,2) =   80;
    lqr_Q(3,3) =   0.2;
    lqr_Q(4,4) =   5;
    lqr_Q(5,5) =   400;
    lqr_Q(6,6) =   3;
    lqr_Q(7,7) =   400;
    lqr_Q(8,8) =   3;
    lqr_Q(9,9) =   8000;
    lqr_Q(10,10) = 1.7;

    lqr_R(1,1) = 3;
    lqr_R(2,2) = 3;
    lqr_R(3,3) = 0.1;
    lqr_R(4,4) = 0.1;
    
    K=dlqr(sys_d.A,sys_d.B,lqr_Q,lqr_R);
else
    disp('Uncontrollable!');
end
fprintf("\n");

% Init_Pos = [0;0;30/180*pi;0;30/180*pi;0];
% Tw_max = 4.7;
% Tg_max = 20;

%% 拟合lqr
clc;
%设定轮杆长度
h_top=0.4;
h_bottom=0.13;
h_step=0.01;
Ts = 0.005;

%取出拟合lqr参数
%[K_s,H_s,C_s] = lqr_fit (Ts,h_top,h_bottom,h_step,lqr_Q,lqr_R,strcat(currentPath,'/../model/BalanceTurn_Model.mat'));
[K_s,P_s,H_s,C_s] = lqr_fit (Ts,h_top,h_bottom,h_step,lqr_Q,lqr_R,'D:/Git_Project/github/WheelBipe_Simulate/matlab/model/BalanceTurn_Model_Simple.mat',robot_type);

xNum = size(lqr_Q,1);%获取行数，也就是状态变量个数
uNum = size(lqr_R,2);%获取列数，也就是输入变量个数
Kp = cell(uNum,xNum);%构建参数元组
Pp = cell(xNum,xNum);

fprintf("MPCFloat model_K[4*4*10]={\n");
for j=1:1:uNum
    for k=1:1:xNum
        for i=C_s:-1:1
            Kp{j,k}(i)=K_s{1,i}(j,k);
        end
        % lqr_fit_feedback(H_s,Kp{j,k},"/*K参数*/",strcat("leg_lqr_params[",num2str(j-1),"][",num2str(k-1),"]"));
        F = fit_feedback(H_s,Kp{j,k},3);
        fprintf(sprintf("%d,%d,%d,%d,\n",F(4),F(3),F(2),F(1)));%打印低次项在前
    end
end
fprintf("};\n");
fprintf("MPCFloat model_P[4*10*10]={\n");
for j=1:1:xNum
    for k=1:1:xNum
        for i=C_s:-1:1
            Pp{j,k}(i)=P_s{1,i}(j,k);
        end
        % figure;
        % scatter(H_s, Pp{j,k}, 40, 'b', 'filled');  % 'filled'表示填充点
        % title('参数散点图');
        % xlabel('h');
        % ylabel('p');
        % grid on;
        % hold on;
        % coefficients = polyfit(H_s, Pp{j,k}, 3);
        % fittedY = polyval(coefficients, H_s);
        % plot(H_s, fittedY, 'r-', 'LineWidth', 2);

        % lqr_fit_feedback(H_s,Pp{j,k},"/*P参数*/",strcat("leg_lqrP_params[",num2str(j-1),"][",num2str(k-1),"]"));
        F = fit_feedback(H_s,Pp{j,k},3);
        fprintf(sprintf("%d,%d,%d,%d,\n",F(4),F(3),F(2),F(1)));%打印低次项在前
    end
end
fprintf("};\n");

%     K11(i)=K_s{1,i}(1,1);%先找到数组头地址再找其内容
%     K12(i)=K_s{1,i}(1,2);
%     K13(i)=K_s{1,i}(1,3);
%     K14(i)=K_s{1,i}(1,4);
%     K15(i)=K_s{1,i}(1,5);
%     K16(i)=K_s{1,i}(1,6);
%     K17(i)=K_s{1,i}(1,7);
%     K18(i)=K_s{1,i}(1,8);
%     K19(i)=K_s{1,i}(1,9);
%     K110(i)=K_s{1,i}(1,10);
% 
%     K21(i)=K_s{1,i}(2,1);%先找到数组头地址再找其内容
%     K22(i)=K_s{1,i}(2,2);
%     K23(i)=K_s{1,i}(2,3);
%     K24(i)=K_s{1,i}(2,4);
%     K25(i)=K_s{1,i}(2,5);
%     K26(i)=K_s{1,i}(2,6);
%     K27(i)=K_s{1,i}(2,7);
%     K28(i)=K_s{1,i}(2,8);
%     K29(i)=K_s{1,i}(2,9);
%     K210(i)=K_s{1,i}(2,10);
% 
%     K31(i)=K_s{1,i}(3,1);
%     K32(i)=K_s{1,i}(3,2);
%     K33(i)=K_s{1,i}(3,3);
%     K34(i)=K_s{1,i}(3,4);
%     K35(i)=K_s{1,i}(3,5);
%     K36(i)=K_s{1,i}(3,6);
%     K37(i)=K_s{1,i}(3,7);
%     K38(i)=K_s{1,i}(3,8);
%     K39(i)=K_s{1,i}(3,9);
%     K310(i)=K_s{1,i}(3,10);
% 
%     K41(i)=K_s{1,i}(4,1);
%     K42(i)=K_s{1,i}(4,2);
%     K43(i)=K_s{1,i}(4,3);
%     K44(i)=K_s{1,i}(4,4);
%     K45(i)=K_s{1,i}(4,5);
%     K46(i)=K_s{1,i}(4,6);
%     K47(i)=K_s{1,i}(4,7);
%     K48(i)=K_s{1,i}(4,8);
%     K49(i)=K_s{1,i}(4,9);
%     K410(i)=K_s{1,i}(4,10);
% end
% % 
% % disp('/*轮子距离环*/');
% % F11=fit_feedback(H_s,K11,3);
% % f111 = sprintf("leg_lqr_params[0][0].a = %d;\n",F11(1,1));
% % f112 = sprintf("leg_lqr_params[0][0].b = %d;\n",F11(1,2));
% % f113 = sprintf("leg_lqr_params[0][0].c = %d;\n",F11(1,3));
% % f114 = sprintf("leg_lqr_params[0][0].d = %d;\n\n",F11(1,4));
% % fprintf(f111);
% % fprintf(f112);
% % fprintf(f113);
% % fprintf(f114);
% 
% disp('/*轮子速度环*/');
% F12=fit_feedback(H_s,K12,3);
% f121 = sprintf(" leg_lqr_params[0][1].a = %d;\n",F12(1,1));
% f122 = sprintf(" leg_lqr_params[0][1].b = %d;\n",F12(1,2));
% f123 = sprintf(" leg_lqr_params[0][1].c = %d;\n",F12(1,3));
% f124 = sprintf(" leg_lqr_params[0][1].d = %d;\n\n",F12(1,4));
% fprintf(f121);
% fprintf(f122);
% fprintf(f123);
% fprintf(f124);
% 
% disp('/*轮子直立环*/');
% disp('//alpha');
% F13=fit_feedback(H_s,K13,3);
% f131 = sprintf(" leg_lqr_params[0][2].a = %d;\n",F13(1,1));
% f132 = sprintf(" leg_lqr_params[0][2].b = %d;\n",F13(1,2));
% f133 = sprintf(" leg_lqr_params[0][2].c = %d;\n",F13(1,3));
% f134 = sprintf(" leg_lqr_params[0][2].d = %d;\n\n",F13(1,4));
% fprintf(f131);
% fprintf(f132);
% fprintf(f133);
% fprintf(f134);
% 
% disp('//dalpha');
% F14=fit_feedback(H_s,K14,3);
% f141 = sprintf(" leg_lqr_params[0][3].a = %d;\n",F14(1,1));
% f142 = sprintf(" leg_lqr_params[0][3].b = %d;\n",F14(1,2));
% f143 = sprintf(" leg_lqr_params[0][3].c = %d;\n",F14(1,3));
% f144 = sprintf(" leg_lqr_params[0][3].d = %d;\n\n",F14(1,4));
% fprintf(f141);
% fprintf(f142);
% fprintf(f143);
% fprintf(f144);
% 
% disp('//beta');
% F15=fit_feedback(H_s,K15,3);
% f151 = sprintf(" leg_lqr_params[0][4].a = %d;\n",F15(1,1));
% f152 = sprintf(" leg_lqr_params[0][4].b = %d;\n",F15(1,2));
% f153 = sprintf(" leg_lqr_params[0][4].c = %d;\n",F15(1,3));
% f154 = sprintf(" leg_lqr_params[0][4].d = %d;\n\n",F15(1,4));
% fprintf(f151);
% fprintf(f152);
% fprintf(f153);
% fprintf(f154);
% 
% disp('//dbeta');
% F16=fit_feedback(H_s,K16,3);
% f161 = sprintf(" leg_lqr_params[0][5].a = %d;\n",F16(1,1));
% f162 = sprintf(" leg_lqr_params[0][5].b = %d;\n",F16(1,2));
% f163 = sprintf(" leg_lqr_params[0][5].c = %d;\n",F16(1,3));
% f164 = sprintf(" leg_lqr_params[0][5].d = %d;\n\n",F16(1,4));
% fprintf(f161);
% fprintf(f162);
% fprintf(f163);
% fprintf(f164);
% 
% % disp('/*轮腿路程环*/');
% % F31=fit_feedback(H_s,K31,3);
% % f311 = sprintf("leg_lqr_params[1][0].a = %d;\n",F31(1,1));
% % f312 = sprintf("leg_lqr_params[1][0].b = %d;\n",F31(1,2));
% % f313 = sprintf("leg_lqr_params[1][0].c = %d;\n",F31(1,3));
% % f314 = sprintf("leg_lqr_params[1][0].d = %d;\n\n",F31(1,4));
% % fprintf(f311);
% % fprintf(f312);
% % fprintf(f313);
% % fprintf(f314);
% 
% disp('/*轮腿速度环*/');
% F32=fit_feedback(H_s,K32,3);
% f321 = sprintf(" leg_lqr_params[1][1].a = %d;\n",F32(1,1));
% f322 = sprintf(" leg_lqr_params[1][1].b = %d;\n",F32(1,2));
% f323 = sprintf(" leg_lqr_params[1][1].c = %d;\n",F32(1,3));
% f324 = sprintf(" leg_lqr_params[1][1].d = %d;\n\n",F32(1,4));
% fprintf(f321);
% fprintf(f322);
% fprintf(f323);
% fprintf(f324);
% 
% disp('/*轮腿直立环*/');
% disp('//alpha');
% F33=fit_feedback(H_s,K33,3);
% f331 = sprintf(" leg_lqr_params[1][2].a = %d;\n",F33(1,1));
% f332 = sprintf(" leg_lqr_params[1][2].b = %d;\n",F33(1,2));
% f333 = sprintf(" leg_lqr_params[1][2].c = %d;\n",F33(1,3));
% f334 = sprintf(" leg_lqr_params[1][2].d = %d;\n\n",F33(1,4));
% fprintf(f331);
% fprintf(f332);
% fprintf(f333);
% fprintf(f334);
% 
% disp('//dalpha');
% F34=fit_feedback(H_s,K34,3);
% f341 = sprintf(" leg_lqr_params[1][3].a = %d;\n",F34(1,1));
% f342 = sprintf(" leg_lqr_params[1][3].b = %d;\n",F34(1,2));
% f343 = sprintf(" leg_lqr_params[1][3].c = %d;\n",F34(1,3));
% f344 = sprintf(" leg_lqr_params[1][3].d = %d;\n\n",F34(1,4));
% fprintf(f341);
% fprintf(f342);
% fprintf(f343);
% fprintf(f344);
% 
% disp('//beta');
% F35=fit_feedback(H_s,K35,3);
% f351 = sprintf(" leg_lqr_params[1][4].a = %d;\n",F35(1,1));
% f352 = sprintf(" leg_lqr_params[1][4].b = %d;\n",F35(1,2));
% f353 = sprintf(" leg_lqr_params[1][4].c = %d;\n",F35(1,3));
% f354 = sprintf(" leg_lqr_params[1][4].d = %d;\n\n",F35(1,4));
% fprintf(f351);
% fprintf(f352);
% fprintf(f353);
% fprintf(f354);
% 
% disp('//dbeta');
% F36=fit_feedback(H_s,K36,3);
% f361 = sprintf(" leg_lqr_params[1][5].a = %d;\n",F36(1,1));
% f362 = sprintf(" leg_lqr_params[1][5].b = %d;\n",F36(1,2));
% f363 = sprintf(" leg_lqr_params[1][5].c = %d;\n",F36(1,3));
% f364 = sprintf(" leg_lqr_params[1][5].d = %d;\n\n",F36(1,4));
% fprintf(f361);
% fprintf(f362);
% fprintf(f363);
% fprintf(f364);

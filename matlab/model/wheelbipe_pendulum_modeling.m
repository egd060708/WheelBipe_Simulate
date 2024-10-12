clc;
clear;

%轮足二级倒立摆理想化拉格朗日建模

syms t real;    %定义微分时间
syms Twl Twr Tjl Tjr;  %定义输入力矩,左轮力矩，右轮力矩，关节力矩
syms Mw R Jw D g real;   %定义轮子相关物理参数
syms Mll Mlr Hll Hlr Jll Jlr Pll Plr real;   %定义轮杆相关物理参数
syms Mb Pb Jbpitch Jbyaw real; %定义车体相关物理参数
syms alpha(t) beta(t) x(t); %定义广义坐标系，定义为时间的函数

%定义微分项
dx = diff(x,t);
dalpha = diff(alpha,t);
dbeta = diff(beta,t);
%定义质点速度
V1x = dx + L1*dalpha*cos(alpha);
V1y = -L1*dalpha*sin(alpha);
V2x = dx + H1*dalpha*cos(alpha) + L2*dbeta*cos(beta);
V2y = -H1*dalpha*sin(alpha) - L2*dbeta*sin(beta);
%定义动能
K1 = Jw*(dx^2)/(R^2); %轮子转动动能
K2 = Mw*(dx^2); %轮子平动动能
K3 = J1*(dalpha^2);
K4 = M1*(V1x^2 + V1y^2);
K5 = 0.5*J2*(dbeta^2);
K6 = 0.5*M2*(V2x^2 + V2y^2);
%定义势能
P1 = 2*M1*g*L1*cos(alpha);
P2 = M2*g*(H1*cos(alpha)+L2*cos(beta));

%拉格朗日函数，并使用符号替代一阶求导项
syms alpha_dot(t) beta_dot(t) x_dot(t);
syms alpha_ddot(t) beta_ddot(t) x_ddot(t);
%构建全局拉格朗日函数
L_global = K1+K2+K3+K4+K5+K6-P1-P2;
L_global = subs(L_global,diff(alpha,t),alpha_dot);
L_global = subs(L_global,diff(beta,t),beta_dot);
L_global = subs(L_global,diff(x,t),x_dot);
L_global = simplify(L_global);
% disp(L_global);
%q=x
L_xdot = functionalDerivative(L_global, x_dot);
L_xdot_t = diff(L_xdot, t);         % 二阶导
L_xdot_t = subs(L_xdot_t, diff(alpha,t), alpha_dot); % 把二阶导中的diff项替换为符号变量
L_xdot_t = subs(L_xdot_t, diff(beta,t), beta_dot);
L_xdot_t = subs(L_xdot_t, diff(x,t), x_dot);
L_xdot_t = subs(L_xdot_t, diff(alpha_dot,t), alpha_ddot);
L_xdot_t = subs(L_xdot_t, diff(beta_dot,t), beta_ddot);
L_xdot_t = subs(L_xdot_t, diff(x_dot,t), x_ddot);
L_xdot_t = simplify(L_xdot_t);
L_xdot_t = collect(L_xdot_t, [x_dot,alpha_dot,beta_dot,x_ddot,alpha_ddot,beta_ddot]);
L_x = functionalDerivative(L_global, x);
% disp(L_xdot);
% disp(L_xdot_t);
% disp(L_x);
%q=alpha
L_alphadot = functionalDerivative(L_global, alpha_dot);
L_alphadot_t = diff(L_alphadot, t);         % 二阶导
L_alphadot_t = subs(L_alphadot_t, diff(alpha,t), alpha_dot); % 把二阶导中的diff项替换为符号变量
L_alphadot_t = subs(L_alphadot_t, diff(beta,t), beta_dot);
L_alphadot_t = subs(L_alphadot_t, diff(x,t), x_dot);
L_alphadot_t = subs(L_alphadot_t, diff(alpha_dot,t), alpha_ddot);
L_alphadot_t = subs(L_alphadot_t, diff(beta_dot,t), beta_ddot);
L_alphadot_t = subs(L_alphadot_t, diff(x_dot,t), x_ddot);
L_alphadot_t = simplify(L_alphadot_t);
L_alphadot_t = collect(L_alphadot_t, [x_dot,alpha_dot,beta_dot,x_ddot,alpha_ddot,beta_ddot]);
L_alpha = functionalDerivative(L_global, alpha);
% disp(L_alphadot);
% disp(L_alphadot_t);
% disp(L_alpha);
%q=beta
L_betadot = functionalDerivative(L_global, beta_dot);
L_betadot_t = diff(L_betadot, t);         % 二阶导
L_betadot_t = subs(L_betadot_t, diff(alpha,t), alpha_dot); % 把二阶导中的diff项替换为符号变量
L_betadot_t = subs(L_betadot_t, diff(beta,t), beta_dot);
L_betadot_t = subs(L_betadot_t, diff(x,t), x_dot);
L_betadot_t = subs(L_betadot_t, diff(alpha_dot,t), alpha_ddot);
L_betadot_t = subs(L_betadot_t, diff(beta_dot,t), beta_ddot);
L_betadot_t = subs(L_betadot_t, diff(x_dot,t), x_ddot);
L_betadot_t = simplify(L_betadot_t);
L_betadot_t = collect(L_betadot_t, [x_dot,alpha_dot,beta_dot,x_ddot,alpha_ddot,beta_ddot]);
L_beta = functionalDerivative(L_global, beta);
% disp(L_betadot);
% disp(L_betadot_t);
% disp(L_beta);

%动力学方程左侧
left_x = L_xdot_t - L_x;    %拉格朗日方程
left_x = simplify(left_x);  %整理一遍
left_x = collect(left_x,[x_dot,alpha_dot,beta_dot,x_ddot,alpha_ddot,beta_ddot]);%整理变量表达式

left_alpha = L_alphadot_t - L_alpha;
left_alpha = simplify(left_alpha);
left_alpha = collect(left_alpha,[x_dot,alpha_dot,beta_dot,x_ddot,alpha_ddot,beta_ddot]);

left_beta = L_betadot_t - L_beta;
left_beta = simplify(left_beta);
left_beta = collect(left_beta,[x_dot,alpha_dot,beta_dot,x_ddot,alpha_ddot,beta_ddot]);

% disp(left_x);
% disp(left_alpha);
% disp(left_beta);
%线性化方程
left_x = subs(left_x,[cos(alpha(t)) cos(beta(t))],[1 1]);%线性化
left_x = subs(left_x,[sin(alpha(t)) sin(beta(t))],[alpha(t) beta(t)]);
left_x = subs(left_x,[alpha_dot(t) beta_dot(t)],[0 0]);
left_alpha = subs(left_alpha,[cos(alpha(t)) cos(alpha(t) - beta(t))],[1 1]);
left_alpha = subs(left_alpha,[sin(alpha(t)) sin(alpha(t) - beta(t))],[alpha(t) (alpha(t)-beta(t))]);
left_alpha = subs(left_alpha,beta_dot(t),0);
left_beta = subs(left_beta,[cos(beta(t)) cos(alpha(t) - beta(t))],[1 1]);
left_beta = subs(left_beta,[sin(beta(t)) sin(alpha(t) - beta(t))],[beta(t) (alpha(t) - beta(t))]);
left_beta = subs(left_beta,alpha_dot(t),0);
disp(left_x);
disp(left_alpha);
disp(left_beta);

%整理动力学方程成如下形式（此形式是根据微分的结果得出最后可能存在的广义坐标项决定的）
%A[q_ddot] + B[q] = J*T
% D矩阵
vars =[x(t),alpha(t),beta(t),x_ddot(t),alpha_ddot(t),beta_ddot(t)];
value =[0 0 0 0 0 0];

%下面通过替换的方法把各种变量一个个筛选出来
%B和A矩阵
syms n;
n = eye(6);
for i = 1:3
    value = n(i,:);%一行一行地按[x(t),alpha(t),beta(t),x_ddot(t),alpha_ddot(t),beta_ddot(t)]顺序把变量取出来
    B(1,i) = subs(left_x,vars,value);
    B(2,i) = subs(left_alpha,vars,value);
    B(3,i) = subs(left_beta,vars,value);
end

%A矩阵
for i = 1:3
    value = n(3+i,:);%一行一行地按[x(t),alpha(t),beta(t),x_ddot(t),alpha_ddot(t),beta_ddot(t)]顺序把变量取出来
    A(1,i) = subs(left_x,vars,value);
    A(2,i) = subs(left_alpha,vars,value);
    A(3,i) = subs(left_beta,vars,value);
end
% 广义力矩阵J，广义力为轮驱动力矩和上下层关节力矩
J = [1/R,1/R,0;...
    -1,-1,1;...
    0,0,-1];
disp(A);
disp(B);

% 获取当前脚本的完整路径
currentFile = mfilename('fullpath');

% 使用fileparts分离出路径部分
[currentPath,~,~] = fileparts(currentFile);
save(strcat(currentPath,'/Wheelbipe_Pendulum_Dynamics.mat'), 'A', 'B', 'J');
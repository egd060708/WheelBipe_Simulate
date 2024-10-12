clc;
clear;

%轮足二级倒立摆结合转向拉格朗日建模

syms t real; %定义微分时间
syms Twl Twr Tjl Tjr; %定义输入力矩，[左右轮，左右关节]
syms Mw R Jw D real; %定义驱动轮相关物理量，[轮质量，半径，惯量，轮距]
syms Ml Mr Hll Hlr Jll Jlr Pll Plr real; %定义腿杆相关物理量，[质量，长度，惯量，重心离轮轴位置]
syms Mb Pb Jbpitch Jbyaw real; %定义机体相关物理量，[质量，重心离关节位置，pitch惯量，yaw惯量]
syms thetawl(t) thetawr(t) thetall(t) thetalr(t) thetab(t); %定义广义坐标系，[位移，旋转角，左腿倾角，右腿倾角，机体倾角]
syms g real; %定义重力加速度

%定义微分项
dthetawl = diff(thetawl,t);
dthetawr = diff(thetawr,t);
dthetall = diff(thetall,t);
dthetalr = diff(thetalr,t);
dthetab = diff(thetab,t);


%运动学引入
xwl = R*thetawl;
xwr = R*thetawr;
zb = Hll/2*cos(thetall) + Hlr/2*cos(thetalr) + Pb*cos(thetab);
xll = xwl + Pll*sin(thetall);
xlr = xwr + Plr*sin(thetalr);
zll = Pll*cos(thetall);
zlr = Plr*cos(thetalr);
fi = (xwr-xwl)/D - Hll/D*sin(thetall) + Hlr/D*sin(thetalr);
xb = (xwl+xwr)/2 + Hll/2*sin(thetall) + Hlr/D*sin(thetalr) + Pb*sin(thetab);

%求运动学的一阶微分
dxwl = R*dthetawl;
dxwr = R*dthetawr;
dzb = -0.5*Hll*dthetall*sin(thetall) - 0.5*Hlr*dthetalr*sin(thetalr) - Pb*dthetab*sin(thetab);
dxll = dxwl + Pll*dthetall*cos(thetall);
dxlr = dxwr + Plr*dthetalr*cos(thetalr);
dzll = -Pll*dthetall*sin(thetall);
dzlr = -Plr*dthetalr*sin(thetalr);
dfi = (dxwr-dxwl)/D-Hll/D*dthetall*cos(thetall)+Hlr/D*dthetalr*cos(thetalr);
dxb = 0.5*(dxwl+dxwr)+0.5*Hll*dthetall*cos(thetall)+0.5*Hlr*dthetalr*cos(thetalr)+Pb*dthetab*cos(thetab);

%定义动能
Kwrl = 0.5*Jw*(dthetawl^2);
Kwrr = 0.5*Jw*(dthetawr^2);
Kwpl = 0.5*Mw*(dxwl^2);
Kwpr = 0.5*Mw*(dxwr^2);
Klrl = 0.5*Jll*(dthetall^2);
Klrr = 0.5*Jlr*(dthetalr^2);
Klpl = 0.5*Ml*(dxll^2+dzll^2);
Klpr = 0.5*Mr*(dxlr^2+dzlr^2);
Kbrpitch = 0.5*Jbpitch*(dthetab^2);
Kbryaw = 0.5*Jbyaw*(dfi^2);
Kbp = 0.5*Mb*(dxb^2+dzb^2);

%定义势能
Pell = Ml*g*zll;
Pelr = Mr*g*zlr;
Peb = Mb*g*zb;

%拉格朗日函数，并使用符号替代一阶求导项
syms thetawl_dot(t) thetawr_dot(t) thetall_dot(t) thetalr_dot(t) thetab_dot(t);
syms thetawl_ddot(t) thetawr_ddot(t) thetall_ddot(t) thetalr_ddot(t) thetab_ddot(t);

%状态变量微分
state_dot = [thetawl_dot,thetawr_dot,thetall_dot,thetalr_dot,thetab_dot,thetawl_ddot,thetawr_ddot,thetall_ddot,thetalr_ddot,thetab_ddot];


%构建全局拉格朗日函数
L_global = Kwrl+Kwrr+Kwpl+Kwpr+Klrl+Klrr+Klpl+Klpr+Kbrpitch+Kbryaw+Kbp - Pell - Pelr - Peb;
L_global = subs(L_global,diff(thetawl,t),thetawl_dot);
L_global = subs(L_global,diff(thetawr,t),thetawr_dot);
L_global = subs(L_global,diff(thetall,t),thetall_dot);
L_global = subs(L_global,diff(thetalr,t),thetalr_dot);
L_global = subs(L_global,diff(thetab,t),thetab_dot);
L_global = simplify(L_global);
pretty(L_global);

%x
L_thetawldot = functionalDerivative(L_global, thetawl_dot);
L_thetawldot_t = diff(L_thetawldot, t);

L_thetawldot_t = subs(L_thetawldot_t,diff(thetawl,t),thetawl_dot);
L_thetawldot_t = subs(L_thetawldot_t,diff(thetawr,t),thetawr_dot);
L_thetawldot_t = subs(L_thetawldot_t,diff(thetall,t),thetall_dot);
L_thetawldot_t = subs(L_thetawldot_t,diff(thetalr,t),thetalr_dot);
L_thetawldot_t = subs(L_thetawldot_t,diff(thetab,t),thetab_dot);

L_thetawldot_t = subs(L_thetawldot_t,diff(thetawl_dot,t),thetawl_ddot);
L_thetawldot_t = subs(L_thetawldot_t,diff(thetawr_dot,t),thetawr_ddot);
L_thetawldot_t = subs(L_thetawldot_t,diff(thetall_dot,t),thetall_ddot);
L_thetawldot_t = subs(L_thetawldot_t,diff(thetalr_dot,t),thetalr_ddot);
L_thetawldot_t = subs(L_thetawldot_t,diff(thetab_dot,t),thetab_ddot);

L_thetawldot_t = simplify(L_thetawldot_t);
L_thetawldot_t = collect(L_thetawldot_t,state_dot);
L_thetawl = functionalDerivative(L_global,thetawl);

%fi
L_thetawrdot = functionalDerivative(L_global, thetawr_dot);
L_thetawrdot_t = diff(L_thetawrdot, t);

L_thetawrdot_t = subs(L_thetawrdot_t,diff(thetawl,t),thetawl_dot);
L_thetawrdot_t = subs(L_thetawrdot_t,diff(thetawr,t),thetawr_dot);
L_thetawrdot_t = subs(L_thetawrdot_t,diff(thetall,t),thetall_dot);
L_thetawrdot_t = subs(L_thetawrdot_t,diff(thetalr,t),thetalr_dot);
L_thetawrdot_t = subs(L_thetawrdot_t,diff(thetab,t),thetab_dot);

L_thetawrdot_t = subs(L_thetawrdot_t,diff(thetawl_dot,t),thetawl_ddot);
L_thetawrdot_t = subs(L_thetawrdot_t,diff(thetawr_dot,t),thetawr_ddot);
L_thetawrdot_t = subs(L_thetawrdot_t,diff(thetall_dot,t),thetall_ddot);
L_thetawrdot_t = subs(L_thetawrdot_t,diff(thetalr_dot,t),thetalr_ddot);
L_thetawrdot_t = subs(L_thetawrdot_t,diff(thetab_dot,t),thetab_ddot);

L_thetawrdot_t = simplify(L_thetawrdot_t);
L_thetawrdot_t = collect(L_thetawrdot_t,state_dot);
L_thetawr = functionalDerivative(L_global,thetawr);

%thetall
L_thetalldot = functionalDerivative(L_global, thetall_dot);
L_thetalldot_t = diff(L_thetalldot, t);

L_thetalldot_t = subs(L_thetalldot_t,diff(thetawl,t),thetawl_dot);
L_thetalldot_t = subs(L_thetalldot_t,diff(thetawr,t),thetawr_dot);
L_thetalldot_t = subs(L_thetalldot_t,diff(thetall,t),thetall_dot);
L_thetalldot_t = subs(L_thetalldot_t,diff(thetalr,t),thetalr_dot);
L_thetalldot_t = subs(L_thetalldot_t,diff(thetab,t),thetab_dot);

L_thetalldot_t = subs(L_thetalldot_t,diff(thetawl_dot,t),thetawl_ddot);
L_thetalldot_t = subs(L_thetalldot_t,diff(thetawr_dot,t),thetawr_ddot);
L_thetalldot_t = subs(L_thetalldot_t,diff(thetall_dot,t),thetall_ddot);
L_thetalldot_t = subs(L_thetalldot_t,diff(thetalr_dot,t),thetalr_ddot);
L_thetalldot_t = subs(L_thetalldot_t,diff(thetab_dot,t),thetab_ddot);

L_thetalldot_t = simplify(L_thetalldot_t);
L_thetalldot_t = collect(L_thetalldot_t,state_dot);
L_thetall = functionalDerivative(L_global,thetall);

%thetalr
L_thetalrdot = functionalDerivative(L_global, thetalr_dot);
L_thetalrdot_t = diff(L_thetalrdot, t);

L_thetalrdot_t = subs(L_thetalrdot_t,diff(thetawl,t),thetawl_dot);
L_thetalrdot_t = subs(L_thetalrdot_t,diff(thetawr,t),thetawr_dot);
L_thetalrdot_t = subs(L_thetalrdot_t,diff(thetall,t),thetall_dot);
L_thetalrdot_t = subs(L_thetalrdot_t,diff(thetalr,t),thetalr_dot);
L_thetalrdot_t = subs(L_thetalrdot_t,diff(thetab,t),thetab_dot);

L_thetalrdot_t = subs(L_thetalrdot_t,diff(thetawl_dot,t),thetawl_ddot);
L_thetalrdot_t = subs(L_thetalrdot_t,diff(thetawr_dot,t),thetawr_ddot);
L_thetalrdot_t = subs(L_thetalrdot_t,diff(thetall_dot,t),thetall_ddot);
L_thetalrdot_t = subs(L_thetalrdot_t,diff(thetalr_dot,t),thetalr_ddot);
L_thetalrdot_t = subs(L_thetalrdot_t,diff(thetab_dot,t),thetab_ddot);

L_thetalrdot_t = simplify(L_thetalrdot_t);
L_thetalrdot_t = collect(L_thetalrdot_t,state_dot);
L_thetalr = functionalDerivative(L_global,thetalr);

%thetab
L_thetabdot = functionalDerivative(L_global, thetab_dot);
L_thetabdot_t = diff(L_thetabdot, t);

L_thetabdot_t = subs(L_thetabdot_t,diff(thetawl,t),thetawl_dot);
L_thetabdot_t = subs(L_thetabdot_t,diff(thetawr,t),thetawr_dot);
L_thetabdot_t = subs(L_thetabdot_t,diff(thetall,t),thetall_dot);
L_thetabdot_t = subs(L_thetabdot_t,diff(thetalr,t),thetalr_dot);
L_thetabdot_t = subs(L_thetabdot_t,diff(thetab,t),thetab_dot);

L_thetabdot_t = subs(L_thetabdot_t,diff(thetawl_dot,t),thetawl_ddot);
L_thetabdot_t = subs(L_thetabdot_t,diff(thetawr_dot,t),thetawr_ddot);
L_thetabdot_t = subs(L_thetabdot_t,diff(thetall_dot,t),thetall_ddot);
L_thetabdot_t = subs(L_thetabdot_t,diff(thetalr_dot,t),thetalr_ddot);
L_thetabdot_t = subs(L_thetabdot_t,diff(thetab_dot,t),thetab_ddot);

L_thetabdot_t = simplify(L_thetabdot_t);
L_thetabdot_t = collect(L_thetabdot_t,state_dot);
L_thetab = functionalDerivative(L_global,thetab);

%动力学方程左侧
left_thetawl = L_thetawldot_t - L_thetawl;
left_thetawl = simplify(left_thetawl);
left_thetawl = collect(left_thetawl,state_dot);

left_thetawr = L_thetawrdot_t - L_thetawr;
left_thetawr = simplify(left_thetawr);
left_thetawr = collect(left_thetawr,state_dot);

left_thetall = L_thetalldot_t - L_thetall;
left_thetall = simplify(left_thetall);
left_thetall = collect(left_thetall,state_dot);

left_thetalr = L_thetalrdot_t - L_thetalr;
left_thetalr = simplify(left_thetalr);
left_thetalr = collect(left_thetalr,state_dot);

left_thetab = L_thetabdot_t - L_thetab;
left_thetab = simplify(left_thetab);
left_thetab = collect(left_thetab,state_dot);

%线性化方程
left_thetawl = subs(left_thetawl,[cos(thetall(t)) cos(thetalr(t)) cos(thetab(t))],[1 1 1]);
left_thetawl = subs(left_thetawl,[sin(thetall(t)) sin(thetalr(t)) sin(thetab(t))],[thetall(t) thetalr(t) thetab(t)]);
left_thetawl = subs(left_thetawl,[thetall_dot(t) thetalr_dot(t) thetab_dot(t)],[0 0 0]);

left_thetawr = subs(left_thetawr,[cos(thetall(t)) cos(thetalr(t)) cos(thetab(t))],[1 1 1]);
left_thetawr = subs(left_thetawr,[sin(thetall(t)) sin(thetalr(t)) sin(thetab(t))],[thetall(t) thetalr(t) thetab(t)]);
left_thetawr = subs(left_thetawr,[thetall_dot(t) thetalr_dot(t) thetab_dot(t)],[0 0 0]);

left_thetall = subs(left_thetall,[cos(thetall(t)) cos(thetalr(t)) cos(thetab(t))],[1 1 1]);
left_thetall = subs(left_thetall,[sin(thetall(t)) sin(thetalr(t)) sin(thetab(t))],[thetall(t) thetalr(t) thetab(t)]);
left_thetall = subs(left_thetall,[thetall_dot(t) thetalr_dot(t) thetab_dot(t)],[0 0 0]);

left_thetalr = subs(left_thetalr,[cos(thetall(t)) cos(thetalr(t)) cos(thetab(t))],[1 1 1]);
left_thetalr = subs(left_thetalr,[sin(thetall(t)) sin(thetalr(t)) sin(thetab(t))],[thetall(t) thetalr(t) thetab(t)]);
left_thetalr = subs(left_thetalr,[thetall_dot(t) thetalr_dot(t) thetab_dot(t)],[0 0 0]);

left_thetab = subs(left_thetab,[cos(thetall(t)) cos(thetalr(t)) cos(thetab(t))],[1 1 1]);
left_thetab = subs(left_thetab,[sin(thetall(t)) sin(thetalr(t)) sin(thetab(t))],[thetall(t) thetalr(t) thetab(t)]);
left_thetab = subs(left_thetab,[thetall_dot(t) thetalr_dot(t) thetab_dot(t)],[0 0 0]);

% disp(L_global);
% disp(left_thetawl);
% disp(left_thetawr);
% disp(left_thetall);
% disp(left_thetalr);
% disp(left_thetab);

% 整理动力学方程成基本形式
% M[q_ddot] + V[q_dot,q] + G[q] = J*T
% 为了线性化，去除离心力和哥氏力矩阵
vars = [thetawl(t) thetawr(t) thetall(t) thetalr(t) thetab(t) thetawl_ddot(t) thetawr_ddot(t) thetall_ddot(t) thetalr_ddot(t) thetab_ddot(t)];
values = [0 0 0 0 0 0 0 0 0 0];

syms n;
n = eye(10);
% G = zeros(5);
% M = zeros(5);
for i = 1:5
    values = n(i,:);
    G(1,i) = subs(left_thetawl,vars,values);
    G(2,i) = subs(left_thetawr,vars,values);
    G(3,i) = subs(left_thetall,vars,values);
    G(4,i) = subs(left_thetalr,vars,values);
    G(5,i) = subs(left_thetab,vars,values);
end

for i = 1:5
    values = n(5+i,:);
    M(1,i) = subs(left_thetawl,vars,values);
    M(2,i) = subs(left_thetawr,vars,values);
    M(3,i) = subs(left_thetall,vars,values);
    M(4,i) = subs(left_thetalr,vars,values);
    M(5,i) = subs(left_thetab,vars,values);
end

J = [1,0,0,0;...
     0,1,0,0;...
     -1,0,1,0;...
     0,-1,0,1;...
     -1,-1,-1,-1];

% 保存动力学方程
% 获取当前脚本的完整路径
currentFile = mfilename('fullpath');

% 使用fileparts分离出路径部分
[currentPath,~,~] = fileparts(currentFile);
save(strcat(currentPath,'/BalanceTurn_Model.mat'), 'M', 'G', 'J');
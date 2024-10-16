% 用于webots仿真的轮足机器人物理参数

% 轮参数
Mw_= 0.308;%kg
R_= 0.07;%m
Jw_= (754.296e-6);%kg/m^2
D_= 0.54;%m

% 腿部参数
Llinks = [0.15,0.15,0.288];%无论是并联腿还是串联腿都可以用三个参数表达
Mlinks = [0.189,0.327];
Ilinks = [608.165e-6,2997.229e-6];
k_ = (4*Mlinks(1)+2*Mlinks(2))/(8*(Mlinks(1)+Mlinks(2)));

% Lleg = 0.3;
% [Pleg,Mleg,Ileg] = leg_param_cal(Lleg,Llinks,Ilinks,Mlinks);
% Ml_= Mleg;%kg
% Mr_= Mleg;
% Hll_= Lleg;
% Hlr_= Lleg;
% Jll_= Ileg;
% Jlr_= Ileg;
% Pll_= Pleg;
% Plr_= Pleg;

% 机体参数
Mb_= 19.919;%kg
Pb_= 0.0315;%m
Jbpitch_ = (424025.762e-6);%kg/m^2
Jbyaw_ = (438071.894e-6);

% 重力参数
g_ = 9.81;


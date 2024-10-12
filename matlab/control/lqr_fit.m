%lqr线性化函数(返回参数：lqr参数，对应轮杆高度，参数组数）
function [K_,H_,C_]=lqr_fit(Ts,top,bottom,step,lqr_Q,lqr_R,model_name)
count=0;
for n=bottom:step:top
    count=count+1;
end

%建立单元数组，开辟空间
K_=cell(1,count);
H_=zeros(1,count);
C_=count;

xNum = size(lqr_Q,1);%获取行数，也就是状态变量个数
uNum = size(lqr_R,2);%获取列数，也就是输入变量个数

%导入模型
load(model_name);
%syms t real; %定义微分时间
%syms Twl Twr Tjl Tjr; %定义输入力矩，[左右轮，左右关节]
syms Mw R Jw D real; %定义驱动轮相关物理量，[轮质量，半径，惯量，轮距]
syms Ml Mr Hll Hlr Jll Jlr Pll Plr real; %定义腿杆相关物理量，[质量，长度，惯量，重心离轮轴位置]
syms Mb Pb Jbpitch Jbyaw real; %定义机体相关物理量，[质量，重心离关节位置，pitch惯量，yaw惯量]
%syms thetawl(t) thetawr(t) thetall(t) thetalr(t) thetab(t); %定义广义坐标系，[位移，旋转角，左腿倾角，右腿倾角，机体倾角]
syms g real; %定义重力加速度

% 导入物理参数
webots_param;

m=1;
for H1_=bottom:step:top
    % 求解部分物理参数
    [Pleg,Mleg,Ileg] = leg_param_cal(H1_,Llinks,Ilinks,Mlinks);
    Ml_= Mleg;%kg
    Mr_= Mleg;
    Hll_= H1_;
    Hlr_= H1_;
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
    lqr_A = zeros(xNum,xNum);
    lqr_B = zeros(xNum,uNum);
    for i=1:1:xNum
        k = rem(i,2);
        if k == 1
            lqr_A(i,i+1) = 1;
        else
            for j=1:1:xNum
                l = rem(j,2);
                if l == 1
                    lqr_A(i,j) = H(i/2,(j+1)/2);
                end
            end
            for n=1:1:uNum
                lqr_B(i,n) = I(i/2,n);
            end
        end
    end
    lqr_C = eye(xNum);
    lqr_D = zeros(xNum,uNum);
    sys_c = ss(lqr_A,lqr_B,lqr_C,lqr_D);
    sys_d = c2d(sys_c, Ts);
    %判断可控性
if (rank(ctrb(sys_d.A,sys_d.B))==xNum)    
    temp=dlqr(sys_d.A,sys_d.B,lqr_Q,lqr_R);
    cell_temp = cell(1,1);
    for i = 1:uNum
        for j = 1:xNum
            cell_temp{1}(i, j) = temp(i, j);
        end
    end
    K_(1,m)=cell_temp;
    H_(1,m)=H1_;
else
    K_(1,m)=0;
    H_(1,m)=H1_;
    disp('Uncontrollable!');
end
    m = m+1;
end
end
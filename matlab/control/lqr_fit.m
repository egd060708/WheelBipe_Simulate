%lqr线性化函数(返回参数：lqr参数，对应轮杆高度，参数组数）
function [K_,P_,H_,C_]=lqr_fit(Ts,top,bottom,step,lqr_Q,lqr_R,model_name,robot_type)
count=0;
for n=bottom:step:top
    count=count+1;
end

%建立单元数组，开辟空间
K_=cell(1,count);
P_=cell(1,count);
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
    if robot_type == 0
        Ml_ = 2*(Mlinks(1)+Mlinks(2));
        Mr_ = 2*(Mlinks(1)+Mlinks(2));
        Hll_ = H1_;
        Hlr_ = H1_;
        Pll_ = H1_*(1-k_);
        Plr_ = H1_*(1-k_);
        Jll_ = (Mlinks(2)*2*((Llinks(3)/2).^2) + Mlinks(1)*(2*(Llinks(1)/2).^2 + Llinks(2).^2 / 2) + Mlinks(2)*(2*(Llinks(1)/2).^2 + 2*Llinks(1).^2))- Ml_ * (k_*Hll_).^2;
        Jlr_ = (Mlinks(2)*2*((Llinks(3)/2).^2) + Mlinks(1)*(2*(Llinks(1)/2).^2 + Llinks(2).^2 / 2) + Mlinks(2)*(2*(Llinks(1)/2).^2 + 2*Llinks(1).^2))- Mr_ * (k_*Hlr_).^2; 
    else if robot_type == 1
        Ml_ = Mlinks(1)+Mlinks(2)+Mlinks(3)+Mlinks(4);
        Mr_ = Ml_;
        Hll_ = H1_;
        Hlr_ = H1_;
        [pim,jm] = IM_for_serialLeg_fact(Llinks,Mlinks,Ilinks,H1_,0);
        Pll_ = H1_+pim(2);
        Plr_ = Pll_;
        Jll_ = jm;
        Jlr_ = Jll_;
    end
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
    [tempK,tempP,cls]=dlqr(sys_d.A,sys_d.B,lqr_Q,lqr_R);
    cell_tempK = cell(1,1);
    for i = 1:uNum
        for j = 1:xNum
            cell_tempK{1}(i, j) = tempK(i, j);
        end
    end
    K_(1,m)=cell_tempK;
    cell_tempP = cell(1,1);
    for i = 1:xNum
        for j = 1:xNum
            cell_tempP{1}(i,j) = tempP(i, j);
        end
    end
    P_(1,m)=cell_tempP;
    H_(1,m)=H1_;
else
    K_(1,m)=cell(1,1);
    P_(1,m)=cell(1,1);
    H_(1,m)=H1_;
    disp('Uncontrollable!');
end
    m = m+1;
end
end
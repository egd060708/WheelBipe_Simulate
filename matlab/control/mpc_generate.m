%%
%拟合mpc模型参数导出
clc;
clear;

x_num = 10;
u_num = 4;

%设定轮杆长度
h_top=0.4;
h_bottom=0.13;
h_step=0.01;
% Ts = 0.005;%qpOASES
Ts = 0.005;%tinyMPC
% Ts = 0.005;%qp++
robot_type = 1;

%取出拟合mpc模型参数
[K_A,K_B,H_s,C_s,cA,cB,cH] = mpc_model_fit (Ts, h_top,h_bottom,h_step,x_num,u_num,'D:/Git_Project/github/WheelBipe_Simulate/matlab/model/BalanceTurn_Model_Simple.mat',robot_type);

for j=1:1:x_num
    for k=1:1:x_num
        for i=C_s:-1:1
            A(j,k,i) = K_A{1,i}(j,k);%取出每一个矩阵元素的组合
        end
    end
end


for j=1:1:x_num
    for k=1:1:u_num
        for i=C_s:-1:1
            B(j,k,i) = K_B{1,i}(j,k);
        end
    end
end

fprintf("double model_A[4*10*10]={\n");
for j=1:1:x_num
    for k=1:1:x_num
        for i = C_s:-1:1
            Y(i) = A(j,k,i);
        end
        % figure;
        % scatter(H_s, Y, 40, 'b', 'filled');  % 'filled'表示填充点
        % title('参数散点图');
        % xlabel('h');
        % ylabel('p');
        % grid on;
        % hold on;
        % coefficients = polyfit(H_s, Y, 3);
        % fittedY = polyval(coefficients, H_s);
        % plot(H_s, fittedY, 'r-', 'LineWidth', 2);

        F = fit_feedback(H_s,Y,3);
        fprintf(sprintf("%d,%d,%d,%d,\n",F(4),F(3),F(2),F(1)));%打印低次项在前
        % disp(cA);
        % disp(cB);
        % disp(cH);
    end
end
fprintf("};\n");

fprintf("double model_B[4*10*4]={\n");
for j=1:1:x_num
    for k=1:1:u_num
        for i = C_s:-1:1
            Y(i) = B(j,k,i);
        end
        % figure;
        % scatter(H_s, Y, 40, 'b', 'filled');  % 'filled'表示填充点
        % title('参数散点图');
        % xlabel('h');
        % ylabel('p');
        % grid on;
        % hold on;
        % coefficients = polyfit(H_s, Y, 3);
        % fittedY = polyval(coefficients, H_s);
        % plot(H_s, fittedY, 'r-', 'LineWidth', 2);

        F = fit_feedback(H_s,Y,3);
        fprintf(sprintf("%d,%d,%d,%d,\n",F(4),F(3),F(2),F(1)));%打印低次项在前
    end
end
fprintf("};\n");
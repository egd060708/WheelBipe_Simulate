function [p] = fit_feedback(x,y,n)
%n为拟合阶数
format long;             %把小数有效位拓展到15位
p= polyfit(x, y, n);    %进行n阶线性拟合，返回p为多项式系数向量
xi=x(1,1):0.01:x(1,end); %构建x轴向量
yi= polyval(p, xi);     %通过x轴向量计算y轴的值
yj=polyval(p,x);
si=y-yj;
si=std(si,0);
% plot(xi,yi,x,y,'r*');   %把曲线打印出来
% for i = n:-1:1          %从高次项开始每次减1，到1为止
%     c=sprintf("第%d次系数%d \n",i,p(1,n-i+1));
%     fprintf(c);
% end
% c=sprintf("常数项系数%d \n",p(1,n+1));
% fprintf(c);
% c=sprintf("样品标准差%d\n",si(1,1));
% fprintf(c);
end
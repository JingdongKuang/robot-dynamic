% [x, fval, exitflag, output] = fmincon(@fun, x0, A, b, Aeq, beq, lb, ub, nonlcon,options)


% fval 目标函数在最优解处的取值。
% exitflag 求解器的退出标志，表示是否收敛、达到最大迭代次数等等。
% output 包含求解器的输出信息，如迭代次数、函数值等等。
% lambda 包含各种约束的拉格朗日乘子
clc
clear all


% fun 匿名函数或句柄函数
f=@(x)0.6224*x(1)*x(2)*x(3)*x(4)+1.7781*x(2)*x(3)^2+3.1661*x(1)^2*x(4)+19.84*x(1)^2*x(3)


% A，b 线性不等式约束 Ax<=b
A=[]
b=[]
% Aeq，beq 线性不等式约束 Aeq x<=beq
Aeq=[]
beq=[]
% lb 取值下限
lb=[0.0625,0.0625,10,10]
% ub 取值上限
ub=[6.1875,6.1875,200,200]
% x0 猜测初始值
x0=lb+1;



[x,f]=fmincon(f,x0,A,b,Aeq,beq,lb,ub,@constraint)


% nonlcon 非线性约束 c(x)<=0 ceq(x)=0
function [c,ceq] = constraint(x)
c=[0.0193*x(3)-x(1);
    0.00954*x(3)-x(2);
    750*1728-pi*x(3)^2*x(4)-4*pi*x(3)^3/3;
    x(4)-240];
ceq=[];
end



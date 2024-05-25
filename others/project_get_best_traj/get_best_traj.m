%用于计算极大线性无关组
clc;
clear;

func=@(x)getCondOfY(x);
nonlcon=@(x)constraint(x);
% There are no linear constraints, so set those arguments to |[]|. 

x0=load("x.txt"); 
A = [];
b = [];
Aeq = [];
beq = [];  
lb=[];
ub=[];
% Solve the problem. 
options = optimoptions('fmincon','Display','notify');
x = fmincon(func,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

writematrix(x,'x.txt');
coefficient_aa=reshape(x(1:30),5,6)';
coefficient_bb=reshape(x(31:end),5,6)';
fprintf('cond of x0:')
getCondOfY(x0)
fprintf('cond of x:')
getCondOfY(x)




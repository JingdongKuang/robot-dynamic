
%%利用标准D-H法建立多轴机器人
clear;
clc;


L1 = Link('d', 0.22,  'a', 0,     'alpha', pi/2, 'offset' ,0);    %Link 类函数;offset建立初始的偏转角
L2 = Link('d', 0,     'a', 0.455, 'alpha', 0,   'offset' ,pi/2);
L3 = Link('d', 0,     'a', 0,     'alpha', pi/2, 'offset' ,0);%offset pi/2
L4 = Link('d', 0.495, 'a', 0,     'alpha', pi/2,'offset', 0);
L5 = Link('d', 0,     'a', 0,     'alpha', -pi/2, 'offset' ,0);%offset 0
L6 = Link('d', 0.1565,'a', 0,     'alpha', 0,    'offset', 0);
L1.qlim = [-pi/2,pi];%利用qlim设置每个关节的旋转角度范围
robot=SerialLink([L1,L2,L3,L4,L5,L6],'name','E05L');   %SerialLink 类函数

%% 普通机器人的示教展示
%J = robot.jacob0([0 0 0 0 0 0]);
%T = robot.fkine([0 0 0 0 0 0]);
theta=deg2rad(15*ones(6,1));
J = robot.jacob0(theta);
robot.display();%展示出机器人的信息
teach(robot);%调出示教滑块
%disp(T)

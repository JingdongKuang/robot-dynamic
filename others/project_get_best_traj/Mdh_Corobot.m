%%利用标准D-H法建立多轴机器人
clear;
clc;
%质心位置
r1=[0.00002,0.02043,-0.051];
r2=[455-260.92,0.01,-112.48]*1e-3;
r3=[0.08,-41.21,-18.1]*1e-3;
r4=[0.03,-37.14,495-218.87]*1e-3;
r5=[0.02, -36.59,15.07]*1e-3;
r6=[-0.00003, -0.00106,-0.03404];

R1=[1,0,0;
    0,0,-1;
    0,1,0];
R3=[1,0,0;
    0,0,-1;
    0,1,0];
R4=[1,0,0;
    0,0,1;
    0,-1,0];
R5=[-1,0,0;
    0,0,-1;
    0,-1,0];
R6=[0,1,0;
    -1,0,0;
    0,0,1];
I1C=R1*[46167039.42,2116.69,26710.49;
    2116.69,21944229.27,8797091.54;
    26710.49,8797091.54, 39652637.03]*1e-9*R1';
I2C=[6131260.18,-4957.82,-2392349.81;
     -4957.82,134685326.69,963.85;
     -2392349.81,963.85,135087058.18]*1e-9;
I3C=R3*[21470885.80,10308.98,1091.43;
    10308.98,17583990.80,-4054753.75;
    1091.43,-4054753.75,10648104.35]*1e-9*R3';
I4C=R4*[49561362.63,13462.22,6498.44;
    13462.22,7464396.28,-12092659.64;
    6498.44,-12092659.64,45441302.73]*1e-9*R4';
I5C=R5*[7956816.25,1803.46,-3094.77;
    1803.46,6275966.23,-1492271.42;
    -3094.77,-1492271.42,3672461.05]*1e-9*R5';
I6C=R6*[483409.42,-276.90,-7019.04;
        -276.90,537684.75,345.58;
        -7019.04,345.58,585171.2]*1e-9*R6';

m1 = 8.30269;
m2 = 3.87623;
m3 = 5.31663;
m4 = 2.58124;
m5 = 2.66205;
m6 = 0.60427;
%m是质量，r是重心位置，I是
L1 = Link('d', 0.22,    'a', 0,       'alpha', 0,      'offset' ,0 ,'modified','m', m1, 'r', r1, 'I', I1C);    %Link 类函数;offset建立初始的偏转角
L2 = Link('d', 0,       'a', 0,       'alpha', pi/2,   'offset' ,pi/2, 'modified','m', m2, 'r', r2, 'I', I2C);
L3 = Link('d', 0,       'a', 0.455,   'alpha',0,       'offset' ,0, 'modified','m', m3, 'r', r3, 'I', I3C);%offset pi/2
L4 = Link('d', 0.495,   'a', 0,       'alpha', pi/2,   'offset', 0, 'modified','m', m4, 'r', r4, 'I', I4C);
L5 = Link('d', 0,       'a', 0,       'alpha', pi/2,   'offset' ,pi, 'modified','m', m5, 'r', r5, 'I', I5C);%offset 0
L6 = Link('d', 0.1565,  'a', 0,       'alpha', pi/2,  'offset', pi, 'modified','m', m6, 'r', r6, 'I', I6C);
%L1.qlim = [-pi/2,pi];%利用qlim设置每个关节的旋转角度范围
robot=SerialLink([L1,L2,L3,L4,L5,L6],'name','E05L');   %SerialLink 类函数


%robot.display();%展示出机器人的信息
%teach(robot);%调出示教滑块

toRad=pi/180;
q=[3, 2, 2, 4, 5, 6];
qdot=[3, 2, 4, 1, 7, 1];
qdotdot=[3, 1, 5, 6, 7, 1];
%robot.inertia(q);
fprintf('ture torque:')
true_torque=robot.rne_mdh(q,qdot,qdotdot)
torque=[];
x=load("x.txt");
coefficient_a=reshape(x(1:30),5,6)';
coefficient_b=reshape(x(31:end),5,6)';
data_size=2500;
for i=1:data_size
    q_qdot_qdotdot = getFourierTrajectory(coefficient_a,coefficient_b,0.004*i);
    torque_=robot.rne_mdh(q_qdot_qdotdot(:,1)',q_qdot_qdotdot(:,2)',q_qdot_qdotdot(:,3)');
    torque=[torque;torque_'];
end
ident_parameters=get_ident_parameters(torque,data_size);%load("ident_parameters.txt");%
%writematrix(ident_parameters,'ident_parameters')
ture_least_parameters=load('ture_least_parameters.txt');
ture_parameters=load('ture_parameters.txt');
tlp=ture_least_parameters(1,5);
tlp_long=[];
for i=1:6
    tlp_long=[tlp_long;ture_parameters(i,:)'];
end

for i=2:6
    tlp=[tlp;ture_least_parameters(i,:)'];
end
Ytilde = getYtilde(q, qdot, qdotdot);
Y=getY(q, qdot, qdotdot);
% fprintf('long torque:')
% long_torque=Y*tlp_long;
% transpose(long_torque)
fprintf('identified torque:')
identified_torque=Ytilde*ident_parameters;
transpose(identified_torque)
%robot.dyn

% robot2=SerialLink([L1,L2],'name','E05L2');
% q2=[0,0];
% T=robot2.fkine(q2);

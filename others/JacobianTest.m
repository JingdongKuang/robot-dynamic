clc;
clear;
format short;
syms q1 q2 q3 q4 q5 q6 para0 para1 para2 para3 para4 para5;
theta=[90,0,0,90,90,0];
toRad=pi/180;
%para=[para0,para1,para2,para3,para4,para5];
para=[0.22,0.455,0.495,0.1565];
% Define joint variables
%q = [q1; q2; q3; q4; q5; q6];
q=[0,0,0,0,0,0]*toRad;
%         a     alpha    d        theta 
DH_Table=[0,    pi/2,    0.22,    theta(1);
          0.455,0   ,    0   ,    theta(2)+pi/2;                       
          0    ,pi/2,    0   ,    theta(3);
          0    ,pi/2,    0.495,   theta(4);
          0    ,-pi/2,   0   ,    theta(5);
          0    ,0   ,    0.1565,  theta(6)];
% Define forward kinematics equations (example: a simple robot with only rotations)
%                  alpha,a,     d,      theta
T01 = dh_transform(pi/2 ,0    , para(1) , q(1));
T12 = dh_transform(0    ,para(2), 0     , q(2)+pi/2);
T23 = dh_transform(pi/2 ,0    , 0       , q(3));
T34 = dh_transform(pi/2 ,0    , para(3) , q(4));
T45 = dh_transform(-pi/2,0    , 0       , q(5));
T56 = dh_transform(0    ,0    , para(4) , q(6));

T02 = T01 * T12;
T03 = T02 * T23;
T04 = T03 * T34;
T05 = T04 * T45;
T06 = T05 * T56;
T67 = [1 0 0 0;];

% Extract position and orientation from the end-effector transformation matrix
p1_ = T01*transl(1,1,1);  % 重心的位置
p2_ = T02*transl(1,1,1);
p3_ = T03*transl(1,1,1);
p4_ = T04*transl(1,1,1);
p5_ = T05*transl(1,1,1);
p6_ = T06*transl(1,1,1);
p1 = p1_(1:3, 4);  % 重心的位置
p2 = p2_(1:3, 4);
p3 = p3_(1:3, 4);
p4 = p4_(1:3, 4);
p5 = p5_(1:3, 4);
p6 = p6_(1:3, 4);



% Linear velocity Jacobian (Jv)
Jv6q1 =diff(p6,q1);
Jv6q2 =diff(p6,q2);
Jv6q3 =diff(p6,q3);
Jv6q4 =diff(p6,q4);
Jv6q5 =diff(p6,q5);
Jv6q6 =diff(p6,q6);

%Jv6 = jacobian(p4,q);
%disp(Jv6q1)
%disp(vpa(Jv6))


% Define a function for DH transformation matrix
function T = dh_transform(alpha, a, d, theta)
    T = [cos(theta), -sin(theta)*cos(alpha), sin(theta)*sin(alpha), a*cos(theta);
         sin(theta), cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta);
         0, sin(alpha), cos(alpha), d;
         0, 0, 0, 1];
end




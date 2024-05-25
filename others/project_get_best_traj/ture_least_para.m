clear;
clc;
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


m=zeros(1,6);
I=cell(1,6);
rm=cell(1,6);
m(1) = 8.30269;
r1=[0.00002;0.02043;-0.051];
rm{1}=m(1)*r1;
I{1}=I1C+m(1)*(r1'*r1*eye(3)-r1*r1');

m(2)=3.87623;
r2=[455-260.92;
    0.01;
    -112.48]*1e-3;
rm{2}=m(2)*r2;
I{2}=I2C+m(2)*(r2'*r2*eye(3)-r2*r2');

m(3)=5.31663;
r3=[0.08;
    -41.21;
    -18.1]*1e-3;
rm{3}=m(3)*r3;
I{3}=I3C+m(3)*(r3'*r3*eye(3)-r3*r3');

m(4)=2.58124;
r4=[0.03;
    -37.14;
    495-218.87]*1e-3;
rm{4}=m(4)*r4;
I{4}=I4C+m(4)*(r4'*r4*eye(3)-r4*r4');

m(5)=2.66205;
r5=[0.02;
    -36.59;
    15.07]*1e-3;
rm{5}=m(5)*r5;
I{5}=I5C+m(5)*(r5'*r5*eye(3)-r5*r5');

m(6)=0.60427;
r6=[-0.00003;
    -0.00106;
    -0.03404];
rm{6}=m(6)*r6;
I{6}=I6C+m(6)*(r6'*r6*eye(3)-r6*r6');
%           theta   d       a       alpha
MDH_Table=[ 0,		0.22,	0,		0;
		    pi / 2, 0,		0,		pi/2;
		    0,		0,		0.455,	0;
		    0,		0.495,	0,		pi / 2;
		    pi,		0,		0,		pi / 2;
		    pi,		0.1565, 0,		pi/2];

ture_parameters=zeros(6,10);
for i=1:6
    ture_parameters(i,:)=[I{i}(1,1),I{i}(1,2),I{i}(1,3),I{i}(2,2),I{i}(2,3),I{i}(3,3),rm{i}',m(i)];
end
writematrix(ture_parameters,'ture_parameters');
ture_least_parameters=zeros(6,7);
parametes_recombin=zeros(6,10);
parametes_recombin(6,:)=ture_parameters(6,:);
for i=5:-1:1
    temp=parametes_recombin(i+1,:);
    MZ=temp(9);
    M=temp(10);
    YY=temp(4);
    d=MDH_Table(i+1,2);
    a=MDH_Table(i+1,3);
    alpha=MDH_Table(i+1,4);
    parametes_recombin(i,1)=ture_parameters(i,1) + temp(4) + 2*d*MZ + d^2*M;
    parametes_recombin(i,2)=ture_parameters(i,2) + a*sin(alpha)*(MZ + d*M);
    parametes_recombin(i,3)=ture_parameters(i,3) - a*cos(alpha)*(MZ + d*M);
    parametes_recombin(i,4)=ture_parameters(i,4) + cos(alpha)^2*YY + 2*d*cos(alpha)^2*MZ + (d^2*cos(alpha)^2 + a^2)*M;
    parametes_recombin(i,5)=ture_parameters(i,5) + cos(alpha)*(YY + 2*d*MZ + d^2*M);
    parametes_recombin(i,6)=ture_parameters(i,6) + sin(alpha)^2*YY + 2*d*sin(alpha)^2*MZ + (d^2*sin(alpha)^2 + a^2)*M;
    parametes_recombin(i,7)=ture_parameters(i,7) + a*M;
    parametes_recombin(i,8)=ture_parameters(i,8) - sin(alpha)*(MZ + d*M);
    parametes_recombin(i,9)=ture_parameters(i,9) + cos(alpha)*(MZ + d*M);
    parametes_recombin(i,10)=ture_parameters(i,10) + M;
end
ture_least_parameters(:,1)=parametes_recombin(:,1)-parametes_recombin(:,4);
ture_least_parameters(:,2)=parametes_recombin(:,2);
ture_least_parameters(:,3)=parametes_recombin(:,3);
ture_least_parameters(:,4)=parametes_recombin(:,5);
ture_least_parameters(:,5)=parametes_recombin(:,6);
ture_least_parameters(:,6)=parametes_recombin(:,7);
ture_least_parameters(:,7)=parametes_recombin(:,8);
writematrix(ture_least_parameters,'ture_least_parameters.txt');
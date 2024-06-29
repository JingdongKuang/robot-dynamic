% nonlcon 非线性约束 c(x)<=0 ceq(x)=0
function [c,ceq] = constraint(x)
coefficient_a=reshape(x(1:30),5,6)';
coefficient_b=reshape(x(31:end),5,6)';
qmax=[pi,pi*3/8,pi/2,pi,pi/2,pi];
qdotmax=[pi,pi,pi,pi,pi,pi];
qdotdotmax=[3*pi,3*pi,3*pi,3*pi,3*pi,3*pi];
c=[];
ceq=[];
wb = 0.1*pi;
for i=1:6
    c1=0;%c1可以用来设置初始角度
    c2=0;
    c3=0;
    ceq1=0;
    ceq2=0;
    ceq3=0;
    for j=1:5
        c1=c1+(coefficient_a(i,j)^2+coefficient_b(i,j)^2)^0.5/(wb*j);
        c2=c2+(coefficient_a(i,j)^2+coefficient_b(i,j)^2)^0.5;
        c3=c3+(coefficient_a(i,j)^2+coefficient_b(i,j)^2)^0.5*(wb*j);
        ceq1=ceq1+coefficient_b(i,j)/(wb*j);
        ceq2=ceq2+coefficient_a(i,j);
        ceq3=ceq3+wb*j*coefficient_b(i,j);
    end
    c1=c1-qmax(i);
    c2=c2-qdotmax(i);
    c3=c3-qdotdotmax(i);
    c=[c;c1;c2;c3];
    ceq=[ceq;ceq1;ceq2;ceq3];
end


end
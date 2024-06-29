function cond_=getCondOfY(coefficient)
    coefficient_a=reshape(coefficient(1:30),5,6)';
    %disp(coefficient_a)
    coefficient_b=reshape(coefficient(31:end),5,6)';
    %disp(coefficient_b)
    data_size=500;
    Ytilde_=[];
    parfor i=1:data_size
        q_qdot_qdotdot = getFourierTrajectory(coefficient_a,coefficient_b,0.04*i);
        Ytilde=getYtilde(q_qdot_qdotdot(:,1)+[0,0,pi/2,0,0,0]',q_qdot_qdotdot(:,2),q_qdot_qdotdot(:,3));%+[0,0,pi/2,0,0,0]'
        Ytilde_=[Ytilde_;Ytilde];
    end
    cond_=cond(Ytilde_);
end
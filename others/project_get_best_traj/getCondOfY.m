function cond_=getCondOfY(coefficient)
    coefficient_a=reshape(coefficient(1:30),5,6)';
    %disp(coefficient_a)
    coefficient_b=reshape(coefficient(31:end),5,6)';
    %disp(coefficient_b)
    data_size=1000;
    Ytilde_=[];
    for i=1:data_size
        q_qdot_qdotdot = getFourierTrajectory(coefficient_a,coefficient_b,0.004*i);
        Ytilde=getYtilde(q_qdot_qdotdot(:,1),q_qdot_qdotdot(:,2),q_qdot_qdotdot(:,3));
        Ytilde_=[Ytilde_;Ytilde];
    end
    cond_=cond(Ytilde_);
end
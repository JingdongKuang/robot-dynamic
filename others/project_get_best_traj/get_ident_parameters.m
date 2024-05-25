function ident_parameters=get_ident_parameters(torque_,data_size)
    x=load("x.txt");
    coefficient_a=reshape(x(1:30),5,6)';
    coefficient_b=reshape(x(31:end),5,6)';

    Ytilde_=[];
    for i=1:data_size
        q_qdot_qdotdot = getFourierTrajectory(coefficient_a,coefficient_b,0.004*i);
        Ytilde=getYtilde(q_qdot_qdotdot(:,1),q_qdot_qdotdot(:,2),q_qdot_qdotdot(:,3));
        Ytilde_=[Ytilde_;Ytilde];
    end
    ident_parameters=inv((Ytilde_'*Ytilde_))*Ytilde_'*torque_;

end
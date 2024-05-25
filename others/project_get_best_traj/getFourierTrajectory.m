function q_qdot_qdotdot = getFourierTrajectory(coefficient_a,coefficient_b,t)
    wb = 0.2*pi; % 2*pi/T
    q_qdot_qdotdot = zeros(6, 3); % 初始化输出矩阵
    
    for j = 1:6
        q_qdot_qdotdot(j, 1) = 0;
        q_qdot_qdotdot(j, 2) = 0;
        q_qdot_qdotdot(j, 3) = 0;

        for k = 1:5
            q_qdot_qdotdot(j, 1) = q_qdot_qdotdot(j, 1) + coefficient_a(j, k) * sin(k * wb * t) / (wb * k) - coefficient_b(j, k) * cos(k * wb * t) / (wb * k);
            q_qdot_qdotdot(j, 2) = q_qdot_qdotdot(j, 2) + coefficient_a(j, k) * cos(k * wb * t) + coefficient_b(j, k) * sin(k * wb * t);
            q_qdot_qdotdot(j, 3) = q_qdot_qdotdot(j, 3) - coefficient_a(j, k) * k * wb * sin(k * wb * t) + coefficient_b(j, k) * k * wb * cos(k * wb * t);
        end
    end
end


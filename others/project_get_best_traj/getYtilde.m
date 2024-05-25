function Ytilde = getYtilde(q, q_dot, q_dot_dot)
    
    MDH_Table=[ 0,		0.22,	0,		0;
			    pi / 2, 0,		0,		pi/2;
			    0,		0,		0.455,	0;
			    0,		0.495,	0,		pi / 2;
			    pi,		0,		0,		pi / 2;
			    pi,		0.1565, 0,		pi/2];
    z = [0; 0; 1];
    w = zeros(3, 7);
    epsilon = zeros(3, 7);
    a = zeros(3, 7);
    a(:,1)=[0;0;9.81];
    linear_independent = [1; 2; 3; 5; 6; 7; 8];
    Ri_iminus1 = cell(1, 6);
    p_i_tilde_star = zeros(3, 6);
    p_i_dash_star = zeros(3, 6);
    A = cell(1, 6);
    U = zeros(6, 10);
    Ti = cell(1, 5);
    Ytilde = zeros(6, 36);

    for i = 2:7
        T = DH2Trans(MDH_Table(i - 1, 1) + q(i - 1), MDH_Table(i - 1, 2), MDH_Table(i - 1, 3), MDH_Table(i - 1, 4));
        Ri_iminus1{i-1} = T(1:3, 1:3)';
        p_i_dash_star(:, i - 1) = T(1:3, 4);
        p_i_tilde_star(:, i - 1) = Ri_iminus1{i-1} * p_i_dash_star(:, i - 1);
        w(:, i) = Ri_iminus1{i-1} * w(:, i - 1) + z * q_dot(i - 1);
        epsilon(:, i) = Ri_iminus1{i-1} * epsilon(:, i - 1) + cross(w(:, i), z * q_dot(i - 1)) + z * q_dot_dot(i - 1);
        a(:, i) = Ri_iminus1{i-1} * (a(:, i - 1) + cross(epsilon(:, i - 1), p_i_dash_star(:, i - 1)) + cross(w(:, i - 1), cross(w(:, i - 1), p_i_dash_star(:, i - 1))));

        % Extracting dynamic parameters to form matrix A, 6x10
        A{i - 1} = zeros(6, 10);
        S_epsilon = Operator_S(epsilon(:, i));
        S_w = Operator_S(w(:, i));
        S_a = Operator_S(a(:, i));
        K_epsilon = Operator_K(epsilon(:, i));
        K_w = Operator_K(w(:, i));

        A{i - 1}(1:3, 7:9) = S_epsilon + S_w * S_w;
        A{i - 1}(1:3, 10) = a(:, i);
        A{i - 1}(4:6, 1:6) = K_epsilon + S_w * K_w;
        A{i - 1}(4:6, 7:9) = -S_a;
    end

    % Calculate Ti
    for i = 1:5
        S_p_iPlus1_tilde_star = Operator_S(p_i_dash_star(:, i + 1));
        T = DH2Trans(MDH_Table(i + 1, 1) + q(i + 1), MDH_Table(i + 1, 2), MDH_Table(i + 1, 3), MDH_Table(i + 1, 4));
        Ri_iplus1 = T(1:3, 1:3);
        Ti{i} = zeros(6, 6);
        Ti{i}(1:3, 1:3) = Ri_iplus1;
        Ti{i}(4:6, 1:3) = S_p_iPlus1_tilde_star * Ri_iplus1;
        Ti{i}(4:6, 4:6) = Ri_iplus1;
    end

    for i = 1:6
        for j = i:6
            kk = [0, 0, 0, 0, 0, 1];
            if j == i
                U = A{j};
            else
                U = A{j};
                for k = j-1:-1:i
                    U = Ti{k} * U;
                end
            end
            Yij = kk * U;
            %disp(Yij);
            if i == 1 && j == 1
                Ytilde(1, 1) = Yij(6);
            else
                %Yij(1)=Yij(1)-Yij(4);
                Ytilde(i, (7 * j - 12):(7 * j-6)) = Yij(linear_independent);

            end
        end
    end
end

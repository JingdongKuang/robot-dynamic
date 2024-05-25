function K = Operator_K(input)
    K = [input(1), input(2), input(3), 0, 0, 0;
         0, input(1), 0, input(2), input(3), 0;
         0, 0, input(1), 0, input(2), input(3)];
end

function H = TransMatrix(a, alpha, d, theta)

% DH2H - 计算DH参数对应的齐次变换矩阵

% 输入参数：
% a: DH参数中的a
% alpha: DH参数中的alpha
% d: DH参数中的d
% theta: DH参数中的theta

% 输出参数：
% H: 4 x 4 的齐次变换矩阵

H = [cos(theta), -sin(theta)*cos(alpha), sin(theta)*sin(alpha), a*cos(theta);
     sin(theta), cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta);
     0,          sin(alpha),             cos(alpha),            d;
     0,          0,                      0,                     1];
end

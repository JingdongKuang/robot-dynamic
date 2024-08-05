#include "POE_kinematics.h"

//指数积求解正解
MatrixXd POE_kinematics::forward_kinematics(VectorXd q) {
	MatrixXd T = MatrixXd::Identity(4, 4);
	return T;
}
//指数积求解逆解
MatrixXd POE_kinematics::inv_kinematics(MatrixXd T) {
	VectorXd q = VectorXd::Zero(6);
	return q;
}
//指数积求解雅克比矩阵
MatrixXd POE_kinematics::Jacobian(MatrixXd T, VectorXd q) {
	MatrixXd J = MatrixXd::Zero(6, 6);
	return J;
}
//罗德里斯公式求旋转矩阵
Matrix3d POE_kinematics::w2R(Vector3d w, double theta) {
	Matrix3d R = Matrix3d::Identity();
	Matrix3d w_hat = Operator_S(w);
	R = R + sin(theta) * w_hat + (1 - cos(theta)) * w_hat * w_hat;
	return R;
}
//旋转矩阵求旋转轴,需要注意，此处当theta为pi的整数倍数时，会出现奇异，这种奇异不可避免
//同时，我们限制theta的范围为0到pi，因此，当theta为0时，旋转轴为任意向量
Vector4d POE_kinematics::R2w_and_theta(Matrix3d R) {
	double theta;
	VectorXd w_and_theta=VectorXd::Zero(4);
	Matrix3d I = Matrix3d::Identity();
	if (R == I) {
		return w_and_theta;
	}
	else if (R.trace() == -1) {
		double temp = 1/pow((R(0, 0) + 1) * 2,0.5);
		w_and_theta(0) = (1 + R(0, 0)) / temp;
		w_and_theta(1) = R(1, 0) / temp;
		w_and_theta(2) = R(2, 0) / temp;
		w_and_theta(3) = PI;
		return w_and_theta;
	}
	else{
		w_and_theta(3) =acos((R.trace() - 1) / 2);
		w_and_theta(0) = (R(2, 1) - R(1, 2)) / (2 * sin(w_and_theta(3)));
		w_and_theta(1) = (R(0, 2) - R(2, 0)) / (2 * sin(w_and_theta(3)));
		w_and_theta(2) = (R(1, 0) - R(0, 1)) / (2 * sin(w_and_theta(3)));
		return w_and_theta;
	}


}
//so3
Matrix3d POE_kinematics::Operator_S(Vector3d w) {

	Matrix3d w_hat;
	w_hat <<	0,		-w(2),	w(1),
				w(2),	0,		-w(0),
				-w(1),	w(0),	0;
	return w_hat;
}
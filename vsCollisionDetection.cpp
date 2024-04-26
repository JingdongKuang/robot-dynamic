// vsCollisionDetection.cpp: 定义应用程序的入口点。
//

#include "vsCollisionDetection.h"
using namespace Eigen;
using namespace std;


//rot 4d


Matrix4d rotx4d(const double angle)

{

	Matrix4d res = Matrix4d::Identity();

	res(1, 1) = cos(angle);

	res(1, 2) = -sin(angle);

	res(2, 1) = sin(angle);

	res(2, 2) = cos(angle);

	return res;

}


Matrix4d roty4d(const double angle)

{

	Matrix4d res = Matrix4d::Identity();

	res(0, 0) = cos(angle);

	res(0, 2) = sin(angle);

	res(2, 0) = -sin(angle);

	res(2, 2) = cos(angle);

	return res;

}


Matrix4d rotz4d(const double angle)

{

	Matrix4d res = Matrix4d::Identity();

	res(0, 0) = cos(angle);

	res(0, 1) = -sin(angle);

	res(1, 0) = sin(angle);

	res(1, 1) = cos(angle);

	return res;

}



Matrix4d trans(const Vector3d trans)
{
	Matrix4d res = Matrix4d::Identity();

	res(0, 3) = trans(0);

	res(1, 3) = trans(1);

	res(2, 3) = trans(2);

	return res;

}

Matrix4d RobotDynamics::DH2Trans(const double theta, const double d, const double a, const double alpha)
{

	Matrix4d T;
	T << cos(theta), -sin(theta) * cos(alpha), sin(theta)* sin(alpha), a* cos(theta),
		sin(theta), cos(theta)* cos(alpha), -cos(theta) * sin(alpha), a* sin(theta),
		0, sin(alpha), cos(alpha), d,
		0, 0, 0, 1;

	return T;
}




/***********************************************
*DESCRIPTION:
*	用来计算新坐标系下的惯性张量，需要注意的是，此处采用的是先旋转，再平移的方式。
*
*
*INPUT:
*	R 为旋转向量，3*1向量，比如（1，2，3）意味着标准单位阵，现在有个新坐标系，x在负x方向，y在负z方向，z在y方向，那么R就是（-1，-3，2）
*
*	旋转的情况I1=R13*I3*R13'
*	平移的情况：I1=I3+m*(P1*P1'*I-P1'*P1) P1为质心在新坐标系下的位置，m为质量
*
*   P 为平移，指的是质心在新坐标系下的位置，3*1向量
*	I 为在质心处的惯性张量，3*3矩阵
*
* OUTPUT:
*	返回新的惯性张量I_new
*
************************************************/

Matrix3d InertialTensor_Trans(const Matrix3d I, const Vector3i R = Vector3i(1, 2, 3), const Vector3d P = Vector3d(0, 0, 0), const double m = 0) {
	//先旋转
	Matrix3d Rot = Matrix3d::Zero(3, 3);
	Rot(abs(R(0)) - 1, 0) = R(1) / abs(R(1));
	Rot(abs(R(1)) - 1, 1) = R(1) / abs(R(1));
	Rot(abs(R(2)) - 1, 2) = R(2) / abs(R(2));

	Matrix3d I_new = Rot * I * Rot.transpose();

	//再平移
	I_new += m * (P.dot(P) * Matrix3d::Identity() - P * P.transpose());

	return I_new;
}
//变换坐标系后，获得新的重心位置，指的是当前坐标系下的重心位置
Vector3d getNewGpositon(const Vector3d p, const Vector3i R = Vector3d(1, 2, 3)) {
	Matrix3d Rot = Matrix3d::Zero(3, 3);
	Rot(abs(R(0)) - 1, 0) = R(1) / abs(R(1));
	Rot(abs(R(1)) - 1, 1) = R(1) / abs(R(1));
	Rot(abs(R(2)) - 1, 2) = R(2) / abs(R(2));

	Vector3d p_new = Rot * p;
	return p_new;
}




Matrix4d RobotDynamics::inverseHomogeneousTransform(const Matrix4d& transform) {
	// Extract rotation matrix and translation vector from homogeneous transform
	Matrix3d rotationMatrix = transform.block<3, 3>(0, 0);
	Vector3d translationVector = transform.block<3, 1>(0, 3);

	// Compute inverse of rotation matrix
	Matrix3d inverseRotationMatrix = rotationMatrix.transpose();

	// Compute inverse of translation vector
	Vector3d inverseTranslationVector = -inverseRotationMatrix * translationVector;

	// Construct the inverse homogeneous transform
	Matrix4d inverseTransform = Eigen::Matrix4d::Identity();
	inverseTransform.block<3, 3>(0, 0) = inverseRotationMatrix;
	inverseTransform.block<3, 1>(0, 3) = inverseTranslationVector;

	return inverseTransform;
}

/***********************************************
* INPUT:
*
*
* OUTPUT:
*
*
* DESCRIPTION:
*	计算质量矩阵
***********************************************/
MatrixXd RobotDynamics::getMassMatrix(const VectorXd& q) {
	getTransMatrix(q);
	Vector4d temp;
	double q1 = q(0);
	double q2 = q(1);
	double q3 = q(2);
	double q4 = q(3);
	double q5 = q(4);
	double q6 = q(5);

	//计算线速度雅可比矩阵
	MatrixXd Jv[6];
	MatrixXd Jw[6];

	// 创建12个全零矩阵
	for (int i = 0; i < 6; i++) {
		Jv[i] = MatrixXd::Zero(3, 6); // 初始化矩阵都是3x6的全零矩阵
		Jw[i] = MatrixXd::Zero(3, 6);
	}
	MatrixXd JwTest(3, 6);

	//Jv6	
	temp = Q * T06 * r[5];
	Jv[5].block(0, 0, 3, 1) = temp.head(3);
	temp = T01 * Q * T12 * T23 * T34 * T45 * T56 * r[5];
	Jv[5].block(0, 1, 3, 1) = temp.head(3);
	temp = T02 * Q * T23 * T34 * T45 * T56 * r[5];
	Jv[5].block(0, 2, 3, 1) = temp.head(3);
	temp = T03 * Q * T34 * T45 * T56 * r[5];
	Jv[5].block(0, 3, 3, 1) = temp.head(3);
	temp = T04 * Q * T45 * T56 * r[5];
	Jv[5].block(0, 4, 3, 1) = temp.head(3);
	temp = T05 * Q * T56 * r[5];
	Jv[5].block(0, 5, 3, 1) = temp.head(3);
	//Jv5
	temp = Q * T05 * r[4];
	Jv[4].block(0, 0, 3, 1) = temp.head(3);
	temp = T01 * Q * T12 * T23 * T34 * T45 * r[4];
	Jv[4].block(0, 1, 3, 1) = temp.head(3);
	temp = T02 * Q * T23 * T34 * T45 * r[4];
	Jv[4].block(0, 2, 3, 1) = temp.head(3);
	temp = T03 * Q * T34 * T45 * r[4];
	Jv[4].block(0, 3, 3, 1) = temp.head(3);
	temp = T04 * Q * T45 * r[4];
	Jv[4].block(0, 4, 3, 1) = temp.head(3);
	//Jv4
	temp = Q * T04 * r[3];
	Jv[3].block(0, 0, 3, 1) = temp.head(3);
	temp = T01 * Q * T12 * T23 * T34 * r[3];
	Jv[3].block(0, 1, 3, 1) = temp.head(3);
	temp = T02 * Q * T23 * T34 * r[3];
	Jv[3].block(0, 2, 3, 1) = temp.head(3);
	temp = T03 * Q * T34 * r[3];
	Jv[3].block(0, 3, 3, 1) = temp.head(3);
	//Jv3
	temp = Q * T03 * r[2];
	Jv[2].block(0, 0, 3, 1) = temp.head(3);
	temp = T01 * Q * T12 * T23 * r[2];
	Jv[2].block(0, 1, 3, 1) = temp.head(3);
	temp = T02 * Q * T23 * r[2];
	Jv[2].block(0, 2, 3, 1) = temp.head(3);
	//Jv2
	temp = Q * T02 * r[1];
	Jv[1].block(0, 0, 3, 1) = temp.head(3);
	temp = T01 * Q * T12 * r[1];
	;	Jv[1].block(0, 1, 3, 1) = temp.head(3);
	//Jv1
	temp = Q * T01 * r[0];
	Jv[0].block(0, 0, 3, 1) = temp.head(3);





	//计算角速度雅可比矩阵
	//Jw1
	Jw[0](2, 0) = 1;
	//Jw2
	Jw[1] = Jw[0];
	Jw[1].block(0, 1, 3, 1) = T01.block(0, 2, 3, 1);
	//Jw3
	Jw[2] = Jw[1];
	Jw[2].block(0, 2, 3, 1) = T02.block(0, 2, 3, 1);
	//Jw4
	Jw[3] = Jw[2];
	Jw[3].block(0, 3, 3, 1) = T03.block(0, 2, 3, 1);
	//Jw5
	Jw[4] = Jw[3];
	Jw[4].block(0, 4, 3, 1) = T04.block(0, 2, 3, 1);
	//Jw6
	Jw[5] = Jw[4];
	Jw[5].block(0, 5, 3, 1) = T05.block(0, 2, 3, 1);



	//质量矩阵
	MatrixXd MassMatrix(6, 6);
	MassMatrix = Jv[0].transpose() * m[0] * Jv[0] + Jw[0].transpose() * I[0] * Jw[0] +
		Jv[1].transpose() * m[1] * Jv[1] + Jw[1].transpose() * I[1] * Jw[1] +
		Jv[2].transpose() * m[2] * Jv[2] + Jw[2].transpose() * I[2] * Jw[2] +
		Jv[3].transpose() * m[3] * Jv[3] + Jw[3].transpose() * I[3] * Jw[3] +
		Jv[4].transpose() * m[4] * Jv[4] + Jw[4].transpose() * I[4] * Jw[4] +
		Jv[5].transpose() * m[5] * Jv[5] + Jw[5].transpose() * I[5] * Jw[5];
	return MassMatrix;
}


/***********************************************
* INPUT:
*	q为关节角度，6*1向量
*	q_dot为关节角速度，6*1向量
*	q_dot_dot为关节角加速度，6*1向量
*	M1为前一时刻的质量矩阵，6*6矩阵
*	M2为当前时刻的质量矩阵，6*6矩阵
*   T为采样周期
*
*
* OUTPUT:
*
*
* DESCRIPTION:
*	基于拉格朗日方程的动力学模型，用于计算关节力矩
***********************************************/
VectorXd RobotDynamics::getTorque(const VectorXd& q, const VectorXd& q_dot, const VectorXd& q_dot_dot, const MatrixXd& M_now) {
	//计算拉格朗日方程中对q求偏导的项，通过数值方法求得
	double delta = 0.00000001;//数值微分的步长
	//每个分量的小变化
	MatrixXd q_ = MatrixXd::Identity(6, 6);
	//计算每个分量的小变化对应的质量矩阵
	MatrixXd dMdq[6];
	for (int i = 0; i < 6; i++) {
		dMdq[i] = MatrixXd::Zero(6, 6); // 初始化矩阵都是3x6的全零矩阵
	}

	dMdq[0] = (getMassMatrix(q + delta * q_.col(0)) - M_now) / delta;
	dMdq[1] = (getMassMatrix(q + delta * q_.col(1)) - M_now) / delta;
	dMdq[2] = (getMassMatrix(q + delta * q_.col(2)) - M_now) / delta;
	dMdq[3] = (getMassMatrix(q + delta * q_.col(3)) - M_now) / delta;
	dMdq[4] = (getMassMatrix(q + delta * q_.col(4)) - M_now) / delta;
	dMdq[5] = (getMassMatrix(q + delta * q_.col(5)) - M_now) / delta;


	VectorXd dKdq(6);
	dKdq(0) = 0.5 * q_dot.transpose() * dMdq[0] * q_dot;
	dKdq(1) = 0.5 * q_dot.transpose() * dMdq[1] * q_dot;
	dKdq(2) = 0.5 * q_dot.transpose() * dMdq[2] * q_dot;
	dKdq(3) = 0.5 * q_dot.transpose() * dMdq[3] * q_dot;
	dKdq(4) = 0.5 * q_dot.transpose() * dMdq[4] * q_dot;
	dKdq(5) = 0.5 * q_dot.transpose() * dMdq[5] * q_dot;

	VectorXd Mdot_multi_qdot(6);
	MatrixXd temp_Matrix = MatrixXd::Zero(6, 6);
	/*
	for (int k = 0; k < 6; k++) {
		for(int i = 0; i < 6; i++){
			temp_Matrix = dMdq[i].row(k);
		}
		Mdot_multi_qdot(k) = q_dot.transpose() * temp_Matrix * q_dot;
	}*/
	for (int i = 0; i < 6; i++) {
		temp_Matrix += dMdq[i] * q_dot(i);
	}
	Mdot_multi_qdot = temp_Matrix * q_dot;

	getTransMatrix(q);
	//计算重力矩
	MatrixXd G0 = m[0] * g * Q * T01 * r[0] + \
		m[1] * g * Q * T01 * (inverseHomogeneousTransform(T01) * T02) * r[1] + \
		m[2] * g * Q * T01 * (inverseHomogeneousTransform(T01) * T03) * r[2] + \
		m[3] * g * Q * T01 * (inverseHomogeneousTransform(T01) * T04) * r[3] + \
		m[4] * g * Q * T01 * (inverseHomogeneousTransform(T01) * T05) * r[4] + \
		m[5] * g * Q * T01 * (inverseHomogeneousTransform(T01) * T06) * r[5];
	G(0) = G0(0, 0);
	MatrixXd G1 = m[1] * g * T01 * Q * T12 * r[1] + \
		m[2] * g * T01 * Q * T12 * T23 * r[2] + \
		m[3] * g * T01 * Q * T12 * (inverseHomogeneousTransform(T02) * T04) * r[3] + \
		m[4] * g * T01 * Q * T12 * (inverseHomogeneousTransform(T02) * T05) * r[4] + \
		m[5] * g * T01 * Q * T12 * (inverseHomogeneousTransform(T02) * T06) * r[5];
	G(1) = G1(0, 0);
	MatrixXd G2 = m[2] * g * T02 * Q * T23 * r[2] + \
		m[3] * g * T02 * Q * T23 * (inverseHomogeneousTransform(T03) * T04) * r[3] + \
		m[4] * g * T02 * Q * T23 * (inverseHomogeneousTransform(T03) * T05) * r[4] + \
		m[5] * g * T02 * Q * T23 * (inverseHomogeneousTransform(T03) * T06) * r[5];
	G(2) = G2(0, 0);
	MatrixXd G3 = m[3] * g * T03 * Q * T34 * r[3] + \
		m[4] * g * T03 * Q * T34 * (inverseHomogeneousTransform(T04) * T05) * r[4] + \
		m[5] * g * T03 * Q * T34 * (inverseHomogeneousTransform(T04) * T06) * r[5];
	G(3) = G3(0, 0);
	MatrixXd G4 = m[4] * g * T04 * Q * T45 * r[4] + \
		m[5] * g * T04 * Q * T45 * T56 * r[5];
	G(4) = G4(0, 0);
	MatrixXd G5 = m[5] * g * T05 * Q * T56 * r[5];
	G(5) = G5(0, 0);




	//将上面计算结果拼接成最终输出的理论力矩
	VectorXd torque(6);

	torque = M_now * q_dot_dot + Mdot_multi_qdot - dKdq + G;

	return torque;
}

//获取齐次变换矩阵的微分
MatrixXd RobotDynamics::get_T_Derivative_of_time() {
	
	
	MatrixXd dTdq(6, 6);
	return dTdq;
}
//获取齐次变换矩阵对于某个关节角度的微分

MatrixXd RobotDynamics::get_T_Derivative_of_qi(int T0i, int qi) {
	MatrixXd dTdq(6, 1);

	return dTdq;
}

VectorXd RobotDynamics::getGravity() {
	return G;
}

void RobotDynamics::getTransMatrix(const VectorXd& q) {

	T01 = DH2Trans(DH_Table(0, 0) + q(0), DH_Table(0, 1), DH_Table(0, 2), DH_Table(0, 3));
	T12 = DH2Trans(DH_Table(1, 0) + q(1), DH_Table(1, 1), DH_Table(1, 2), DH_Table(1, 3));
	T23 = DH2Trans(DH_Table(2, 0) + q(2), DH_Table(2, 1), DH_Table(2, 2), DH_Table(2, 3));
	T34 = DH2Trans(DH_Table(3, 0) + q(3), DH_Table(3, 1), DH_Table(3, 2), DH_Table(3, 3));
	T45 = DH2Trans(DH_Table(4, 0) + q(4), DH_Table(4, 1), DH_Table(4, 2), DH_Table(4, 3));
	T56 = DH2Trans(DH_Table(5, 0) + q(5), DH_Table(5, 1), DH_Table(5, 2), DH_Table(5, 3));
	T02 = T01 * T12;
	T03 = T02 * T23;
	T04 = T03 * T34;
	T05 = T04 * T45;
	T06 = T05 * T56;
}
//输出齐次变换矩阵
void RobotDynamics::coutTransMatrix() {
	cout << "T01: " << endl << T01 << endl;
	cout << "T02: " << endl << T02 << endl;
	cout << "T03: " << endl << T03 << endl;
	cout << "T04: " << endl << T04 << endl;
	cout << "T05: " << endl << T05 << endl;
	cout << "T06: " << endl << T06 << endl;
}
//见机器人动力学与控制p.78
VectorXd RobotDynamics::getTorque_Newton_Euler(const VectorXd& q, const VectorXd& q_dot, const VectorXd& q_dot_dot){
	Vector3d w[7];//公式1-53
	w[0] <<0,0,0; // 初始化全零向量
	Vector3d epsilon[7];//公式1-74
	epsilon[0] << 0, 0, 0; // 初始化全零向量
	Vector3d a[7];//公式2-79
	a[0] << 0, 0, 0; // 初始化全零向量
	Matrix3d Ri_iminus1;

	Vector3d z;
	z<< 0, 0, 1;

	for (int i = 1; i < 7; i++) {
		Matrix4d T = DH2Trans(DH_Table(i-1, 0) + q(i-1), DH_Table(i-1, 1), DH_Table(i-1, 2), DH_Table(i-1, 3));
		Ri_iminus1 = T.block(0, 0, 3, 3).transpose();
		w[i] = Ri_iminus1*(w[i - 1] + z*q_dot(i-1));//公式1-53
		epsilon[i] = Ri_iminus1 * (epsilon[i - 1] + w[i - 1].cross(z*q_dot(i-1)) + z*q_dot_dot(i-1));//公式1-74
		//a[i] = ;//公式2-79
	}


	VectorXd torque(6);
	return torque;
}


int main()
{
	RobotDynamics objRobotDynamics;
	ChebyshefFilter objChebyshefFilter;


	MatrixXd M_past(6, 6);
	MatrixXd M_now(6, 6);
	VectorXd q_now(6);
	q_now << 0.004, 0.004, 0.004, 0.004, 0.004, 0.004;

	VectorXd v(6);
	v << 2, 2, 2, 2, 2, 2;

	VectorXd a(6);
	a << 2, 2, 2, 2, 2, 2;

	M_now = objRobotDynamics.getMassMatrix(q_now);

	VectorXd torque(6);
	torque = objRobotDynamics.getTorque(q_now, v, a, M_now);
	VectorXd G = objRobotDynamics.getGravity();
	cout << "torque: " << endl << torque << endl;

	double data = objChebyshefFilter.dofilter(1);
	cout << "data: " << data << endl;
	return 0;
}
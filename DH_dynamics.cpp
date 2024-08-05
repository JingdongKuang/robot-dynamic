// vsCollisionDetection.cpp: 定义应用程序的入口点。
//

#include "DH_dynamics.h"
using namespace Eigen;
using namespace std;
//求矩阵条件数
/*
double matrix_condition_number(const MatrixXd& mat)
{
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(mat);

	// 获取奇异值
	Eigen::VectorXd singularValues = svd.singularValues();

	// 计算条件数
	double conditionNumber = singularValues.maxCoeff() / singularValues.minCoeff();
	return conditionNumber;
}*/
//输出矩阵txt
void matrix_save_txt(const MatrixXd& mat, string filename)
{
	ofstream outfile(filename, ios::trunc);
	outfile << mat;
	outfile.close();
	cout<< "Save " << filename << " successfully!" << endl;
}

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

Matrix4d RobotDynamics::DH2Trans(const double theta, const double d, const double a, const double alpha,const int table_type)
{
	Matrix4d T;//SDH
	if (table_type == 0) {
		T <<	cos(theta),		-sin(theta) * cos(alpha),		sin(theta)* sin(alpha),		a* cos(theta),
				sin(theta),		cos(theta)* cos(alpha),			-cos(theta) * sin(alpha),	a* sin(theta),
				0,				sin(alpha),						cos(alpha),					d,
				0,				0,								0,							1;
	}
	else {//MDH
		T <<	cos(theta),					-sin(theta),				0,					a,
				sin(theta)*cos(alpha),		cos(theta)*cos(alpha),		-sin(alpha),		-sin(alpha)*d,
				sin(theta)*sin(alpha),		cos(theta)*sin(alpha),		cos(alpha),			cos(alpha)*d,
				0,							0,							0,					1;
	}
	

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
	double delta = 0.000000000001;//数值微分的步长
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

	for (int i = 0; i < 6; i++) {
		temp_Matrix += dMdq[i] * q_dot(i);
	}
	Mdot_multi_qdot = temp_Matrix * q_dot;

	getTransMatrix(q);
	//计算重力矩
	MatrixXd G0 = m[0] * g.transpose() * Q * T01 * r[0] + \
		m[1] * g.transpose() * Q * T01 * (inverseHomogeneousTransform(T01) * T02) * r[1] + \
		m[2] * g.transpose() * Q * T01 * (inverseHomogeneousTransform(T01) * T03) * r[2] + \
		m[3] * g.transpose() * Q * T01 * (inverseHomogeneousTransform(T01) * T04) * r[3] + \
		m[4] * g.transpose() * Q * T01 * (inverseHomogeneousTransform(T01) * T05) * r[4] + \
		m[5] * g.transpose() * Q * T01 * (inverseHomogeneousTransform(T01) * T06) * r[5];
	G(0) = G0(0, 0);
	MatrixXd G1 = m[1] * g.transpose() * T01 * Q * T12 * r[1] + \
		m[2] * g.transpose() * T01 * Q * T12 * T23 * r[2] + \
		m[3] * g.transpose() * T01 * Q * T12 * (inverseHomogeneousTransform(T02) * T04) * r[3] + \
		m[4] * g.transpose() * T01 * Q * T12 * (inverseHomogeneousTransform(T02) * T05) * r[4] + \
		m[5] * g.transpose() * T01 * Q * T12 * (inverseHomogeneousTransform(T02) * T06) * r[5];
	G(1) = G1(0, 0);
	MatrixXd G2 = m[2] * g.transpose() * T02 * Q * T23 * r[2] + \
		m[3] * g.transpose() * T02 * Q * T23 * (inverseHomogeneousTransform(T03) * T04) * r[3] + \
		m[4] * g.transpose() * T02 * Q * T23 * (inverseHomogeneousTransform(T03) * T05) * r[4] + \
		m[5] * g.transpose() * T02 * Q * T23 * (inverseHomogeneousTransform(T03) * T06) * r[5];
	G(2) = G2(0, 0);
	MatrixXd G3 = m[3] * g.transpose() * T03 * Q * T34 * r[3] + \
		m[4] * g.transpose() * T03 * Q * T34 * (inverseHomogeneousTransform(T04) * T05) * r[4] + \
		m[5] * g.transpose() * T03 * Q * T34 * (inverseHomogeneousTransform(T04) * T06) * r[5];
	G(3) = G3(0, 0);
	MatrixXd G4 = m[4] * g.transpose() * T04 * Q * T45 * r[4] + \
		m[5] * g.transpose() * T04 * Q * T45 * T56 * r[5];
	G(4) = G4(0, 0);
	MatrixXd G5 = m[5] * g.transpose() * T05 * Q * T56 * r[5];
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
//算子K
MatrixXd RobotDynamics::Operator_K(const Vector3d input) {
	MatrixXd Operator_K(3, 6);
	Operator_K <<	input(0),	input(1),	input(2),	0,			0,			0,
					0,			input(0),	0,			input(1),	input(2),	0,
					0,			0,			input(0),	0,			input(1),	input(2);
	return Operator_K;
}
//算子S
Matrix3d RobotDynamics::Operator_S(const Vector3d input) {
	Matrix3d Operator_S;
	Operator_S <<	0,			-input(2),	input(1),
					input(2),	0,			-input(0),
					-input(1),	input(0),	0;
	return Operator_S;
}

VectorXd RobotDynamics::getGravity() {
	return G;
}
//
void RobotDynamics::getTransMatrix(const VectorXd& q,const int table_type) {
	MatrixXd Table;
	if (table_type == 0) {
		Table = DH_Table;
	}
	else {
		Table = MDH_Table;
	}
	T01 = DH2Trans(Table(0, 0) + q(0), Table(0, 1), Table(0, 2), Table(0, 3), table_type);
	T12 = DH2Trans(Table(1, 0) + q(1), Table(1, 1), Table(1, 2), Table(1, 3), table_type);
	T23 = DH2Trans(Table(2, 0) + q(2), Table(2, 1), Table(2, 2), Table(2, 3), table_type);
	T34 = DH2Trans(Table(3, 0) + q(3), Table(3, 1), Table(3, 2), Table(3, 3), table_type);
	T45 = DH2Trans(Table(4, 0) + q(4), Table(4, 1), Table(4, 2), Table(4, 3), table_type);
	T56 = DH2Trans(Table(5, 0) + q(5), Table(5, 1), Table(5, 2), Table(5, 3), table_type);
	T02 = T01 * T12;
	T03 = T02 * T23;
	T04 = T03 * T34;
	T05 = T04 * T45;
	T06 = T05 * T56;
}


//输出齐次变换矩阵
void RobotDynamics::coutTransMatrix() {
	cout << "T01: " << endl << T01 << endl;
	cout << "T12: " << endl << T12 << endl;
	cout << "T23: " << endl << T23 << endl;
	cout << "T34: " << endl << T34 << endl;
	cout << "T45: " << endl << T45 << endl;
	cout << "T56: " << endl << T56 << endl;
	cout << "T02: " << endl << T02 << endl;
	cout << "T03: " << endl << T03 << endl;
	cout << "T04: " << endl << T04 << endl;
	cout << "T05: " << endl << T05 << endl;
	cout << "T06: " << endl << T06 << endl;
}
/*
*见机器人动力学与控制p.78
*newton-euler方法
*/
VectorXd RobotDynamics::getTorque_Newton_Euler(const VectorXd& q, const VectorXd& q_dot, const VectorXd& q_dot_dot){
	VectorXd torque(6);
	Vector3d z;
	z << 0, 0, 1;
	Vector3d w[7];//公式1-53
	w[0] <<0,0,0; // 初始化全零向量
	Vector3d epsilon[7];//公式1-74
	epsilon[0] << 0, 0, 0; // 初始化全零向量
	Vector3d a[7];//公式2-79
	a[0] = g.head(3); // 初始化全零向量
	Vector3d a_C_i[7];//公式2-76
	Matrix3d Ri_iminus1;
	Vector3d p_i_tilde_star[6];
	Vector3d rc;
	for (int i = 1; i < 7; i++) {
		Matrix4d T = DH2Trans(DH_Table(i-1, 0) + q(i-1), DH_Table(i-1, 1), DH_Table(i-1, 2), DH_Table(i-1, 3));
		Ri_iminus1 = T.block(0, 0, 3, 3).transpose();
		p_i_tilde_star[i - 1]<< DH_Table(i - 1, 2), DH_Table(i - 1, 1)*sin(DH_Table(i - 1, 3)), DH_Table(i - 1, 1)*cos(DH_Table(i - 1, 3));
		w[i] = Ri_iminus1*(w[i - 1] + z*q_dot(i-1));
		epsilon[i] = Ri_iminus1 * (epsilon[i - 1] + w[i - 1].cross(z * q_dot(i - 1))  + z*q_dot_dot(i-1));
		a[i] = Ri_iminus1 * a[i - 1] + epsilon[i].cross(p_i_tilde_star[i-1]) + w[i].cross(w[i].cross(p_i_tilde_star[i-1]));
		rc = r[i - 1].head(3);
		a_C_i[i] = a[i] + epsilon[i].cross(rc) + w[i].cross(w[i].cross(rc));
	}
	Vector3d n[7];
	n[6] << 0, 0, 0;
	Vector3d N[6];
	Vector3d f[7];
	f[6] << 0, 0, 0;
	Vector3d F[6];
	Matrix3d Ri_iplus1;
	Matrix4d T;
	for (int i = 5; i >= 0; i--) {
		rc = r[i].head(3);
		F[i] = m[i] * a_C_i[i + 1];
		N[i] = I_C[i] * epsilon[i+1] + w[i+1].cross(I_C[i] * w[i+1]);
		if (i == 5) {
			f[i] = F[i];
			n[i] = N[i]+ p_i_tilde_star[i].cross(f[i])+rc.cross(F[i]);
		}
		else {
			T = DH2Trans(DH_Table(i+1, 0) + q(i+1), DH_Table(i+1, 1), DH_Table(i+1, 2), DH_Table(i+1, 3));
			Ri_iplus1 = T.block(0, 0, 3, 3);
			f[i] = Ri_iplus1 * f[i + 1] + F[i];
			n[i] = Ri_iplus1 * n[i + 1] + N[i] + p_i_tilde_star[i].cross(f[i]) + rc.cross(F[i]);
		}
		Matrix4d T1 = DH2Trans(DH_Table(i, 0) + q(i), DH_Table(i, 1), DH_Table(i, 2), DH_Table(i, 3));
		Matrix3d Ri_iminus1_ = T1.block(0, 0, 3, 3).transpose();
		torque(i) = (Ri_iminus1_ * z).transpose()* n[i];
	}
	return torque;
}

//p100,2-151至2-156
VectorXd RobotDynamics::getTorque_Newton_Euler_MDH(const VectorXd& q, const VectorXd& q_dot, const VectorXd& q_dot_dot) {
	VectorXd torque(6);
	Vector3d z;
	z << 0, 0, 1;
	Vector3d w[7];//公式1-53
	w[0] << 0, 0, 0; // 初始化全零向量
	Vector3d epsilon[7];//公式1-74
	epsilon[0] << 0, 0, 0; // 初始化全零向量
	Vector3d a[7];//公式2-79
	a[0] = g.head(3); // 初始化全零向量
	Vector3d a_C_i[7];//公式2-76
	Matrix3d Ri_iminus1;
	Vector3d p_i_tilde_star[6];
	Vector3d p_i_dash_star[6];
	Vector3d rc;
	for (int i = 1; i < 7; i++) {
		Matrix4d T = DH2Trans(MDH_Table(i - 1, 0) + q(i - 1), MDH_Table(i - 1, 1), MDH_Table(i - 1, 2), MDH_Table(i - 1, 3),1);
		Ri_iminus1 = T.block(0, 0, 3, 3).transpose();
		p_i_dash_star[i - 1] = T.block(0,3,3,1);
		//cout<<"p_i_dash_star["<<i-1<<"]:"<<p_i_dash_star[i-1].transpose()<<endl;
		p_i_tilde_star[i - 1] = Ri_iminus1* p_i_dash_star[i-1];
		//cout << "p_i_tilde_star[" << i - 1 << "]:" << p_i_tilde_star[i - 1].transpose() << endl;
		w[i] = Ri_iminus1 * w[i - 1] + z * q_dot(i - 1);
		//cout << "w[" << i << "]:" << w[i].transpose() << endl;
		epsilon[i] = Ri_iminus1 * epsilon[i - 1] + w[i].cross(z * q_dot(i - 1)) + z * q_dot_dot(i - 1);
		//cout << "epsilon[" << i << "]:" << epsilon[i].transpose() << endl;
		a[i] = Ri_iminus1 * (a[i - 1] + epsilon[i-1].cross(p_i_dash_star[i - 1]) + w[i-1].cross(w[i-1].cross(p_i_dash_star[i - 1])));
		//cout << "a[" << i << "]:" << a[i].transpose() << endl;
		rc = r_MDH[i - 1].head(3);
		a_C_i[i] = a[i] + epsilon[i].cross(rc) + w[i].cross(w[i].cross(rc));
		//cout << "a_C_i[" << i << "]:" << a_C_i[i].transpose() << endl<<endl;
	}
	Vector3d n[7];
	n[6] << 0, 0, 0;
	Vector3d N[6];
	Vector3d f[7];
	f[6] << 0, 0, 0;
	Vector3d F[6];
	Matrix3d Ri_iplus1;
	Matrix4d T;
	for (int i = 5; i >= 0; i--) {
		rc = r_MDH[i].head(3);
		F[i] = m[i] * a_C_i[i + 1];
		N[i] = I_C[i] * epsilon[i + 1] + w[i + 1].cross(I_C[i] * w[i + 1]);
		if (i == 5) {
			f[i] = F[i];
			n[i] = N[i] + rc.cross(F[i]);
		}
		else {
			T = DH2Trans(MDH_Table(i + 1, 0) + q(i + 1), MDH_Table(i + 1, 1), MDH_Table(i + 1, 2), MDH_Table(i + 1, 3),1);
			Ri_iplus1 = T.block(0, 0, 3, 3);
			f[i] = Ri_iplus1 * f[i + 1] + F[i];
			n[i] = N[i] + rc.cross(F[i]) + Ri_iplus1*(p_i_tilde_star[i+1].cross(f[i+1]) + n[i + 1]);
		}
		//cout << "F[" << i << "]:" << F[i].transpose() << endl;
		//cout << "f[" << i << "]:" << f[i].transpose() << endl;
		//cout << "N[" << i << "]:" << N[i].transpose() << endl;
		//cout<< "n[" << i << "]:" << n[i].transpose() << endl;
		
		torque(i) = z.transpose()*n[i];
	}
	return torque;
}


MatrixXd RobotDynamics::getYtilde(const VectorXd& q, const VectorXd& q_dot, const VectorXd& q_dot_dot) {
	Vector3d z;
	z << 0, 0, 1;
	Vector3d w[7];//公式1-53
	w[0] << 0, 0, 0; // 初始化全零向量
	Vector3d epsilon[7];//公式1-74
	epsilon[0] << 0, 0, 0; // 初始化全零向量
	Vector3d a[7];//公式2-79
	a[0] = g.head(3); // 初始化全零向量
	VectorXi linear_independent(7);//线性无关项的列的索引
	linear_independent << 0,1,2,4,5,6,7;
	Matrix3d Ri_iminus1[6];
	Matrix3d Ri_iplus1;
	Vector3d p_i_tilde_star[6];
	Vector3d p_i_dash_star[6];
	MatrixXd A[6];
	MatrixXd U(6, 10);
	MatrixXd Ti[5];
	//MatrixXd Y=MatrixXd::Zero(6,60);
	MatrixXd Ytilde=MatrixXd::Zero(6, 36);
	for (int i = 1; i < 7; i++) {
		Matrix4d T = DH2Trans(MDH_Table(i - 1, 0) + q(i - 1), MDH_Table(i - 1, 1), MDH_Table(i - 1, 2), MDH_Table(i - 1, 3), 1);
		Ri_iminus1[i-1] = T.block(0, 0, 3, 3).transpose();
		p_i_dash_star[i - 1] = T.block(0, 3, 3, 1);//p_i-1_i
		p_i_tilde_star[i - 1] = Ri_iminus1[i-1] * p_i_dash_star[i - 1];
		w[i] = Ri_iminus1[i-1] * w[i - 1] + z * q_dot(i - 1);
		epsilon[i] = Ri_iminus1[i-1] * epsilon[i - 1] + w[i].cross(z * q_dot(i - 1)) + z * q_dot_dot(i - 1);
		a[i] = Ri_iminus1[i-1] * (a[i - 1] + epsilon[i - 1].cross(p_i_dash_star[i - 1]) + w[i - 1].cross(w[i - 1].cross(p_i_dash_star[i - 1])));


		//将动力学参数从耦合中提取出来，构成矩阵A，形状为6*10，公式2-137
		A[i - 1] = MatrixXd::Zero(6, 10);
		//Astar[i - 1] = MatrixXd::Zero(6, 6);
		//将epsilon和w,p*还有a转化为算子形式
		Matrix3d S_epsilon = Operator_S(epsilon[i]);
		Matrix3d S_w = Operator_S(w[i]);
		Matrix3d S_a = Operator_S(a[i]);
		MatrixXd K_epsilon = Operator_K(epsilon[i]);
		MatrixXd K_w = Operator_K(w[i]);

		//将这些算子填入矩阵A
		A[i - 1].block(0, 6, 3, 3) = S_epsilon+S_w*S_w;
		A[i - 1].block(0, 9, 3, 1) = a[i];
		A[i - 1].block(3, 0, 3, 6) = K_epsilon+ S_w* K_w;
		A[i - 1].block(3, 6, 3, 3) = -S_a;
		//A[i - 1].block(3, 9, 3, 1) = MatrixXd::Zero(3,1);
		//cout<<"A"<<i-1<<endl<<A[i - 1] << endl;
		
	}
	//计算Ti
	for (int i = 0; i < 5;i++) {
		Matrix3d S_p_iPlus1_tilde_star = Operator_S(p_i_dash_star[i+1]);
		MatrixXd T = DH2Trans(MDH_Table(i+1, 0) + q(i+1), MDH_Table(i+1, 1), MDH_Table(i+1, 2), MDH_Table(i+1, 3), 1);
		Ri_iplus1 = T.block(0, 0, 3, 3);
		Ti[i] = MatrixXd::Zero(6, 6);//2-143
		Ti[i].block(0, 0, 3, 3) = Ri_iplus1;
		Ti[i].block(3, 0, 3, 3) = S_p_iPlus1_tilde_star * Ri_iplus1;
		Ti[i].block(3, 3, 3, 3) = Ri_iplus1;
		//cout << "T" << i << endl << Ti[i] << endl;
	}
	for (int i = 0; i < 6; i++) {
		//求取大Y矩阵
		for (int j = i; j < 6; j++) {
			RowVectorXd kk(6);
			kk<<0, 0, 0, 0,0,1;//2-145
			if (j == i) {
				U = A[j];
			}
			else {
				U = A[j];
				for (int k = j-1; k >=i; k--) {
					U = Ti[k] * U;
				}
			}
			RowVectorXd Yij = kk * U;//2-145*/
			//cout << "kk*Ustar:"<<endl << kk*Ustar << endl;
			if (i == 0&&j==0) {
				Ytilde(0,0) = Yij(5);
			}
			else {
				Ytilde.block(i, 7 * j-6, 1, 7) = Yij(Eigen::placeholders::all, linear_independent);//2-145*/
			}
			//cout << "Ytilde"<<endl<<Ytilde << endl;
		}
	}
	

	//matrix_save_txt(Ytilde, "Ytilde.txt");
	return Ytilde;
}

void RobotDynamics::identifyDynamicsParameters(const MatrixXd& Ytilde_, const VectorXd& torque_) {
	//int num=torque_.size();
	VectorXd DynamicsParametersparameters(36);
	//MatrixXd Ytilde__;
	//VectorXd torque__;
	DynamicsParametersparameters= (Ytilde_.transpose()*Ytilde_).inverse()* Ytilde_.transpose()*torque_;
	/*//2-147将Ytilde_转化为上三角模式，同样的，将torque_转化为上三角模式
	MatrixXd TransRowMatrix=MatrixXd::Zero(num, num);//用来将输入矩阵转化为上三角模式
	for (int j = 0; j < 6; j++) {
		for (int i = 0; i < num / 6; i++) {
			TransRowMatrix(j * num / 6+i,6*i+j) = 1;
		}
	}
	Ytilde__= TransRowMatrix*Ytilde_;
	torque__ = TransRowMatrix*torque_;
	VectorXd torque_tilde(6);
	for (int i = num / 6-1; i > -1; i--) {
		if (i == 5) {
			torque_tilde=torque__.segment(6 * i, 6);
		}else {
			torque_tilde = torque__.segment(6 * i, 6);
		}
		
	}*/
	cout<<"DynamicsParametersparameters: "<<endl<<DynamicsParametersparameters<<endl;
}

/*
*基于傅里叶级数的激励轨迹
*输入为时间t
*返回值为关节角度，关节角速度，关节角加速度，形状为6*3的矩阵
*/
MatrixXd RobotDynamics::getFourierTrajectory(const double&t) {
	double wb=0.2*PI;//2*PI/T
	MatrixXd q_qdot_qdotdot(6, 3);//y
	for (int j = 0; j < 6; j++) {
		q_qdot_qdotdot(j, 0) = 0;
		q_qdot_qdotdot(j, 1) = 0;
		q_qdot_qdotdot(j, 2) = 0;
		for (int k = 0; k < 5; k++) {
			q_qdot_qdotdot(j, 0) += coefficient_a(j, k ) * sin((k + 1) * wb * t) / (wb * (k + 1)) - coefficient_b(j, k) * cos((k + 1) * wb * t) / (wb * (k + 1));
			q_qdot_qdotdot(j, 1) += coefficient_a(j, k) * cos((k + 1) * wb * t) + coefficient_b(j, k) * sin((k + 1) * wb * t);
			q_qdot_qdotdot(j, 2) += -coefficient_a(j, k) * (k + 1) * wb * sin((k + 1) * wb * t) + coefficient_b(j, k) * (k + 1) * wb * cos((k + 1) * wb * t);
		}
	}
	return q_qdot_qdotdot;
}

VectorXd ZeroPhaseAverageFilter(int width, const VectorXd& data) {
	VectorXd filter = VectorXd::Ones(width) / width;
	double floorWidthdivedTwo= floor(width/2);
	VectorXd filteredData_(data.head(floorWidthdivedTwo));
	for (int i = floorWidthdivedTwo; i < data.size() - floorWidthdivedTwo; i++) {
		double temp = data.segment(i - floorWidthdivedTwo, width).transpose() * filter;
		filteredData_.conservativeResize(filteredData_.size() + 1);
		filteredData_(filteredData_.size() - 1) = temp;
	}
	filteredData_.conservativeResize(filteredData_.size() + floorWidthdivedTwo);
	filteredData_.tail(floorWidthdivedTwo) = data.tail(floorWidthdivedTwo);

	VectorXd data_ = filteredData_.reverse();
	VectorXd filteredData(data_.head(floorWidthdivedTwo));
	for (int i = floorWidthdivedTwo; i < data_.size() - floorWidthdivedTwo; ++i) {
		VectorXd temp = data_.segment(i - floorWidthdivedTwo, width).transpose() * filter;
		filteredData.conservativeResize(filteredData.size() + 1);
		filteredData(filteredData.size() - 1) = temp(0);
	}
	filteredData.conservativeResize(filteredData.size() + floorWidthdivedTwo);
	filteredData.tail(floorWidthdivedTwo) = data_.tail(floorWidthdivedTwo);

	return filteredData.reverse();
}

VectorXd RobotDynamics::least_sq_fit(VectorXd x, VectorXd y, int degree) {

	VectorXd k(degree + 1);

	int len = y.size();

	MatrixXd A = MatrixXd::Ones(len, degree + 1);

	for (int i = 0; i < len; i++) {

		for (int j = 0; j < degree; j++) {

			A(i, j) = pow(x(i), degree - j);

		}

	}

	k = (A.transpose() * A).inverse() * A.transpose() * y;

	return k;

}

double RobotDynamics::polynomial_Curve(double x, VectorXd k) {

	double y = 0;

	int len = k.size();

	for (int i = 0; i < len; i++) {

		y += k(i) * pow(x, len - i - 1);

	}

	return y;

}
//int main()
//{
//	RobotDynamics objRobotDynamics;
//	ChebyshefFilter objChebyshefFilter;
//	/***/
//	VectorXd q_now(6);
//	q_now << 3, 2, 2, 2, 2, 2;
//	VectorXd v(6);
//	v << 3, 2, 3, 3, 3, 3;
//	VectorXd a(6);
//	a << 3, 3, 4, 3, 3, 3;
//	VectorXd id_para(36);
//	id_para << 2.67879856071727,
//		- 2.58582663613946,
//		- 1.24808071877980e-05,
//		0.126011316400632,
//		5.32383350523857e-06,
//		2.59235962858469,
//		5.83200516840000,
//		3.87622999990234e-05,
//		2.40053551817202,
//		1.64364357836122e-05,
//		1.80074602398606e-05,
//		8.90741163689975e-05,
//		2.40729679557488,
//		0.000425330399999091,
//		- 3.82639832350000,
//		0.0119575940894542,
//		9.37445760748612e-06,
//		- 3.48449540342605e-05,
//		0.0385644843765673,
//		0.0153019351156814,
//		7.74372000001619e-05,
//		- 0.135984347100000,
//		0.0173936972509635,
//		- 1.14668181056324e-06,
//		1.00111813058912e-06,
//		- 2.43869688341877e-05,
//		0.0193853104106020,
//		5.32410000006756e-05,
//		- 0.171403313700000,
//		5.49537439283182e-05,
//		2.57684213678805e-07,
//		- 2.71500523660317e-07,
//		- 1.47844718488277e-05,
//		0.000585850701614120,
//		- 1.81281000005556e-05,
//		- 0.000640526199999256;
//	
//	VectorXd torque(6);
//	torque = objRobotDynamics.getTorque_Newton_Euler_MDH(q_now, v, a);
//	MatrixXd Ytilde = objRobotDynamics.getYtilde(q_now, v, a);
//
//
//
//
//
//	/**/
//	int imax = 1000;
//	MatrixXd q_qdot_qdotdot(6, 3);
//	MatrixXd Ytilde_(6*imax, 36);
//	VectorXd torque_(6*imax);
//	for (int i = 0; i < imax; i++) {
//		q_qdot_qdotdot = objRobotDynamics.getFourierTrajectory(i*0.004);//MatrixXd::Random(6, 3);
//		Ytilde_.block(6*i, 0, 6, 36) = objRobotDynamics.getYtilde(q_qdot_qdotdot.col(0), q_qdot_qdotdot.col(1), q_qdot_qdotdot.col(2));
//		torque_.segment(6*i,6) = objRobotDynamics.getTorque_Newton_Euler_MDH(q_qdot_qdotdot.col(0), q_qdot_qdotdot.col(1), q_qdot_qdotdot.col(2));
//		//cout << "torque_" << i << ": " << endl << torque_.segment(6 * i, 6) << endl;
//	}
//	objRobotDynamics.identifyDynamicsParameters(Ytilde_, torque_);
//	/**/
//	
//	return 0;
//}
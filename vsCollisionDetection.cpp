// vsCollisionDetection.cpp: 定义应用程序的入口点。
//

#include "vsCollisionDetection.h"
using namespace Eigen;
using namespace std;

#define PI 3.1415926535897931
#define toRad 0.01745329252
//rot 4d

Matrix4d inverseHomogeneousTransform(const Matrix4d& transform) {
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

Matrix4d DH2Trans(const double theta, const double d, const double a, const double alpha) 
{

	Matrix4d T;
	T << cos(theta), -sin(theta) * cos(alpha), sin(theta)* sin(alpha),   a* cos(theta),
		 sin(theta), cos(theta)* cos(alpha),   -cos(theta) * sin(alpha), a* sin(theta),
		 0,          sin(alpha),               cos(alpha),               d,
		 0,          0,                        0,                        1;

	return T;
}

/***********************************************
* INPUT:
*	theta=[0,0,0,0,0,0],为1*6的向量
*	DH_table为6*4的矩阵，列分别为offset，d，a，alpha，此处仿照matlab设计
* 
* OUTPUT:
*	返回正向计算矩阵
* 
************************************************/
Matrix4d forward_kine(const VectorXd theta, const MatrixXd DH_Table) 
{
	
	Matrix4d T01 = DH2Trans(DH_Table(0, 0) + theta(0), DH_Table(0, 1), DH_Table(0, 2), DH_Table(0, 3));
	Matrix4d T12 = DH2Trans(DH_Table(1, 0) + theta(1), DH_Table(1, 1), DH_Table(1, 2), DH_Table(1, 3));
	Matrix4d T23 = DH2Trans(DH_Table(2, 0) + theta(2), DH_Table(2, 1), DH_Table(2, 2), DH_Table(2, 3));
	Matrix4d T34 = DH2Trans(DH_Table(3, 0) + theta(3), DH_Table(3, 1), DH_Table(3, 2), DH_Table(3, 3));
	Matrix4d T45 = DH2Trans(DH_Table(4, 0) + theta(4), DH_Table(4, 1), DH_Table(4, 2), DH_Table(4, 3));
	Matrix4d T56 = DH2Trans(DH_Table(5, 0) + theta(5), DH_Table(5, 1), DH_Table(5, 2), DH_Table(5, 3));

	Matrix4d kine = T01*T12*T23*T34*T45*T56;
	return kine;
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

Matrix3d InertialTensor_Trans(const Matrix3d I, const Vector3i R = Vector3i (1,2,3), const Vector3d P = Vector3d (0,0,0),const double m=0) {
	//先旋转
	Matrix3d Rot = Matrix3d::Zero(3, 3);
	Rot(abs(R(0))-1, 0) = R(1) / abs(R(1));
	Rot(abs(R(1))-1, 1) = R(1) / abs(R(1));
	Rot(abs(R(2))-1, 2) = R(2) / abs(R(2));
	
	Matrix3d I_new = Rot * I * Rot.transpose();

	//再平移
	I_new+=m*(P.dot(P)*Matrix3d::Identity()-P*P.transpose());

	return I_new;
}
//变换坐标系后，获得新的重心位置，指的是当前坐标系下的重心位置
Vector3d getNewGpositon(const Vector3d p,const Vector3i R=Vector3d(1,2,3)) {
	Matrix3d Rot = Matrix3d::Zero(3, 3);
	Rot(abs(R(0)) - 1, 0) = R(1) / abs(R(1));
	Rot(abs(R(1)) - 1, 1) = R(1) / abs(R(1));
	Rot(abs(R(2)) - 1, 2) = R(2) / abs(R(2));
	
	Vector3d p_new = Rot * p;
	return p_new;
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
MatrixXd getMassMatrix(const VectorXd q) {
	
	double q1 = q(0);
	double q2 = q(1);
	double q3 = q(2);
	double q4 = q(3);
	double q5 = q(4);
	double q6 = q(5);

	//计算线速度雅可比矩阵
	MatrixXd Jv1(3, 6);
	MatrixXd Jv2(3, 6);
	MatrixXd Jv3(3, 6);
	MatrixXd Jv4(3, 6);
	MatrixXd Jv5(3, 6);
	MatrixXd Jv6(3, 6);
	
	Jv6 << cos(q1) * (sin(q4) * (-0.12246 * sin(q5) + cos(q5) * (-0.00003 * cos(q6) + 0.00106 * sin(q6))) + cos(q4) * (-0.00106 * cos(q6) - 0.00003 * sin(q6))) - sin(q1) * (sin(q2) * (-0.455 + sin(q3) * (-0.495 - 0.12246 * cos(q5) +
		0.00003 * cos(q6) * sin(q5) - 0.00106 * sin(q5) * sin(q6)) + cos(q3) * (sin(q4) * (-0.00106 * cos(q6) - 0.00003 * sin(q6)) + cos(q4) * (0.00003 * cos(q5) * cos(q6) + 0.12246 * sin(q5) - 0.00106 * cos(q5) * sin(q6)))) +
		cos(q2) * (cos(q3) * (0.495 + 0.12246 * cos(q5) - 0.00003 * cos(q6) * sin(q5) + 0.00106 * sin(q5) * sin(q6)) + sin(q3) * (sin(q4) * (-0.00106 * cos(q6) - 0.00003 * sin(q6)) + cos(q4) * (0.00003 * cos(q5) * cos(q6) +
			0.12246 * sin(q5) - 0.00106 * cos(q5) * sin(q6))))), cos(q1)* (cos(q2) * (-0.455 + sin(q3) * (-0.495 - 0.12246 * cos(q5) + 0.00003 * cos(q6) * sin(q5) - 0.00106 * sin(q5) * sin(q6)) +
				cos(q3) * (sin(q4) * (-0.00106 * cos(q6) - 0.00003 * sin(q6)) + cos(q4) * (0.00003 * cos(q5) * cos(q6) + 0.12246 * sin(q5) - 0.00106 * cos(q5) * sin(q6)))) - sin(q2) * (cos(q3) * (0.495 +
					0.12246 * cos(q5) - 0.00003 * cos(q6) * sin(q5) + 0.00106 * sin(q5) * sin(q6)) + sin(q3) * (sin(q4) * (-0.00106 * cos(q6) - 0.00003 * sin(q6)) + cos(q4) * (0.00003 * cos(q5) * cos(q6) + 0.12246 * sin(q5) -
						0.00106 * cos(q5) * sin(q6))))), cos(q1)* (cos(q2) * (-sin(q3) * (0.495 + 0.12246 * cos(q5) - 0.00003 * cos(q6) * sin(q5) + 0.00106 * sin(q5) * sin(q6)) + cos(q3) * (sin(q4) * (-0.00106 * cos(q6) - 0.00003 * sin(q6)) +
							cos(q4) * (0.00003 * cos(q5) * cos(q6) + 0.12246 * sin(q5) - 0.00106 * cos(q5) * sin(q6)))) + sin(q2) * (cos(q3) * (-0.495 - 0.12246 * cos(q5) + 0.00003 * cos(q6) * sin(q5) - 0.00106 * sin(q5) * sin(q6)) -
								sin(q3) * (sin(q4) * (-0.00106 * cos(q6) - 0.00003 * sin(q6)) + cos(q4) * (0.00003 * cos(q5) * cos(q6) + 0.12246 * sin(q5) - 0.00106 * cos(q5) * sin(q6))))), sin(q1)* (cos(q4) * (-0.12246 * sin(q5) +
									cos(q5) * (-0.00003 * cos(q6) + 0.00106 * sin(q6))) - sin(q4) * (-0.00106 * cos(q6) - 0.00003 * sin(q6))) + cos(q1) * (cos(q3) * sin(q2) * (cos(q4) * (-0.00106 * cos(q6) - 0.00003 * sin(q6)) -
										sin(q4) * (0.00003 * cos(q5) * cos(q6) + 0.12246 * sin(q5) - 0.00106 * cos(q5) * sin(q6))) + cos(q2) * sin(q3) * (cos(q4) * (-0.00106 * cos(q6) - 0.00003 * sin(q6)) - sin(q4) * (0.00003 * cos(q5) * cos(q6) +
											0.12246 * sin(q5) - 0.00106 * cos(q5) * sin(q6)))), sin(q1)* sin(q4)* (-0.12246 * cos(q5) - sin(q5) * (-0.00003 * cos(q6) + 0.00106 * sin(q6))) + cos(q1) * (sin(q2) * (sin(q3) * (0.00003 * cos(q5) * cos(q6) +
												0.12246 * sin(q5) - 0.00106 * cos(q5) * sin(q6)) + cos(q3) * cos(q4) * (0.12246 * cos(q5) - 0.00003 * cos(q6) * sin(q5) + 0.00106 * sin(q5) * sin(q6))) + cos(q2) * (cos(q3) * (-0.00003 * cos(q5) * cos(q6) -
													0.12246 * sin(q5) + 0.00106 * cos(q5) * sin(q6)) + cos(q4) * sin(q3) * (0.12246 * cos(q5) - 0.00003 * cos(q6) * sin(q5) + 0.00106 * sin(q5) * sin(q6)))), sin(q1)* (cos(q5) * sin(q4) * (0.00106 * cos(q6) +
														0.00003 * sin(q6)) + cos(q4) * (-0.00003 * cos(q6) + 0.00106 * sin(q6))) + cos(q1) * (sin(q2) * (sin(q3) * (-0.00106 * cos(q6) * sin(q5) - 0.00003 * sin(q5) * sin(q6)) +
															cos(q3) * (sin(q4) * (-0.00003 * cos(q6) + 0.00106 * sin(q6)) + cos(q4) * (-0.00106 * cos(q5) * cos(q6) - 0.00003 * cos(q5) * sin(q6)))) + cos(q2) * (cos(q3) * (0.00106 * cos(q6) * sin(q5) +
																0.00003 * sin(q5) * sin(q6)) + sin(q3) * (sin(q4) * (-0.00003 * cos(q6) + 0.00106 * sin(q6)) + cos(q4) * (-0.00106 * cos(q5) * cos(q6) - 0.00003 * cos(q5) * sin(q6))))),

		-sin(q1) * (sin(q4) * (0.12246 * sin(q5) + cos(q5) * (0.00003 * cos(q6) - 0.00106 * sin(q6))) + cos(q4) * (0.00106 * cos(q6) + 0.00003 * sin(q6))) + cos(q1) * (sin(q2) * (-0.455 + sin(q3) * (-0.495 - 0.12246 * cos(q5) + 0.00003 * cos(q6) * sin(q5) - 0.00106 * sin(q5) * sin(q6)) + cos(q3) * (sin(q4) * (-0.00106 * cos(q6) - 0.00003 * sin(q6)) + cos(q4) * (0.00003 * cos(q5) * cos(q6) + 0.12246 * sin(q5) - 0.00106 * cos(q5) * sin(q6)))) + cos(q2) * (cos(q3) * (0.495 + 0.12246 * cos(q5) - 0.00003 * cos(q6) * sin(q5) + 0.00106 * sin(q5) * sin(q6)) + sin(q3) * (sin(q4) * (-0.00106 * cos(q6) - 0.00003 * sin(q6)) + cos(q4) * (0.00003 * cos(q5) * cos(q6) + 0.12246 * sin(q5) - 0.00106 * cos(q5) * sin(q6))))),
		sin(q1)* (cos(q2) * (-0.455 + sin(q3) * (-0.495 - 0.12246 * cos(q5) + 0.00003 * cos(q6) * sin(q5) - 0.00106 * sin(q5) * sin(q6)) + cos(q3) * (sin(q4) * (-0.00106 * cos(q6) - 0.00003 * sin(q6)) + cos(q4) * (0.00003 * cos(q5) * cos(q6) + 0.12246 * sin(q5) - 0.00106 * cos(q5) * sin(q6)))) - sin(q2) * (cos(q3) * (0.495 + 0.12246 * cos(q5) - 0.00003 * cos(q6) * sin(q5) + 0.00106 * sin(q5) * sin(q6)) + sin(q3) * (sin(q4) * (-0.00106 * cos(q6) - 0.00003 * sin(q6)) + cos(q4) * (0.00003 * cos(q5) * cos(q6) + 0.12246 * sin(q5) - 0.00106 * cos(q5) * sin(q6))))),
		sin(q1)* (cos(q2) * (-sin(q3) * (0.495 + 0.12246 * cos(q5) - 0.00003 * cos(q6) * sin(q5) + 0.00106 * sin(q5) * sin(q6)) + cos(q3) * (sin(q4) * (-0.00106 * cos(q6) - 0.00003 * sin(q6)) + cos(q4) * (0.00003 * cos(q5) * cos(q6) + 0.12246 * sin(q5) - 0.00106 * cos(q5) * sin(q6)))) + sin(q2) * (cos(q3) * (-0.495 - 0.12246 * cos(q5) + 0.00003 * cos(q6) * sin(q5) - 0.00106 * sin(q5) * sin(q6)) - sin(q3) * (sin(q4) * (-0.00106 * cos(q6) - 0.00003 * sin(q6)) + cos(q4) * (0.00003 * cos(q5) * cos(q6) + 0.12246 * sin(q5) - 0.00106 * cos(q5) * sin(q6))))),
		cos(q1)* (cos(q4) * (0.12246 * sin(q5) + cos(q5) * (0.00003 * cos(q6) - 0.00106 * sin(q6))) - sin(q4) * (0.00106 * cos(q6) + 0.00003 * sin(q6))) + sin(q1) * (cos(q3) * sin(q2) * (cos(q4) * (-0.00106 * cos(q6) - 0.00003 * sin(q6)) - sin(q4) * (0.00003 * cos(q5) * cos(q6) + 0.12246 * sin(q5) - 0.00106 * cos(q5) * sin(q6))) + cos(q2) * sin(q3) * (cos(q4) * (-0.00106 * cos(q6) - 0.00003 * sin(q6)) - sin(q4) * (0.00003 * cos(q5) * cos(q6) + 0.12246 * sin(q5) - 0.00106 * cos(q5) * sin(q6)))),
		cos(q1)* sin(q4)* (0.12246 * cos(q5) - sin(q5) * (0.00003 * cos(q6) - 0.00106 * sin(q6))) + sin(q1) * (sin(q2) * (sin(q3) * (0.00003 * cos(q5) * cos(q6) + 0.12246 * sin(q5) - 0.00106 * cos(q5) * sin(q6)) + cos(q3) * cos(q4) * (0.12246 * cos(q5) - 0.00003 * cos(q6) * sin(q5) + 0.00106 * sin(q5) * sin(q6))) + cos(q2) * (cos(q3) * (-0.00003 * cos(q5) * cos(q6) - 0.12246 * sin(q5) + 0.00106 * cos(q5) * sin(q6)) + cos(q4) * sin(q3) * (0.12246 * cos(q5) - 0.00003 * cos(q6) * sin(q5) + 0.00106 * sin(q5) * sin(q6)))),
		cos(q1)* (cos(q4) * (0.00003 * cos(q6) - 0.00106 * sin(q6)) + cos(q5) * sin(q4) * (-0.00106 * cos(q6) - 0.00003 * sin(q6))) + sin(q1) * (sin(q2) * (sin(q3) * (-0.00106 * cos(q6) * sin(q5) - 0.00003 * sin(q5) * sin(q6)) + cos(q3) * (sin(q4) * (-0.00003 * cos(q6) + 0.00106 * sin(q6)) + cos(q4) * (-0.00106 * cos(q5) * cos(q6) - 0.00003 * cos(q5) * sin(q6)))) + cos(q2) * (cos(q3) * (0.00106 * cos(q6) * sin(q5) + 0.00003 * sin(q5) * sin(q6)) + sin(q3) * (sin(q4) * (-0.00003 * cos(q6) + 0.00106 * sin(q6)) + cos(q4) * (-0.00106 * cos(q5) * cos(q6) - 0.00003 * cos(q5) * sin(q6))))),


		0,0.00003 * cos(q2) * cos(q4) * cos(q5) * cos(q6) * sin(q3) - 0.00106 * cos(q2) * cos(q6) * sin(q3) * sin(q4) + 0.12246 * cos(q2) * cos(q4) * sin(q3) * sin(q5) - 0.00106 * cos(q2) * cos(q4) * cos(q5) * sin(q3) * sin(q6) - 0.00003 * cos(q2) * sin(q3) * sin(q4) * sin(q6) + cos(q2) * cos(q3) * (0.495 + 0.12246 * cos(q5) - 0.00003 * cos(q6) * sin(q5) + 0.00106 * sin(q5) * sin(q6)) - sin(q2) * (0.455 + sin(q3) * (0.495 + 0.12246 * cos(q5) - 0.00003 * cos(q6) * sin(q5) + 0.00106 * sin(q5) * sin(q6)) + cos(q3) * (sin(q4) * (0.00106 * cos(q6) + 0.00003 * sin(q6)) + cos(q4) * (-0.00003 * cos(q5) * cos(q6) - 0.12246 * sin(q5) + 0.00106 * cos(q5) * sin(q6)))),
		0.00003 * cos(q3) * cos(q4) * cos(q5) * cos(q6) * sin(q2) - 0.00106 * cos(q3) * cos(q6) * sin(q2) * sin(q4) + 0.12246 * cos(q3) * cos(q4) * sin(q2) * sin(q5) - 0.00106 * cos(q3) * cos(q4) * cos(q5) * sin(q2) * sin(q6) - 0.00003 * cos(q3) * sin(q2) * sin(q4) * sin(q6) - sin(q2) * sin(q3) * (0.495 + 0.12246 * cos(q5) - 0.00003 * cos(q6) * sin(q5) + 0.00106 * sin(q5) * sin(q6)) + cos(q2) * (cos(q3) * (0.495 + 0.12246 * cos(q5) - 0.00003 * cos(q6) * sin(q5) + 0.00106 * sin(q5) * sin(q6)) - sin(q3) * (sin(q4) * (0.00106 * cos(q6) + 0.00003 * sin(q6)) + cos(q4) * (-0.00003 * cos(q5) * cos(q6) - 0.12246 * sin(q5) + 0.00106 * cos(q5) * sin(q6)))),
		-0.00106 * cos(q4) * cos(q6) * sin(q2) * sin(q3) - 0.00003 * cos(q5) * cos(q6) * sin(q2) * sin(q3) * sin(q4) - 0.12246 * sin(q2) * sin(q3) * sin(q4) * sin(q5) - 0.00003 * cos(q4) * sin(q2) * sin(q3) * sin(q6) + 0.00106 * cos(q5) * sin(q2) * sin(q3) * sin(q4) * sin(q6) + cos(q2) * cos(q3) * (cos(q4) * (0.00106 * cos(q6) + 0.00003 * sin(q6)) - sin(q4) * (-0.00003 * cos(q5) * cos(q6) - 0.12246 * sin(q5) + 0.00106 * cos(q5) * sin(q6))),
		0.12246 * cos(q4) * cos(q5) * sin(q2) * sin(q3) - 0.00003 * cos(q4) * cos(q6) * sin(q2) * sin(q3) * sin(q5) + 0.00106 * cos(q4) * sin(q2) * sin(q3) * sin(q5) * sin(q6) + cos(q3) * sin(q2) * (-0.00003 * cos(q5) * cos(q6) - 0.12246 * sin(q5) + 0.00106 * cos(q5) * sin(q6)) + cos(q2) * (sin(q3) * (-0.00003 * cos(q5) * cos(q6) - 0.12246 * sin(q5) + 0.00106 * cos(q5) * sin(q6)) + cos(q3) * cos(q4) * (-0.12246 * cos(q5) + 0.00003 * cos(q6) * sin(q5) - 0.00106 * sin(q5) * sin(q6))),
		-0.00106 * cos(q4) * cos(q5) * cos(q6) * sin(q2) * sin(q3) - 0.00003 * cos(q6) * sin(q2) * sin(q3) * sin(q4) - 0.00003 * cos(q4) * cos(q5) * sin(q2) * sin(q3) * sin(q6) + 0.00106 * sin(q2) * sin(q3) * sin(q4) * sin(q6) + cos(q3) * sin(q2) * (0.00106 * cos(q6) * sin(q5) + 0.00003 * sin(q5) * sin(q6)) + cos(q2) * (sin(q3) * (0.00106 * cos(q6) * sin(q5) + 0.00003 * sin(q5) * sin(q6)) + cos(q3) * (sin(q4) * (0.00003 * cos(q6) - 0.00106 * sin(q6)) + cos(q4) * (0.00106 * cos(q5) * cos(q6) + 0.00003 * cos(q5) * sin(q6))));


		
		



	Jv5 << -sin(q1) * (cos(q2) * (sin(q3) * (-0.01507 * sin(q4) + cos(q4) * (0.0002 * cos(q5) + 0.03659 * sin(q5))) + cos(q3) * (0.495 + 0.03659 * cos(q5) - 0.0002 * sin(q5))) + sin(q2) * (-0.455 + cos(q3) * (-0.01507 * sin(q4) + cos(q4) * (0.0002 * cos(q5) + 0.03659 * sin(q5))) + sin(q3) * (-0.495 - 0.03659 * cos(q5) + 0.0002 * sin(q5)))) + cos(q1) * (-0.01507 * cos(q4) + sin(q4) * (-0.0002 * cos(q5) - 0.03659 * sin(q5))),
		cos(q1)* (-sin(q2) * (sin(q3) * (-0.01507 * sin(q4) + cos(q4) * (0.0002 * cos(q5) + 0.03659 * sin(q5))) + cos(q3) * (0.495 + 0.03659 * cos(q5) - 0.0002 * sin(q5))) + cos(q2) * (-0.455 + cos(q3) * (-0.01507 * sin(q4) + cos(q4) * (0.0002 * cos(q5) + 0.03659 * sin(q5))) + sin(q3) * (-0.495 - 0.03659 * cos(q5) + 0.0002 * sin(q5)))),
		cos(q1)* (cos(q2) * (cos(q3) * (-0.01507 * sin(q4) + cos(q4) * (0.0002 * cos(q5) + 0.03659 * sin(q5))) - sin(q3) * (0.495 + 0.03659 * cos(q5) - 0.0002 * sin(q5))) + sin(q2) * (-sin(q3) * (-0.01507 * sin(q4) + cos(q4) * (0.0002 * cos(q5) + 0.03659 * sin(q5))) + cos(q3) * (-0.495 - 0.03659 * cos(q5) + 0.0002 * sin(q5)))),
		cos(q1)* (cos(q3) * sin(q2) * (-0.01507 * cos(q4) - sin(q4) * (0.0002 * cos(q5) + 0.03659 * sin(q5))) + cos(q2) * sin(q3) * (-0.01507 * cos(q4) - sin(q4) * (0.0002 * cos(q5) + 0.03659 * sin(q5)))) + sin(q1) * (0.01507 * sin(q4) + cos(q4) * (-0.0002 * cos(q5) - 0.03659 * sin(q5))),
		cos(q1)* (cos(q2) * (cos(q3) * (-0.0002 * cos(q5) - 0.03659 * sin(q5)) + cos(q4) * sin(q3) * (0.03659 * cos(q5) - 0.0002 * sin(q5))) + sin(q2) * (cos(q3) * cos(q4) * (0.03659 * cos(q5) - 0.0002 * sin(q5)) + sin(q3) * (0.0002 * cos(q5) + 0.03659 * sin(q5)))) + sin(q1) * sin(q4) * (-0.03659 * cos(q5) + 0.0002 * sin(q5)), 0,

		cos(q1)* (cos(q2) * (sin(q3) * (-0.01507 * sin(q4) + cos(q4) * (0.0002 * cos(q5) + 0.03659 * sin(q5))) + cos(q3) * (0.495 + 0.03659 * cos(q5) - 0.0002 * sin(q5))) + sin(q2) * (-0.455 + cos(q3) * (-0.01507 * sin(q4) + cos(q4) * (0.0002 * cos(q5) + 0.03659 * sin(q5))) + sin(q3) * (-0.495 - 0.03659 * cos(q5) + 0.0002 * sin(q5)))) - sin(q1) * (0.01507 * cos(q4) + sin(q4) * (0.0002 * cos(q5) + 0.03659 * sin(q5))),
		sin(q1)* (-sin(q2) * (sin(q3) * (-0.01507 * sin(q4) + cos(q4) * (0.0002 * cos(q5) + 0.03659 * sin(q5))) + cos(q3) * (0.495 + 0.03659 * cos(q5) - 0.0002 * sin(q5))) + cos(q2) * (-0.455 + cos(q3) * (-0.01507 * sin(q4) + cos(q4) * (0.0002 * cos(q5) + 0.03659 * sin(q5))) + sin(q3) * (-0.495 - 0.03659 * cos(q5) + 0.0002 * sin(q5)))),
		sin(q1)* (cos(q2) * (cos(q3) * (-0.01507 * sin(q4) + cos(q4) * (0.0002 * cos(q5) + 0.03659 * sin(q5))) - sin(q3) * (0.495 + 0.03659 * cos(q5) - 0.0002 * sin(q5))) + sin(q2) * (-sin(q3) * (-0.01507 * sin(q4) + cos(q4) * (0.0002 * cos(q5) + 0.03659 * sin(q5))) + cos(q3) * (-0.495 - 0.03659 * cos(q5) + 0.0002 * sin(q5)))),
		sin(q1)* (cos(q3) * sin(q2) * (-0.01507 * cos(q4) - sin(q4) * (0.0002 * cos(q5) + 0.03659 * sin(q5))) + cos(q2) * sin(q3) * (-0.01507 * cos(q4) - sin(q4) * (0.0002 * cos(q5) + 0.03659 * sin(q5)))) + cos(q1) * (-0.01507 * sin(q4) + cos(q4) * (0.0002 * cos(q5) + 0.03659 * sin(q5))),
		sin(q1)* (cos(q2) * (cos(q3) * (-0.0002 * cos(q5) - 0.03659 * sin(q5)) + cos(q4) * sin(q3) * (0.03659 * cos(q5) - 0.0002 * sin(q5))) + sin(q2) * (cos(q3) * cos(q4) * (0.03659 * cos(q5) - 0.0002 * sin(q5)) + sin(q3) * (0.0002 * cos(q5) + 0.03659 * sin(q5)))) + cos(q1) * sin(q4) * (0.03659 * cos(q5) - 0.0002 * sin(q5)), 0,0,

		0.495 * cos(q2 + q3) + 0.03659 * cos(q2 + q3) * cos(q5) - 0.01507 * cos(q2) * sin(q3) * sin(q4) + cos(q2) * cos(q4) * sin(q3) * (0.0002 * cos(q5) + 0.03659 * sin(q5)) - 0.0002 * cos(q2 + q3) * sin(q5) - sin(q2) * (0.455 + cos(q3) * (-0.0002 * cos(q4) * cos(q5) + 0.01507 * sin(q4) - 0.03659 * cos(q4) * sin(q5))),
		0.495 * cos(q2 + q3) + 0.03659 * cos(q2 + q3) * cos(q5) - 0.01507 * cos(q3) * sin(q2) * sin(q4) + cos(q3) * cos(q4) * sin(q2) * (0.0002 * cos(q5) + 0.03659 * sin(q5)) - 0.0002 * cos(q2 + q3) * sin(q5) - cos(q2) * sin(q3) * (-0.0002 * cos(q4) * cos(q5) + 0.01507 * sin(q4) - 0.03659 * cos(q4) * sin(q5)),
		-0.01507 * cos(q4) * sin(q2) * sin(q3) - sin(q2) * sin(q3) * sin(q4) * (0.0002 * cos(q5) + 0.03659 * sin(q5)) + cos(q2) * cos(q3) * (0.01507 * cos(q4) + 0.0002 * cos(q5) * sin(q4) + 0.03659 * sin(q4) * sin(q5)),
		-0.0002 * cos(q5) * sin(q2 + q3) + cos(q4) * sin(q2) * sin(q3) * (0.03659 * cos(q5) - 0.0002 * sin(q5)) - 0.03659 * sin(q2 + q3) * sin(q5) + cos(q2) * cos(q3) * (-0.03659 * cos(q4) * cos(q5) + 0.0002 * cos(q4) * sin(q5)),0;






	
	Jv4 <<		cos(q1) * (-0.03714 * cos(q4) + 0.00003 * sin(q4)) - 
				sin(q1) * (sin(q2) * (-0.455 - 0.00003 * cos(q3) * cos(q4) - 0.27613 * sin(q3) - 
				0.03714 * cos(q3) * sin(q4)) + 
				cos(q2) * (0.27613 * cos(q3) - 0.00003 * cos(q4) * sin(q3) - 
				0.03714 * sin(q3) * sin(q4))), 
				cos(q1) * (cos(q2) * (-0.455 - 0.00003 * cos(q3) * cos(q4) - 0.27613 * sin(q3) - 
				0.03714 * cos(q3) * sin(q4)) - 
				sin(q2) * (0.27613 * cos(q3) - 0.00003 * cos(q4) * sin(q3) - 
				0.03714 * sin(q3) * sin(q4))), 
				cos(q1) * (cos(q2) * (-0.00003 * cos(q3) * cos(q4) - 0.27613 * sin(q3) - 
				0.03714 * cos(q3) * sin(q4)) + 
				sin(q2) * (-0.27613 * cos(q3) + 0.00003 * cos(q4) * sin(q3) + 
				0.03714 * sin(q3) * sin(q4))), 
				sin(q1) * (0.00003 * cos(q4) + 0.03714 * sin(q4)) + 
				cos(q1) * (sin(q2) * (-0.03714 * cos(q3) * cos(q4) + 0.00003 * cos(q3) * sin(q4)) + 
				cos(q2) * (-0.03714 * cos(q4) * sin(q3) + 0.00003 * sin(q3) * sin(q4))), 0, 0,

				- sin(q1) * (0.03714 * cos(q4) - 0.00003 * sin(q4)) +
				cos(q1) * (sin(q2) * (-0.455 - 0.00003 * cos(q3) * cos(q4) - 0.27613 * sin(q3) -
				0.03714 * cos(q3) * sin(q4)) +
				cos(q2) * (0.27613 * cos(q3) - 0.00003 * cos(q4) * sin(q3) -
				0.03714 * sin(q3) * sin(q4))),
				sin(q1) * (cos(q2) * (-0.455 - 0.00003 * cos(q3) * cos(q4) - 0.27613 * sin(q3) -
				0.03714 * cos(q3) * sin(q4)) -
				sin(q2) * (0.27613 * cos(q3) - 0.00003 * cos(q4) * sin(q3) -
				0.03714 * sin(q3) * sin(q4))),
				sin(q1) * (cos(q2) * (-0.00003 * cos(q3) * cos(q4) - 0.27613 * sin(q3) -
				0.03714 * cos(q3) * sin(q4)) +
				sin(q2) * (-0.27613 * cos(q3) + 0.00003 * cos(q4) * sin(q3) +
				0.03714 * sin(q3) * sin(q4))),
				cos(q1) * (-0.00003 * cos(q4) - 0.03714 * sin(q4)) +
				sin(q1) * (sin(q2) * (-0.03714 * cos(q3) * cos(q4) + 0.00003 * cos(q3) * sin(q4)) +
				cos(q2) * (-0.03714 * cos(q4) * sin(q3) + 0.00003 * sin(q3) * sin(q4))), 0, 0,

				0, 0.27613 * cos(q2) * cos(q3) - 0.00003 * cos(q2) * cos(q4) * sin(q3) -
				sin(q2) * (0.455 + 0.27613 * sin(q3) +
				cos(q3) * (0.00003 * cos(q4) + 0.03714 * sin(q4))) -
				0.03714 * cos(q2) * sin(q3) * sin(q4), -0.00003 * cos(q3) * cos(q4) * sin(q2) -
				0.27613 * sin(q2) * sin(q3) +
				cos(q2) * (0.27613 * cos(q3) -
				sin(q3) * (0.00003 * cos(q4) + 0.03714 * sin(q4))) -
				0.03714 * cos(q3) * sin(q2) * sin(q4), -0.03714 * cos(q4) * sin(q2) * sin(q3) +
				cos(q2) * cos(q3) * (0.03714 * cos(q4) - 0.00003 * sin(q4)) +
				0.00003 * sin(q2) * sin(q3) * sin(q4), 0, 0;




	Jv3 <<	-0.0181 * cos(q1) - 
				sin(q1) * (sin(q2) * (-0.455 - 0.00008 * cos(q3) - 0.04121 * sin(q3)) + 
				cos(q2) * (0.04121 * cos(q3) - 0.00008 * sin(q3))), 
				cos(q1) * (cos(q2) * (-0.455 - 0.00008 * cos(q3) - 0.04121 * sin(q3)) - 
				sin(q2) * (0.04121 * cos(q3) - 0.00008 * sin(q3))), 
				cos(q1) * (cos(q2) * (-0.00008 * cos(q3) - 0.04121 * sin(q3)) + 
				sin(q2) * (-0.04121 * cos(q3) + 0.00008 * sin(q3))), 0, 0, 0,

				-0.0181 * sin(q1) +
				cos(q1) * (sin(q2) * (-0.455 - 0.00008 * cos(q3) - 0.04121 * sin(q3)) +
				cos(q2) * (0.04121 * cos(q3) - 0.00008 * sin(q3))),
				sin(q1) * (cos(q2) * (-0.455 - 0.00008 * cos(q3) - 0.04121 * sin(q3)) -
				sin(q2) * (0.04121 * cos(q3) - 0.00008 * sin(q3))),
				sin(q1) * (cos(q2) * (-0.00008 * cos(q3) - 0.04121 * sin(q3)) +
				sin(q2) * (-0.04121 * cos(q3) + 0.00008 * sin(q3))), 0, 0, 0,

				0, 0.04121 * cos(q2) * cos(q3) -
				sin(q2) * (0.455 + 0.00008 * cos(q3) + 0.04121 * sin(q3)) -
				0.00008 * cos(q2) * sin(q3), -0.00008 * cos(q3) * sin(q2) +
				cos(q2) * (0.04121 * cos(q3) - 0.00008 * sin(q3)) -
				0.04121 * sin(q2) * sin(q3), 0, 0, 0;


;



	Jv2 <<	-0.11248 * cos(q1) - sin(q1) * (-0.00001 * cos(q2) - 0.19408 * sin(q2)), 
				cos(q1) * (-0.19408 * cos(q2) + 0.00001 * sin(q2)), 0, 0, 0, 0,
				
				-0.11248 * sin(q1) + cos(q1) * (-0.00001 * cos(q2) - 0.19408 * sin(q2)), 
				sin(q1) * (-0.19408 * cos(q2) + 0.00001 * sin(q2)), 0, 0, 0, 0,

				0, -0.00001 * cos(q2) - 0.19408 * sin(q2), 0, 0, 0, 0;

;
	Jv1 << -0.02043 * cos(q1) - 0.00002 * sin(q1), 0, 0, 0, 0, 0, 0.00002 * cos(q1) - 0.02043 * sin(q1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

	




	//计算角速度雅可比矩阵
	MatrixXd Jw6(3, 6);
	MatrixXd Jw5(3, 6);
	MatrixXd Jw4(3, 6);
	MatrixXd Jw3(3, 6);
	MatrixXd Jw2(3, 6);
	MatrixXd Jw1(3, 6);
	Jw6 <<	0, sin(q1), sin(q1), cos(q1) * cos(q2) * cos(q3) - cos(q1) * sin(q2) * sin(q3), -cos(q4) * sin(q1) + (-cos(q1) * cos(q3) * sin(q2) - 
		cos(q1) * cos(q2) * sin(q3)) * sin(q4), cos(q5) * (cos(q1) * cos(q2) * cos(q3) - cos(q1) * sin(q2) *sin(q3)) - (cos(q4) * (-cos(q1) * cos(q3) * sin(q2) - cos(q1) * cos(q2) * sin(q3)) + sin(q1) * sin(q4)) * sin(q5),

				0, -cos(q1), -cos(q1), cos(q2) * cos(q3) * sin(q1) - sin(q1) * sin(q2) * sin(q3), cos(q1) * cos(q4) + (-cos(q3) * sin(q1) * sin(q2) - cos(q2) * sin(q1) * sin(q3)) * sin(q4), cos(q5) * (cos(q2) * cos(q3) * sin(q1) - sin(q1) * sin(q2) * sin(q3)) - (cos(q4) * (-cos(q3) * sin(q1) * sin(q2) - cos(q2) * sin(q1) * sin(q3)) - cos(q1) * sin(q4)) * sin(q5),

				1, 0, 0, cos(q3) * sin(q2) + cos(q2) * sin(q3), (cos(q2) * cos(q3) - sin(q2) * sin(q3)) * sin(q4), cos(q5) * (cos(q3) * sin(q2) + cos(q2) * sin(q3)) - cos(q4) * (cos(q2) * cos(q3) - sin(q2) * sin(q3)) * sin(q5);

	// 下面的Jw[0]~Jw[4]只需要将Jw[5]一列一列的替换成0即可
	Jw5 <<		0, sin(q1), sin(q1), cos(q1)* cos(q2)* cos(q3) - cos(q1) * sin(q2) * sin(q3), -cos(q4) * sin(q1) + (-cos(q1) * cos(q3) * sin(q2) - cos(q1) * cos(q2) * sin(q3)) * sin(q4),0,

				0, -cos(q1), -cos(q1), cos(q2)* cos(q3)* sin(q1) - sin(q1) * sin(q2) * sin(q3), cos(q1)* cos(q4) + (-cos(q3) * sin(q1) * sin(q2) - cos(q2) * sin(q1) * sin(q3)) * sin(q4), 0,

				1, 0, 0, cos(q3)* sin(q2) + cos(q2) * sin(q3), (cos(q2) * cos(q3) - sin(q2) * sin(q3))* sin(q4),0;


	Jw4 <<		0, sin(q1), sin(q1), cos(q1)* cos(q2)* cos(q3) - cos(q1) * sin(q2) * sin(q3), 0, 0,

				0, -cos(q1), -cos(q1), cos(q2)* cos(q3)* sin(q1) - sin(q1) * sin(q2) * sin(q3), 0, 0,

				1, 0, 0, cos(q3)* sin(q2) + cos(q2) * sin(q3), 0, 0;


	Jw3 <<		0, sin(q1), sin(q1), 0, 0, 0,

				0, -cos(q1), -cos(q1), 0, 0, 0,

				1, 0, 0, 0, 0, 0;
	Jw2 <<		0, sin(q1), 0, 0, 0, 0,

				0, -cos(q1), 0, 0, 0, 0,

				1, 0, 0, 0, 0, 0;

	Jw1 <<		0, 0, 0, 0, 0, 0,

				0, 0, 0, 0, 0, 0,

				1, 0, 0, 0, 0, 0;


	//质量
	double m[6];
	m[0] = 8.30269;
	m[1] = 3.87623;
	m[2] = 5.31663;
	m[3] = 2.58124;
	m[4] = 2.66205;
	m[5] = 0.60427;




	//惯性张量
	Matrix3d I[6];
	/*
	I[0] <<     0.046167,		-2.67105e-05,	2.11669e-06,
				-2.67105e-05,	0.0396526,		-0.00879709,
				2.11669e-06,	-0.00879709,    0.0219442;
	I[1] <<		0.134685,		-4.95782e-06,   9.6385e-07,
				-4.95782e-06,   0.00613126,		-0.00239235,
				9.6385e-07,		-0.00239235,    0.135087;
	I[2] <<		0.0214709,		-1.09143e-06,   1.0309e-05,
				-1.09143e-06,   0.0106481,		0.00405475,
				1.0309e-05,		0.00405475,     0.017584;
	I[3] <<		0.0495614,		-6.49844e-06,	1.34622e-05,
				-6.49844e-06,   0.0454413,		0.0120927,
				1.34622e-05,    0.0120927,		0.0074644;
	I[4] <<		0.00795682,		3.09477e-06,	1.80346e-06,
				3.09477e-06,	0.00367246,		0.00149227,
				1.80346e-06,	0.00149227,		0.00627597;
	I[5] <<		0.000537685,	-2.769e-07,		-3.4558e-07,
				-2.769e-07,		0.000483409,	7.01904e-06,
				-3.4558e-07,	7.01904e-06,	0.000585171;
	*/
	I[0] << 0.0464076,		-2.67139e-05,	2.08863e-06,
			-2.67139e-05,	0.0398898,		-0.00882576,
			2.08863e-06,	-0.00882576,	0.0219477;

	I[1] << 0.13488,		-4.9503e-06,	-9.6821e-07,
			-4.9503e-06,   0.0061803,		0.00247697,
			-9.6821e-07,	0.00247697,		0.135233;

	I[2] << 0.0214817,		-1.0739e-06,	1.03167e-05,
			-1.0739e-06,	0.0106498,		0.00405079,
			1.03167e-05,	0.00405079,		0.017593;

	I[3] << 0.0497617,		-6.49556e-06,	1.34408e-05,
			-6.49556e-06,	0.0456381,		0.0121191,
			1.34408e-05,    0.0121191,		0.00746796;

	I[4] << 0.0497617,		-6.49556e-06,  1.34408e-05,
			-6.49556e-06,	0.0456381,		0.0121191,
			1.34408e-05,    0.0121191,		0.00746796;

	I[5] << 0.000538386,	-2.76881e-07,	3.46197e-07,
			-2.76881e-07,	0.00048411,		-7.04084e-06,
			3.46197e-07,	-7.04084e-06,  0.000585172;


	
	
	
	
	
	//质量矩阵
	MatrixXd MassMatrix(6,6);
	MassMatrix=Jv1.transpose()*m[0]*Jv1+Jw1.transpose()*I[0]*Jw1+
			   Jv2.transpose()*m[1]*Jv2+Jw2.transpose()*I[1]*Jw2+
			   Jv3.transpose()*m[2]*Jv3+Jw3.transpose()*I[2]*Jw3+
			   Jv4.transpose()*m[3]*Jv4+Jw4.transpose()*I[3]*Jw4+
			   Jv5.transpose()*m[4]*Jv5+Jw5.transpose()*I[4]*Jw5+
			   Jv6.transpose()*m[5]*Jv6+Jw6.transpose()*I[5]*Jw6;




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
VectorXd getTorque(const VectorXd q,const VectorXd q_dot, const VectorXd q_dot_dot,const MatrixXd M_past, const MatrixXd M_now,const double T) {
	//质量
	double m[6];
	m[0] = 8.30269;
	m[1] = 3.87623;
	m[2] = 5.31663;
	m[3] = 2.58124;
	m[4] = 2.66205;
	m[5] = 0.60427;
	//重力
	MatrixXd g(1, 4);
	g<< 0, 0, 9.8, 0;
	//重心位置
	Vector4d r[6];
	r[0] << 0.00002,	-0.051,		-0.02043,	1;
	r[1] << -0.26092,	0.00001,	-0.11248 ,	1;
	r[2] << 0.00008,	-0.0181,	0.04121,	1;
	r[3] << 0.00003,	-0.21887 ,	0.03714 ,	1;
	r[4] << -0.0002,	-0.01507 ,  0.03659,	1;
	r[5] << -0.00003,	-0.00106,	-0.03404 ,	1;

	double q1 = q(0);
	double q2 = q(1);
	double q3 = q(2);
	double q4 = q(3);
	double q5 = q(4);
	double q6 = q(5);
	

	//DH参数表
	MatrixXd DH_Table(6, 4);
	//         offst	d		a		alpha
	DH_Table << 0,		0.22,	0,		PI / 2,
				PI / 2, 0,		0.455,	0,
				0,		0,		0,		PI / 2,
				0,		0.495,	0,		PI / 2,
				0,		0,		0,		-PI / 2,
				0,		0.1565, 0,		0;


	Matrix4d T01 = DH2Trans(DH_Table(0, 0) + q1, DH_Table(0, 1), DH_Table(0, 2), DH_Table(0, 3));
	Matrix4d T12 = DH2Trans(DH_Table(1, 0) + q2, DH_Table(1, 1), DH_Table(1, 2), DH_Table(1, 3));
	Matrix4d T23 = DH2Trans(DH_Table(2, 0) + q3, DH_Table(2, 1), DH_Table(2, 2), DH_Table(2, 3));
	Matrix4d T34 = DH2Trans(DH_Table(3, 0) + q4, DH_Table(3, 1), DH_Table(3, 2), DH_Table(3, 3));
	Matrix4d T45 = DH2Trans(DH_Table(4, 0) + q5, DH_Table(4, 1), DH_Table(4, 2), DH_Table(4, 3));
	Matrix4d T56 = DH2Trans(DH_Table(5, 0) + q6, DH_Table(5, 1), DH_Table(5, 2), DH_Table(5, 3));
	Matrix4d T02 = T01 * T12;
	Matrix4d T03 = T02 * T23;
	Matrix4d T04 = T03 * T34;
	Matrix4d T05 = T04 * T45;
	Matrix4d T06 = T05 * T56;
	//此处获得T01~T56各个齐次矩阵的导数，用于后续计算
	Matrix4d dTdq[6];
	dTdq[0] << -sin(q1), 0, cos(q1), 0, cos(q1), 0, sin(q1), 0, 0, 0, 0, 0, 0, 0, 0, 0;
	dTdq[1] << -cos(q2), sin(q2), 0, DH_Table(1, 2)* cos(q2), -sin(q2), -cos(q2), 0, DH_Table(1, 2)* sin(q2), 0, 0, 0, 0, 0, 0, 0, 0;
	dTdq[2] << -sin(q3), 0, cos(q3), 0, cos(q3), 0, sin(q3), 0, 0, 0, 0, 0, 0, 0, 0, 0;
	dTdq[2] << -sin(q4), 0, cos(q4), 0, cos(q4), 0, sin(q4), 0, 0, 0, 0, 0, 0, 0, 0, 0;
	dTdq[4] << -sin(q5), 0, -cos(q5), 0, cos(q5), 0, -sin(q5), 0, 0, 0, 0, 0, 0, 0, 0, 0;
	dTdq[5] << -sin(q6), -cos(q6), 0, 0, cos(q6), -sin(q6), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
	//计算拉格朗日方程中对q求偏导的项，通过数值方法求得
	double delta = 0.000000001;//数值微分的步长
	//每个分量的小变化
	MatrixXd q_ = MatrixXd::Identity(6, 6);
	//计算每个分量的小变化对应的质量矩阵
	MatrixXd dMdq0(6, 6);
	MatrixXd dMdq1(6, 6);
	MatrixXd dMdq2(6, 6);
	MatrixXd dMdq3(6, 6);
	MatrixXd dMdq4(6, 6);
	MatrixXd dMdq5(6, 6);
	dMdq0 = getMassMatrix(q + delta * q_.col(0));
	dMdq1 = getMassMatrix(q + delta * q_.col(1));
	dMdq2 = getMassMatrix(q + delta * q_.col(2));
	dMdq3 = getMassMatrix(q + delta * q_.col(3));
	dMdq4 = getMassMatrix(q + delta * q_.col(4));
	dMdq5 = getMassMatrix(q + delta * q_.col(5));
	

	VectorXd dKdq(6);
	dKdq(0) = 0.5 * q_dot.transpose() * (dMdq0 - M_now) / delta * q_dot;
	dKdq(1) = 0.5 * q_dot.transpose() * (dMdq1 - M_now) / delta * q_dot;
	dKdq(2) = 0.5 * q_dot.transpose() * (dMdq2 - M_now) / delta * q_dot;
	dKdq(3) = 0.5 * q_dot.transpose() * (dMdq3 - M_now) / delta * q_dot;
	dKdq(4) = 0.5 * q_dot.transpose() * (dMdq4 - M_now) / delta * q_dot;
	dKdq(5) = 0.5 * q_dot.transpose() * (dMdq5 - M_now) / delta * q_dot;


	

	//计算重力矩
	VectorXd G(6);
	MatrixXd G0 =	m[0] * g * dTdq[0] * r[0] + \
					m[1] * g * dTdq[0] * (inverseHomogeneousTransform(T01) * T02) * r[1] + \
					m[2] * g * dTdq[0] * (inverseHomogeneousTransform(T01) * T03) * r[2] + \
					m[3] * g * dTdq[0] * (inverseHomogeneousTransform(T01) * T04) * r[3] + \
					m[4] * g * dTdq[0] * (inverseHomogeneousTransform(T01) * T05) * r[4] + \
					m[5] * g * dTdq[0] * (inverseHomogeneousTransform(T01) * T06) * r[5];
	G(0) = G0(0,0);
	MatrixXd G1 =	m[1] * g * T01 * dTdq[1] * r[1] + \
					m[2] * g * T01 * dTdq[1] * T23 * r[2] + \
					m[3] * g * T01 * dTdq[1] * (inverseHomogeneousTransform(T02) * T04) * r[3] + \
					m[4] * g * T01 * dTdq[1] * (inverseHomogeneousTransform(T02) * T05) * r[4] + \
					m[5] * g * T01 * dTdq[1] * (inverseHomogeneousTransform(T02) * T06) * r[5];
	G(1) = G1(0, 0);
	MatrixXd G2 =	m[2] * g * T02 * dTdq[2] * r[2] + \
					m[3] * g * T02 * dTdq[2] * (inverseHomogeneousTransform(T03) * T04) * r[3] + \
					m[4] * g * T02 * dTdq[2] * (inverseHomogeneousTransform(T03) * T05) * r[4] + \
					m[5] * g * T02 * dTdq[2] * (inverseHomogeneousTransform(T03) * T06) * r[5];
	G(2) = G2(0, 0);
	MatrixXd G3 =	m[3] * g * T03 * dTdq[3] * r[3] + \
					m[4] * g * T03 * dTdq[3] * (inverseHomogeneousTransform(T04) * T05) * r[4] + \
					m[5] * g * T03 * dTdq[3] * (inverseHomogeneousTransform(T04) * T06) * r[5];
	G(3) = G3(0, 0);
	MatrixXd G4 =	m[4] * g * T04 * dTdq[4] * r[4] + \
					m[5] * g * T04 * dTdq[4] * T56 * r[5];
	G(4) = G4(0,0);
	MatrixXd G5 =	m[5] * g * T05 * dTdq[5] * r[5];//OK
	G(5) = G5(0,0);
	/*
	G(1) = 37.9871 * (-0.00001 * cos(q2) - 0.19408 * sin(q2)) +
		52.103 * (0.04121 * cos(q2) * cos(q3) -
			sin(q2) * (0.455 + 0.00008 * cos(q3) + 0.04121 * sin(q3)) -
			0.00008 * cos(q2) * sin(q3)) +
		25.2962 * (0.27613 * cos(q2) * cos(q3) - 0.00003 * cos(q2) * cos(q4) * sin(q3) -
			sin(q2) * (0.455 + 0.27613 * sin(q3) +
				cos(q3) * (0.00003 * cos(q4) + 0.03714 * sin(q4))) -
			0.03714 * cos(q2) * sin(q3) * sin(q4)) +
		26.0881 * (0.495 * cos(q2 + q3) + 0.03659 * cos(q2 + q3) * cos(q5) -
			0.01507 * cos(q2) * sin(q3) * sin(q4) +
			cos(q2) * cos(q4) * sin(q3) * (0.0002 * cos(q5) + 0.03659 * sin(q5)) -
			0.0002 * cos(q2 + q3) * sin(q5) -
			sin(q2) * (0.455 +
				cos(q3) * (-0.0002 * cos(q4) * cos(q5) + 0.01507 * sin(q4) -
					0.03659 * cos(q4) * sin(q5)))) +
		5.92185 * (0.00003 * cos(q2) * cos(q4) * cos(q5) * cos(q6) * sin(q3) -
			0.00106 * cos(q2) * cos(q6) * sin(q3) * sin(q4) +
			0.12246 * cos(q2) * cos(q4) * sin(q3) * sin(q5) -
			0.00106 * cos(q2) * cos(q4) * cos(q5) * sin(q3) * sin(q6) -
			0.00003 * cos(q2) * sin(q3) * sin(q4) * sin(q6) +
			cos(q2) * cos(
				q3) * (0.495 + 0.12246 * cos(q5) - 0.00003 * cos(q6) * sin(q5) +
					0.00106 * sin(q5) * sin(q6)) -
			sin(q2) * (0.455 +
				sin(q3) * (0.495 + 0.12246 * cos(q5) - 0.00003 * cos(q6) * sin(q5) +
					0.00106 * sin(q5) * sin(q6)) +
				cos(q3) * (sin(q4) * (0.00106 * cos(q6) + 0.00003 * sin(q6)) +
					cos(q4) * (-0.00003 * cos(q5) * cos(q6) - 0.12246 * sin(q5) +
						0.00106 * cos(q5) * sin(q6)))));





	G(2) = 52.103 * (-0.00008 * cos(q3) * sin(q2) + cos(q2) * (0.04121 * cos(q3) - 0.00008 * sin(q3)) - 0.04121 * sin(q2) * sin(q3)) + 25.2962 * (-0.00003 * cos(q3) * cos(q4) * sin(q2) - 0.27613 * sin(q2) * sin(q3) + cos(q2) * (0.27613 * cos(q3) - sin(q3) * (0.00003 * cos(q4) + 0.03714 * sin(q4))) - 0.03714 * cos(q3) * sin(q2) * sin(q4)) + 26.0881 * (0.495 * cos(q2 + q3) + 0.03659 * cos(q2 + q3) * cos(q5) - 0.01507 * cos(q3) * sin(q2) * sin(q4) + cos(q3) * cos(q4) * sin(q2) * (0.0002 * cos(q5) + 0.03659 * sin(q5)) - 0.0002 * cos(q2 + q3) * sin(q5) - cos(q2) * sin(q3) * (-0.0002 * cos(q4) * cos(q5) + 0.01507 * sin(q4) - 0.03659 * cos(q4) * sin(q5))) + 5.92185 * (0.00003 * cos(q3) * cos(q4) * cos(q5) * cos(q6) * sin(q2) - 0.00106 * cos(q3) * cos(q6) * sin(q2) * sin(q4) + 0.12246 * cos(q3) * cos(q4) * sin(q2) * sin(q5) - 0.00106 * cos(q3) * cos(q4) * cos(q5) * sin(q2) * sin(q6) - 0.00003 * cos(q3) * sin(q2) * sin(q4) * sin(q6) - sin(q2) * sin(q3) * (0.495 + 0.12246 * cos(q5) - 0.00003 * cos(q6) * sin(q5) + 0.00106 * sin(q5) * sin(q6)) + cos(q2) * (cos(q3) * (0.495 + 0.12246 * cos(q5) - 0.00003 * cos(q6) * sin(q5) + 0.00106 * sin(q5) * sin(q6)) - sin(q3) * (sin(q4) * (0.00106 * cos(q6) + 0.00003 * sin(q6)) + cos(q4) * (-0.00003 * cos(q5) * cos(q6) - 0.12246 * sin(q5) + 0.00106 * cos(q5) * sin(q6)))));



	G(3) = 25.2962 * (-0.03714 * cos(q4) * sin(q2) * sin(q3) + cos(q2) * cos(q3) * (0.03714 * cos(q4) - 0.00003 * sin(q4)) + 0.00003 * sin(q2) * sin(q3) * sin(q4)) + 26.0881 * (-0.01507 * cos(q4) * sin(q2) * sin(q3) - sin(q2) * sin(q3) * sin(q4) * (0.0002 * cos(q5) + 0.03659 * sin(q5)) + cos(q2) * cos(q3) * (0.01507 * cos(q4) + 0.0002 * cos(q5) * sin(q4) + 0.03659 * sin(q4) * sin(q5))) + 5.92185 * (-0.00106 * cos(q4) * cos(q6) * sin(q2) * sin(q3) - 0.00003 * cos(q5) * cos(q6) * sin(q2) * sin(q3) * sin(q4) - 0.12246 * sin(q2) * sin(q3) * sin(q4) * sin(q5) - 0.00003 * cos(q4) * sin(q2) * sin(q3) * sin(q6) + 0.00106 * cos(q5) * sin(q2) * sin(q3) * sin(q4) * sin(q6) + cos(q2) * cos(q3) * (cos(q4) * (0.00106 * cos(q6) + 0.00003 * sin(q6)) - sin(q4) * (-0.00003 * cos(q5) * cos(q6) - 0.12246 * sin(q5) + 0.00106 * cos(q5) * sin(q6))));



	G(4) = 26.0881 * (-0.0002 * cos(q5) * sin(q2 + q3) + cos(q4) * sin(q2) * sin(q3) * (0.03659 * cos(q5) - 0.0002 * sin(q5)) - 0.03659 * sin(q2 + q3) * sin(q5) + cos(q2) * cos(q3) * (-0.03659 * cos(q4) * cos(q5) + 0.0002 * cos(q4) * sin(q5))) + 5.92185 * (0.12246 * cos(q4) * cos(q5) * sin(q2) * sin(q3) - 0.00003 * cos(q4) * cos(q6) * sin(q2) * sin(q3) * sin(q5) + 0.00106 * cos(q4) * sin(q2) * sin(q3) * sin(q5) * sin(q6) + cos(q3) * sin(q2) * (-0.00003 * cos(q5) * cos(q6) - 0.12246 * sin(q5) + 0.00106 * cos(q5) * sin(q6)) + cos(q2) * (sin(q3) * (-0.00003 * cos(q5) * cos(q6) - 0.12246 * sin(q5) + 0.00106 * cos(q5) * sin(q6)) + cos(q3) * cos(q4) * (-0.12246 * cos(q5) + 0.00003 * cos(q6) * sin(q5) - 0.00106 * sin(q5) * sin(q6))));



	G(5) = 5.92185 *(-0.00106 *cos(q4) *cos(q5) *cos(q6) *sin(q2) *sin(q3) -
		0.00003 *cos(q6)* sin(q2) *sin(q3) *sin(q4) -
		0.00003 *cos(q4) *cos(q5)* sin(q2)* sin(q3)* sin(q6) +
		0.00106* sin(q2)* sin(q3)* sin(q4)* sin(q6) +
		cos(q3) *sin(q2)*(0.00106* cos(q6)* sin(q5) + 0.00003* sin(q5)* sin(q6)) +
				cos(q2)*(sin(q3)*(0.00106 *cos(q6) *sin(q5) + 0.00003* sin(q5)* sin(q6)) +
						cos(q3)*(sin(q4)*(0.00003 *cos(q6) - 0.00106* sin(q6)) +
							cos(q4)*(0.00106 *cos(q5) *cos(q6) +
								0.00003 *cos(q5)* sin(q6)))));
*/



	//将上面计算结果拼接成最终输出的理论力矩
	VectorXd torque(6);

	torque = M_now * q_dot_dot + (M_now - M_past) / T * q_dot - dKdq + G;




	return torque;
}




int main()
{

	
	
	
	MatrixXd M_past(6, 6);
	MatrixXd M_now(6, 6);
	VectorXd q_now(6);
	q_now << 0.004, 0.004, 0.004, 0.004, 0.004, 0.004;

	VectorXd q_past(6);
	q_past << 0, 0, 0, 0, 0, 0;


	VectorXd v(6);
	v << 1, 1, 1, 1, 1, 1;

	VectorXd a(6);
	a << 0,      0,      0    ,  0 ,     0     , 0;
	
	M_past = getMassMatrix(q_past);
	M_now = getMassMatrix(q_now);

	VectorXd torque(6);
	torque = getTorque(q_now, v, a, M_past, M_now, 0.004);
	
	cout << "torque: " << endl << torque << endl;

	
	


	
	return 0;
}
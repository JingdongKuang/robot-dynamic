// vsCollisionDetection.h: 标准系统包含文件的包含文件
// 或项目特定的包含文件。

#pragma once

#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Geometry>
#define PI 3.1415926535897931
#define toRad 0.01745329252
using namespace Eigen;

class RobotDynamics {
public:
	RobotDynamics() {//初始化
		//质量
		m[0] = 8.30269;
		m[1] = 3.87623;
		m[2] = 5.31663;
		m[3] = 2.58124;
		m[4] = 2.66205;
		m[5] = 0.60427;

		//重力
		g.resize(1, 4);
		g << 0, 0, 9.8, 0;
		//重力矩
		G.resize(6, 1);
		G << 0, 0, 0, 0, 0, 0;
		//重心位置
		r[0] << 0.00002, -0.051, -0.02043, 1;
		r[1] << -0.26092, 0.00001, -0.11248, 1;
		r[2] << 0.00008, -0.0181, 0.04121, 1;
		r[3] << 0.00003, -0.21887, 0.03714, 1;
		r[4] << -0.0002, -0.01507, 0.03659, 1;
		r[5] << -0.00003, -0.00106, -0.03404, 1;

		//DH参数表      offst	d		a		alpha
		DH_Table.resize(6, 4);
		DH_Table << 0, 0.22, 0, PI / 2,
			PI / 2, 0, 0.455, 0,
			0, 0, 0, PI / 2,
			0, 0.495, 0, PI / 2,
			0, 0, 0, -PI / 2,
			0, 0.1565, 0, 0;

		//惯性张量
		I[0] << 0.0464076, -2.67139e-05, 2.08863e-06,
			-2.67139e-05, 0.0398898, -0.00882576,
			2.08863e-06, -0.00882576, 0.0219477;

		I[1] << 0.13488, -4.9503e-06, -9.6821e-07,
			-4.9503e-06, 0.0061803, 0.00247697,
			-9.6821e-07, 0.00247697, 0.135233;

		I[2] << 0.0214817, -1.0739e-06, 1.03167e-05,
			-1.0739e-06, 0.0106498, 0.00405079,
			1.03167e-05, 0.00405079, 0.017593;

		I[3] << 0.0497617, -6.49556e-06, 1.34408e-05,
			-6.49556e-06, 0.0456381, 0.0121191,
			1.34408e-05, 0.0121191, 0.00746796;

		I[4] << 0.0497617, -6.49556e-06, 1.34408e-05,
			-6.49556e-06, 0.0456381, 0.0121191,
			1.34408e-05, 0.0121191, 0.00746796;

		I[5] << 0.000538386, -2.76881e-07, 3.46197e-07,
			-2.76881e-07, 0.00048411, -7.04084e-06,
			3.46197e-07, -7.04084e-06, 0.000585172;

		//用于求导的矩阵
		Q.resize(4, 4);
		Q <<	0, -1, 0, 0,
				1,  0, 0, 0,
				0,  0, 0, 0,
				0,  0, 0, 0;
	}
	
	// 成员函数声明
	MatrixXd getMassMatrix(const VectorXd& q);//获取质量矩阵
	VectorXd getTorque(const VectorXd& q, const VectorXd& q_dot, const VectorXd& q_dot_dot, const MatrixXd& M_now);//获取关节力矩
	VectorXd getGravity();	//获取重力矩
	void getTransMatrix(const VectorXd &q);//获取齐次变换矩阵


private:
	/*变量声明*/
	//质量
	double m[6];
	//重力
	MatrixXd g;
	//重心位置
	Vector4d r[6];
	//DH参数表
	MatrixXd DH_Table;
	//惯性张量
	Matrix3d I[6];
	//重力矩
	VectorXd G;
	//用于求导的矩阵，属于机器人的特性dTdq=Q*T
	MatrixXd Q;
	//齐次变换矩阵
	Matrix4d T01;
	Matrix4d T12;
	Matrix4d T23;
	Matrix4d T34;
	Matrix4d T45;
	Matrix4d T56;
	Matrix4d T02;
	Matrix4d T03;
	Matrix4d T04;
	Matrix4d T05;
	Matrix4d T06;
	/*函数声明*/
	Matrix4d inverseHomogeneousTransform(const Matrix4d& transform);
	MatrixXd get_T_Derivative_of_time();//获取齐次变换矩阵对于时间的微分
	MatrixXd get_T_Derivative_of_qi(int T0i, int qi);//获取齐次变换矩阵对于qi的微分
};

﻿// vsCollisionDetection.h: 标准系统包含文件的包含文件
// 或项目特定的包含文件。

#pragma once

#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "Chebyshef_40hz_Order3.h"
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
		DH_Table <<		0,		0.22,	0,		PI / 2,
						PI / 2,	0,		0.455,	0,
						0,		0,		0,		PI / 2,
						0,		0.495,	0,		PI / 2,
						0,		0,		0,		-PI / 2,
						0,		0.1565, 0,		0;

		//惯性张量
		I[0] << 0.286765585945381,		2.33180108660000e-05,	-3.01797822000000e-05,
				2.33180108660000e-05,	0.276785769441076,		-0.0374635402223000,
				-3.01797822000000e-05,	-0.0374635402223000,	0.0254096500264570;

		I[1] << 0.329732556410064,		1.24808071840000e-05,	-3.39613350400000e-06,
				1.24808071840000e-05,	0.0551723550206150,		0.0870109096556320,
				-3.39613350400000e-06,	0.0870109096556320,		0.281093193834695;

		I[2] << 0.0322417088162830,		1.86192957840000e-05,	-2.61049976000000e-06,
				1.86192957840000e-05,	0.0123899195307320,		8.90741163699993e-05,
				-2.61049976000000e-06,	8.90741163699993e-05,	0.0266130666884150;

		I[3] << 0.249935684074060,		-3.62242239200000e-06,	-7.92051403600000e-06,
				-3.62242239200000e-06,	0.242255116698472,		0.0385644843765680,
				-7.92051403600000e-06,	0.0385644843765680,		0.0110249084018200;

		I[4] << 0.0121254081926500,		5.04285819000000e-06,	2.60580187000000e-06,
				5.04285819000000e-06,	0.00427702671386500,	2.43869688349998e-05,
				2.60580187000000e-06,	2.43869688349998e-05,	0.00983999463842500;

		I[5] << 0.00960026951610400,	2.57684214000000e-07,	1.87438712600000e-06,
				2.57684214000000e-07,	0.00954531577217500,	7.14197984520000e-05,
				1.87438712600000e-06,	7.14197984520000e-05,	0.000585850701615000;

		//用于求导的矩阵
		Q.resize(4, 4);
		Q <<	0, -1, 0, 0,
				1, 0, 0, 0,
				0, 0, 0, 0,
				0, 0, 0, 0;
	}

	// 成员函数声明
	MatrixXd getMassMatrix(const VectorXd& q);//获取质量矩阵
	VectorXd getTorque(const VectorXd& q, const VectorXd& q_dot, const VectorXd& q_dot_dot, const MatrixXd& M_now);//获取关节力矩
	VectorXd getGravity();	//获取重力矩
	void getTransMatrix(const VectorXd& q);//获取齐次变换矩阵
	void coutTransMatrix();//输出齐次变换矩阵
	VectorXd getTorque_Newton_Euler(const VectorXd& q, const VectorXd& q_dot, const VectorXd& q_dot_dot);//获取关节力矩牛顿欧拉法

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
	Matrix4d DH2Trans(const double theta, const double d, const double a, const double alpha);//DH参数转换为齐次变换矩阵
};
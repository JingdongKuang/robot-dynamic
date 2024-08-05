//指数积求解逆解

#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Geometry>
#define PI 3.1415926535897931
#define toRad 0.01745329252
using namespace Eigen;


class POE_kinematics {//基于指数积的运动学
public:
	POE_kinematics() {//初始化


	}
	//指数积求解正解
	MatrixXd forward_kinematics(VectorXd q);
	//指数积求解逆解
	MatrixXd inv_kinematics(MatrixXd T);
	//指数积求解雅克比矩阵
	MatrixXd Jacobian(MatrixXd T, VectorXd q);
	//通过转轴获得旋转矩阵
	Matrix3d w2R(Vector3d w, double theta);
//通过旋转矩阵获得转轴
	Vector4d R2w_and_theta(Matrix3d R);
private:
	//获得so3矩阵
	Matrix3d Operator_S(Vector3d w);

};
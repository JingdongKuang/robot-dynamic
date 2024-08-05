//ָ����������

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


class POE_kinematics {//����ָ�������˶�ѧ
public:
	POE_kinematics() {//��ʼ��


	}
	//ָ�����������
	MatrixXd forward_kinematics(VectorXd q);
	//ָ����������
	MatrixXd inv_kinematics(MatrixXd T);
	//ָ��������ſ˱Ⱦ���
	MatrixXd Jacobian(MatrixXd T, VectorXd q);
	//ͨ��ת������ת����
	Matrix3d w2R(Vector3d w, double theta);
//ͨ����ת������ת��
	Vector4d R2w_and_theta(Matrix3d R);
private:
	//���so3����
	Matrix3d Operator_S(Vector3d w);

};
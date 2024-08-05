
#include "mian_dynamics.h"
using namespace std;
int main()
{
	RobotDynamics objRobotDynamics;
	POE_kinematics objPOE_kinematics;
	ChebyshefFilter objChebyshefFilter;
	/***/
	VectorXd q_now(6);
	q_now << 3, 2, 2, 2, 2, 2;
	VectorXd v(6);
	v << 3, 2, 3, 3, 3, 3;
	VectorXd a(6);
	a << 3, 3, 4, 3, 3, 3;
	VectorXd id_para(36);
	id_para << 2.67879856071727,
		-2.58582663613946,
		-1.24808071877980e-05,
		0.126011316400632,
		5.32383350523857e-06,
		2.59235962858469,
		5.83200516840000,
		3.87622999990234e-05,
		2.40053551817202,
		1.64364357836122e-05,
		1.80074602398606e-05,
		8.90741163689975e-05,
		2.40729679557488,
		0.000425330399999091,
		-3.82639832350000,
		0.0119575940894542,
		9.37445760748612e-06,
		-3.48449540342605e-05,
		0.0385644843765673,
		0.0153019351156814,
		7.74372000001619e-05,
		-0.135984347100000,
		0.0173936972509635,
		-1.14668181056324e-06,
		1.00111813058912e-06,
		-2.43869688341877e-05,
		0.0193853104106020,
		5.32410000006756e-05,
		-0.171403313700000,
		5.49537439283182e-05,
		2.57684213678805e-07,
		-2.71500523660317e-07,
		-1.47844718488277e-05,
		0.000585850701614120,
		-1.81281000005556e-05,
		-0.000640526199999256;

	VectorXd torque(6);
	torque = objRobotDynamics.getTorque_Newton_Euler_MDH(q_now, v, a);
	MatrixXd Ytilde = objRobotDynamics.getYtilde(q_now, v, a);





	/**/
	int imax = 1000;
	MatrixXd q_qdot_qdotdot(6, 3);
	MatrixXd Ytilde_(6 * imax, 36);
	VectorXd torque_(6 * imax);
	for (int i = 0; i < imax; i++) {
		q_qdot_qdotdot = objRobotDynamics.getFourierTrajectory(i * 0.004);//MatrixXd::Random(6, 3);
		Ytilde_.block(6 * i, 0, 6, 36) = objRobotDynamics.getYtilde(q_qdot_qdotdot.col(0), q_qdot_qdotdot.col(1), q_qdot_qdotdot.col(2));
		torque_.segment(6 * i, 6) = objRobotDynamics.getTorque_Newton_Euler_MDH(q_qdot_qdotdot.col(0), q_qdot_qdotdot.col(1), q_qdot_qdotdot.col(2));
		//cout << "torque_" << i << ": " << endl << torque_.segment(6 * i, 6) << endl;
	}
	objRobotDynamics.identifyDynamicsParameters(Ytilde_, torque_);
	/**/
	VectorXd axis(3);
	axis << 0, 0, 1;
	double theta = 0.5;
	MatrixXd T = objPOE_kinematics.w2R(axis, theta);
	cout << "T: " << endl << T << endl;
	Vector4d w_and_theta = objPOE_kinematics.R2w_and_theta(T);
	cout << "w_and_theta: " << endl << w_and_theta << endl;
	return 0;
}
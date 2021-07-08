#include "tools.h"
#include <iostream>
#include <cmath>
using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::endl;
using std::cout;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth)
{
    /**
    * TODO: Calculate the RMSE here.
    */
    VectorXd rmse(4);
    rmse << 0,0,0,0;
    int n = estimations.size();
    // check the validity of the following inputs:
    //  * the estimation vector size should not be zero
    //  * the estimation vector size should equal ground truth vector size
    if (n == 0 || n != ground_truth.size())
    {
        cout << "Invalid estimation or ground_truth data" << endl;
        return rmse;
    }

    for (int i = 0; i < n; ++i)
    {
        VectorXd residual = estimations[i] - ground_truth[i];
        residual = residual.array()*residual.array();
        rmse += residual;
    }

    rmse = rmse/n;
    rmse = rmse.array().sqrt();
    return rmse;
}

VectorXd Tools::ConvertToPolarCoordinate(const VectorXd& x_state)
{
	VectorXd x_ = x_state;
	double px = x_state(0);
    double py = x_state(1);
    double vx = x_state(2);
    double vy = x_state(3);
	double pxpx = px * px;
	double pypy = py * py;
	double pxvx_pyvy = px * vx + py * vy; 
		
	double rho = sqrt(pxpx + pypy);
    double theta = atan2(py, px);
	double rho_dot = 0;

	if (fabs(rho) > 0.0001) 
	{
		rho_dot = pxvx_pyvy/rho;
	}
	VectorXd z_pred(3);
    z_pred << rho, theta, rho_dot;
	return z_pred;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state)
{
    /**
     * TODO:
     * Calculate a Jacobian here.
     */
    MatrixXd Hj(3,4);
    // recover state parameters
    double px = x_state(0);
    double py = x_state(1);
    double vx = x_state(2);
    double vy = x_state(3);

    // pre-compute
    double pxpx = px * px;
    double pypy = py * py;
    double pxpy = px * py;
    double c1 = pxpx + pypy;
    double c2 = sqrt(c1);
    double c3 = (c1 * c2);

    // check division by zero
    if (fabs(c1) < 0.0001)
    {
        cout << "CalculateJacobian () - Error - Division by Zero" << endl;
        return Hj;
    }

    double pxDc2 = px/c2;
    double pyDc2 = py/c2;


    // compute the Jacobian matrix
    Hj << pxDc2                ,  pyDc2                    , 0    , 0,
          -(py / c1)           , (px / c1)                 , 0    , 0,
    (pypy * vx - vy * pxpy)/c3 , (pxpx * vy - vx * pxpy)/c3, pxDc2, pyDc2;

    return Hj;
}


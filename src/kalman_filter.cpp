#include "kalman_filter.h"
#include "tools.h"
#include <cmath>
using Eigen::MatrixXd;
using Eigen::VectorXd;

/*
 * Please note that the Eigen library does not initialize
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in)
{
    x_ = x_in;
    P_ = P_in;
    F_ = F_in;
    H_ = H_in;
    R_ = R_in;
    Q_ = Q_in;
}

void KalmanFilter::Predict()
{
    /**
     * TODO: predict the state
     */
	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z)
{
    /**
     * TODO: update the state by using Kalman Filter equations
     */
    VectorXd z_pred = H_ * x_;
    VectorXd y = z - z_pred;

    //new estimate
    UpdateEstimate(y);
	
}

void KalmanFilter::UpdateEstimate(const VectorXd &y)
{
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd K = P_ * Ht * Si;
	
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z)
{
    /**
     * TODO: update the state by using Extended Kalman Filter equations
     */
	 
	Tools tools; 
	VectorXd z_pred = tools.ConvertToPolarCoordinate(x_);
	VectorXd y = z - z_pred;
	
	while(y(1)> M_PI)
	{
		y(1) -= M_PI;
	}
	while(y(1)< -M_PI)
	{
		y(1) += M_PI;
	}
	
    UpdateEstimate(y);

}

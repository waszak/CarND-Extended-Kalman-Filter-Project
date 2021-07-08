#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF()
{
    is_initialized_ = false;

    previous_timestamp_ = 0;
	previous_diff = 0;
    // initializing matrices
    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);
    H_laser_ = MatrixXd(2, 4);
    Hj_ = MatrixXd(3, 4);

    //measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
             0, 0.0225;

    //measurement covariance matrix - radar
    R_radar_ << 0.09, 0, 0,
             0, 0.0009, 0,
             0, 0, 0.09;

    /**
     * TODO: Finish initializing the FusionEKF.
     * TODO: Set the process and measurement noises
     */
    noise_ax = 9;
    noise_ay = 9;

    ekf_.H_ = MatrixXd(2, 4);
    ekf_.H_ << 1, 0, 0, 0,
               0, 1, 0, 0;
			
	H_laser_ = MatrixXd(2, 4);
    H_laser_<< 1, 0, 0, 0,
               0, 1, 0, 0;

    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ << 1, 0, 1, 0,
               0, 1, 0, 1,
               0, 0, 1, 0,
               0, 0, 0, 1;

    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0, 0,
               0, 1, 0, 0,
               0, 0, 1000, 0,
               0, 0, 0, 1000;

}

void FusionEKF::UpdateQ(long long diff, double dt)
{

	//Another idea would be calculate avg diff and cache some Q matrices
	//For simulator this is sufficient because diff is 50000
	if (diff == previous_diff)
	{
		return;
	}
	previous_diff = diff;
	
	double dt_2 = dt * dt; 
	double dt_3 = dt_2 * dt;
	double dt_4 = dt_3 * dt; 
	double dt_3_2 = dt_3 / 2;
	double ax_dt_3_2 = dt_3_2 * noise_ax;
	double ay_dt_3_2 = dt_3_2 * noise_ay;
	ekf_.Q_ = MatrixXd(4, 4);
	ekf_.Q_ << dt*ax_dt_3_2/2, 0          , ax_dt_3_2      , 0,
	         0            , dt*ay_dt_3_2/2, 0              , ay_dt_3_2,
	         ax_dt_3_2    , 0          , dt_2 * noise_ax, 0,
 	         0            , ay_dt_3_2  , 0              , dt_2 * noise_ay;
			 	 
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack)
{
    /**
     * Initialization
     */
    if (!is_initialized_)
    {
        /**
         * TODO: Initialize the state ekf_.x_ with the first measurement.
         * TODO: Create the covariance matrix.
         * You'll need to convert radar from polar to cartesian coordinates.
         */

        // first measurement
        cout << "EKF: " << endl;
        ekf_.x_ = VectorXd(4);
        ekf_.x_ << 1, 1, 1, 1;

        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
        {
            // TODO: Convert radar from polar to cartesian coordinates
            //         and initialize state.
			double rho = measurement_pack.raw_measurements_[0];
			double phi = measurement_pack.raw_measurements_[1];
			double x = rho * cos(phi);
			double y = rho * sin(phi);

			ekf_.x_ << x, y, 0, 0;
            
        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
        {
            // TODO: Initialize state.
            ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
        }

        // done initializing, no need to predict or update
        is_initialized_ = true;
		previous_timestamp_ = measurement_pack.timestamp_;
        return;
    }

    /**
     * Prediction
     */

    /**
     * TODO: Update the state transition matrix F according to the new elapsed time.
     * Time is measured in seconds.
     * TODO: Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
     */
    long long diff = measurement_pack.timestamp_ - previous_timestamp_;
    float dt = (diff) / 1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;

    // TODO: YOUR CODE HERE
    UpdateQ(diff, dt);

    // Modify the F matrix so that the time is integrated
	
    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;


    cout<<"Predict\n";

    ekf_.Predict();

    /**
     * Update
     */

    /**
     * TODO:
     * - Use the sensor type to perform the update step.
     * - Update the state and covariance matrices.
     */

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
    {
		cout<<"Radar"<<endl;
		Hj_ = tools.CalculateJacobian(ekf_.x_);
		ekf_.R_ = R_radar_;
		ekf_.H_ = Hj_;
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);

    }
    else
    {
		cout<<"Laser"<<endl;
        // TODO: Laser updates
		ekf_.R_ = R_laser_;
		ekf_.H_ = H_laser_;
        ekf_.Update(measurement_pack.raw_measurements_);
    }

    // print the output
	//cout<<"diff = "<<diff<<endl;
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}

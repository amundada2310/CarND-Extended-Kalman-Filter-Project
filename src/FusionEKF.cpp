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
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during initializing)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during initializing)
  use_radar_ = true;

  previous_timestamp_ = 0;

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

  //1. initialize variables and matrices (x, F, H_laser, H_jacobian, P, etc.)

   // initialize state vector - x
  ekf_.x_ = VectorXd(4);

  // initialize covariance matrix - P
  ekf_.P_ = MatrixXd(4, 4);

  ekf_.P_ << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1000, 0,
            0, 0, 0, 1000;

  //initialize transition matrix F_
  ekf_.F_ = MatrixXd(4, 4);

  ekf_.F_ << 1, 0, 1, 0,
      0, 1, 0, 1,
      0, 0, 1, 0,
      0, 0, 0, 1;

  //initialize process noise
  ekf_.Q_ = MatrixXd(4, 4);
      
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
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
      std::cout << "Kalman Filter Initialization " << std::endl;


    // first measurement
    cout << "EKF - First Measurement: " << endl;
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) 
    {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
        double rho = measurement_pack.raw_measurements_[0];//first radar measurement in rho- range
        double phi = measurement_pack.raw_measurements_[1];//second is the phi - bearing
        double rho_dot = measurement_pack.raw_measurements_[2];

        ekf_.x_ << rho * cos(phi), //calculate and input the value  Px
            rho* sin(phi),//calculate and input the value Py
            0, //we do not have details regarding this value vx
            0;//we do not have details regarding this value vy

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) 
    {
      // TODO: Initialize state.
        ekf_.x_ << measurement_pack.raw_measurements_[0],
                    measurement_pack.raw_measurements_[1],
                    0, //we do not have details regarding this value vx
                    0;//we do not have details regarding this value vy
    }

    else
    {
        std::cout << "NO MEASUREMENT RECORDED" << std::endl;// if nothing recorded from sensor
    }

    // put the current time stamp into previous 
    previous_timestamp_ = measurement_pack.timestamp_;//input the time value
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  // this is not a initial measurement -

  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;//calculate the time between the previous and current measurements
  previous_timestamp_ = measurement_pack.timestamp_;//input current time as the previous time 

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */


  //update the state transition matrix F

  ekf_.F_ << 1, 0, dt, 0,
        0, 1, 0, dt,
        0, 0, 1, 0,
        0, 0, 0, 1;

  //update the process noise matrix
  double dt_2 = dt * dt;
  double dt_3 = dt_2 * dt;
  double dt_4 = dt_3 * dt;

  double noise_ax = 9;
  double noise_ay = 9;

  ekf_.Q_ << (dt_4 / 4) * noise_ax, 0, (dt_3 / 2) * noise_ax, 0,
      0, (dt_4 / 4) * noise_ay, 0, (dt_3 / 2) * noise_ay,
      (dt_3 / 2) * noise_ax, 0, (dt_2* noise_ax), 0,
      0, (dt_3 / 2) * noise_ay, 0, (dt_2* noise_ay);


  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR && use_radar_) 
  {
    // TODO: Radar updates
      // we need to calculate jacobian as it is a nonlinear function
      Hj_ = tools.CalculateJacobian(ekf_.x_);
      ekf_.H_ = Hj_;// input new H matrix
      ekf_.R_ = R_radar_;// input the R matrix
      ekf_.UpdateEKF(measurement_pack.raw_measurements_);// call Update function of Radar
  } 
  
  if (measurement_pack.sensor_type_ == MeasurementPackage::LASER && use_laser_)
  {
    // TODO: Laser updates
      H_laser_ << 1, 0, 0, 0,
                  0, 1, 0, 0;

      ekf_.H_ = H_laser_; // input the H matrix values
      ekf_.R_ = R_laser_;// input the R matrix values
      ekf_.Update(measurement_pack.raw_measurements_); // call Update function of laser
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}

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

  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;  
  
  // the initial transition matrix F_
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
            0, 1, 0, 1,
            0, 0, 1, 0,
            0, 0, 0, 1;
  
    // state covariance matrix P_
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1000, 0,
             0, 0, 0, 1000;
//   process noise covariance matrix
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
  if (!is_initialized_) {

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);


    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
//       convert from polar to cartesian
      float ro = measurement_pack.raw_measurements_[0];
      float fi = measurement_pack.raw_measurements_[1];
//       float ro_dot = measurement_pack.raw_measurements_[2];
      ekf_.x_ <<ro*cos(fi), 
              ro*sin(fi), 
//               ro_dot*cos(fi), 
//               ro_dot*sin(fi); 
              0,0;

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      ekf_.x_ << measurement_pack.raw_measurements_[0], 
              measurement_pack.raw_measurements_[1], 
              0, 
              0;    
    }
    
    // done initializing, no need to predict or update
    previous_timestamp_ = measurement_pack.timestamp_;      
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  // compute the time elapsed between the current and previous measurements
  // dt - expressed in seconds
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
   
  // Modify the F matrix so that the time is integrated
  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;
  // Set the process covariance matrix Q
  float t2 = dt*dt;
  float t3 = t2*dt/2;
  float t4 = t2*t2/4;
  float noise_ax = 9;
  float noise_ay = 9;
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << t4*noise_ax, 0, t3*noise_ax, 0,
            0, t4*noise_ay, 0, t3*noise_ay,
            t3*noise_ax, 0, t2*noise_ax, 0,
            0, t3*noise_ay, 0, t2*noise_ay;

  ekf_.Predict();

  /**
   * Update
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
     ekf_.R_ = R_radar_;
//     calculate Jacobian for predicted state
     Hj_ = tools.CalculateJacobian(ekf_.x_);
     ekf_.H_ = Hj_;
     VectorXd z = VectorXd(3);
     z << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], measurement_pack.raw_measurements_[2];
     ekf_.UpdateEKF(z);

  } else {     
      ekf_.R_ = R_laser_;
      ekf_.H_ = H_laser_;
      VectorXd z = VectorXd(2);
      z << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1];
      ekf_.Update(z);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}

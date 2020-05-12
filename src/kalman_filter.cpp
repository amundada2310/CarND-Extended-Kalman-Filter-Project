#include "kalman_filter.h"
#include<math.h>//added
#include<fstream>//added

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
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
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;

    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;

    //Calculate NIS- Normalized Innovation squared

    //the NIS meaurement value initialize for LIdar
    double NIS_laser;
    VectorXd z2 = z - z_pred;
    NIS_laser = z2.transpose() * S.inverse() * z2;

    std::ofstream log("NIS_laser.txt", std::ios_base::app | std::ios_base::out);//creating a new text file called NIS_laser

    log << NIS_laser;//add NIS calculated value 
    log << "\n";

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */

    //New value that come in are in rho, phi and rhodot form
    // predicted value are in px, py, vx, vy form - cartisian coordinates
    // we need to first convert the previous value in rho, phi and rhodot form - polar form

    //load the values from x state matrix
    double px, py, vx, vy;
    px = x_[0];
    py = x_[1];
    vx = x_[2];
    vy = x_[3];

    // convert the values into polar coordinates
    double rho, phi, rho_dot;
    rho = sqrt((px * px) + (py * py));
    phi = atan2(py, px);  // returns values between -pi and pi

    //avoid division by zero '0'
    if (fabs(rho < 0.000001))
    {
      rho_dot = 0.0;
    }
    else
    {
      rho_dot = ((px * vx) + (py * vy)) / rho;
    }

    // input the polar state x predicted values
    VectorXd z_pred = VectorXd(3);
    z_pred << rho, phi, rho_dot;

    // create the measurement error
    VectorXd y = z - z_pred;

    //normalize the angles
    while (y(1) > M_PI)
    {
        y(1) -= (2. * M_PI);
    }

    while (y(1) < -M_PI)
    {
        y(1) -= (2. * M_PI);
    }

    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;

    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;

    //3. Calculate NIS normalized innovation squared to check consistency of the designed filter
    //the NIS meaurement value initialize for LIdar
    double NIS_radar;
    VectorXd z2 = z - z_pred;

    NIS_radar = z2.transpose() * S.inverse() * z2;

    std::ofstream log("NIS_radar.txt", std::ios_base::app | std::ios_base::out);//creating a new text file called NIS_radar

    log << NIS_radar;//add the NIS calculated values
    log << "\n";
}

#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.8;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.6;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  // Set state dimension (px, py, speed v (magnitude of the velocty), Si (angle of the orientation toward which the tracking object move
  // and yaw rate Si dot)
  n_x_ = 5;

  // Radar measurement dimension can measure r, phi, and r_dot
  n_z_ = 3;

  // Lidar can only measure px and py
  n_z_lidar = 2;

  // Set augmented dimension
  n_aug_ = 7;

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // initial state vector
  x_ = VectorXd(n_x_);

  // Augmented mean vector
  x_aug = VectorXd(n_aug_);

  // Sigma point matrix
  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // Sensor matrices
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << std_laspx_ * std_laspx_, 0,
          0, std_laspy_ * std_laspy_;

  R_radar_ = MatrixXd(3, 3);
  R_radar_ <<	std_radr_ * std_radr_, 0, 0,
          0, std_radphi_ * std_radphi_, 0,
          0, 0, std_radrd_ * std_radrd_;

  // Vector for weights
  weights = VectorXd(2 * n_aug_ + 1);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  /**
   * Initialization
   */

  if(!is_initialized_){
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
      double rho = meas_package.raw_measurements_[0];
      // phi is also yaw
      double phi = meas_package.raw_measurements_[1];
      double v = meas_package.raw_measurements_[2];

      double px = rho*cos(phi);
      double py = rho*sin(phi);
      double yawd = 0.0;

      x_ << px, py, v, phi, yawd;
    }

    else if(meas_package.sensor_type_ == MeasurementPackage::LASER){
      double px = meas_package.raw_measurements_[0];
      double py = meas_package.raw_measurements_[1];

      x_ << px, py, 0, 0, 0;
    }

    // Init covariance matrix
    P_ <<	1, 0, 0, 0, 0,
            0, 1, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;

    previous_timestamp_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */
  double dt = (meas_package.timestamp_ - previous_timestamp_)/1000000.0;

  Prediction(dt);
  previous_timestamp_ = meas_package.timestamp_;

  /**
   * Update
   */

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  }

  else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  // Generate augmented sigma points Xsig_aug_
  GenerateAugmentedSigmaPoints();

  // Predict the future Xsig_pred_
  PredictAugmentedSigmaPoints(delta_t);

  // Predict mean and covariance
  PredictMeanAndCovariance();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  Calculate the lidar NIS.
  */
  double px = meas_package.raw_measurements_[0];
  double py = meas_package.raw_measurements_[1];

  VectorXd z = VectorXd(n_z_lidar);
  z << px, py;

  MatrixXd Z_sig = MatrixXd(n_z_lidar, 2 * n_aug_ + 1);
  VectorXd z_pred = VectorXd(n_z_lidar);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      // measurement model
      Z_sig(0,i) = Xsig_pred_(0,i);
      Z_sig(1,i) = Xsig_pred_(1,i);
  }

  // mean predicted measurement
  z_pred.fill(0.0);
  for (int i=0; i < 2 * n_aug_+1; i++) {
    z_pred = z_pred + weights(i) * Z_sig.col(i);
  }

  // Covariance Matrix
  MatrixXd S = MatrixXd(n_z_lidar, n_z_lidar);

  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd z_diff = Z_sig.col(i) - z_pred;

    while (z_diff(1) > M_PI) z_diff(1)-=2. *M_PI;
    while (z_diff(1) < -M_PI) z_diff(1)+=2. *M_PI;

    S = S + weights(i) * z_diff * z_diff.transpose();
  }

  S = S + R_laser_;

  /*
   * update the state matrix x_ and the covariance matrix P_
   */

  // cross-correlation matrix Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_lidar);

  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd z_diff = Z_sig.col(i) - z_pred;

    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // normalize angle
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  VectorXd z_diff = z - z_pred;

  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  Calculate the radar NIS.
  */

  double rho = meas_package.raw_measurements_[0];
  double phi = meas_package.raw_measurements_[1];
  double v = meas_package.raw_measurements_[2];

  VectorXd z = VectorXd(n_z_);
  z << rho, phi, v;

  // predict radar measurement
  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);
  VectorXd z_pred = VectorXd(n_z_);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double vel  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v_x = cos(yaw) * vel;
    double v_y = sin(yaw) * vel;

    // rho
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);

    //phi
    Zsig(1,i) = atan2(p_y,p_x);

	Zsig(2, i) = (Zsig(0,i) < 0.0001 ? (p_x * v_x + p_y * v_y) / 0.0001 : (p_x * v_x + p_y * v_y) / Zsig(0, i));

  }

  //mean predicted measurement
  z_pred.fill(0.0);
  for (int i=0; i < 2 * n_aug_+1; i++) {
    z_pred = z_pred + weights(i) * Zsig.col(i);
  }

  //measurement covariance matrix
  MatrixXd S = MatrixXd(n_z_, n_z_);

  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred; //residual
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  S = S + R_radar_;

  // Cross correlation matrix Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);

  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // normalize angle
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // anormalize angle
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // esidual
  VectorXd z_diff = z - z_pred;

  // normalize angle
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

}

void UKF::GenerateAugmentedSigmaPoints() {
  // Augmented mean vector
  x_aug = VectorXd(n_aug_);

  // Augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  // Augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  // Augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5, 5) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;

  // Square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  // Augmented sigma points
  Xsig_aug_.col(0) = x_aug;
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug_.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug_.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }
}


void UKF::PredictAugmentedSigmaPoints(double dt){

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    double p_x = Xsig_aug_(0, i);
    double p_y = Xsig_aug_(1, i);
    double v = Xsig_aug_(2, i);
    double yaw = Xsig_aug_(3, i);
    double yawd = Xsig_aug_(4, i);
    double a_ = Xsig_aug_(5, i);
    double yawdd_n = Xsig_aug_(6, i);

    //predicted state values
    double pred_px, pred_py;

    // check division by zero
    if (fabs(yawd) > 0.0001) {
      pred_px = p_x + v / yawd * (sin(yaw + yawd * dt) - sin(yaw));
      pred_py = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * dt));
    } else {
      pred_px = p_x + v * dt * cos(yaw);
      pred_py = p_y + v * dt * sin(yaw);
    }

    double pred_v = v;
    double yaw_p = yaw + yawd * dt;
    double yawd_p = yawd;

    // noise
    pred_px = pred_px + 0.5 * a_ * dt * dt * cos(yaw);
    pred_py = pred_py + 0.5 * a_ * dt * dt * sin(yaw);
    pred_v = pred_v + a_ * dt;

    yaw_p = yaw_p + 0.5 * yawdd_n * dt * dt;
    yawd_p = yawd_p + yawdd_n * dt;

    Xsig_pred_(0, i) = pred_px;
    Xsig_pred_(1, i) = pred_py;
    Xsig_pred_(2, i) = pred_v;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
  }
}

void UKF::PredictMeanAndCovariance(){

  // weights
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights(0) = weight_0;

  for (int i=1; i < 2 * n_aug_ + 1; i++) {
    double weight = 0.5 / (n_aug_ + lambda_);
    weights(i) = weight;
  }

  // Predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    x_ = x_+ weights(i) * Xsig_pred_.col(i);
  }

  // Predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // normalize angle
    while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3) <- M_PI) x_diff(3) += 2.*M_PI;
    P_ = P_ + weights(i) * x_diff * x_diff.transpose();
  }
}

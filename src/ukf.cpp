#include <iostream>
#include "ukf.h"

using namespace std;
/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  
  is_initialized_ = false;
  n_x_ =5; 
  n_aug_ =7;
  n_z_ = 3;
  lambda_ = 3 - n_aug_; //n_aug_; 
 
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // predicted sigma points
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);
  
  //create matrix for sigma points in measurement space
  Zsig_ = MatrixXd(n_z_, 2 * n_aug_ + 1);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.1; //30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5; //30;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  //radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  //radar measurement noise standard deviation angle in rad
  std_radphi_ =  0.02;

  //radar measurement noise standard deviation radius change in m/s
  std_radrd_ =  0.3;
/**

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  time_us_ = 0;
  weights_ = VectorXd (2*n_aug_+1);
  //mean predicted measurement
  z_pred_ = VectorXd(n_z_);
  //measurement covariance matrix S
  S_ = MatrixXd(n_z_,n_z_);
  NIS_radar_ = 0;
  NIS_laser_ = 0;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
bool UKF::Initialization(MeasurementPackage meas_package) {
   float px = 0;
   float py = 0;
   float vel = 0;
   float yaw_angle = 0;
   float yaw_rate = 0;
    x_ << px, py, vel, yaw_angle, yaw_rate;
    Xsig_pred_.fill(0.0);
    z_pred_.fill(0.0);
    S_.fill(0.0);
    Zsig_.fill(0.0);
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
          float rho = meas_package.raw_measurements_[0];
          float phi = meas_package.raw_measurements_[1];
          float rho_dot = meas_package.raw_measurements_[2];
	  px = rho*cos(phi);
          py = rho*sin(phi);
          float vx = rho_dot * cos(phi);
          float vy = rho_dot * sin(phi);
	  vel = fabs(rho_dot);
	  yaw_angle = atan(phi);
          if(vx !=0){
                yaw_angle = vy/vx;
            }
	  x_ << px, py, vel, yaw_angle, 0;
          P_ <<  1, 0, 0, 0, 0,
                 0, 1, 0, 0, 0,
                 0, 0, 1, 0, 0,
                 0, 0, 0, 0.5, 0,
                 0, 0, 0, 0, 0.5;
        }
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER){
           px = meas_package.raw_measurements_[0];
           py = meas_package.raw_measurements_[1];
	   x_ << px, py, 0, 0, 0;
           P_ << 1, 0, 0, 0, 0,
                 0, 1, 0, 0, 0,
                 0, 0, 1, 0, 0,
                 0, 0, 0, 10, 0,
                 0, 0, 0, 0, 10;
         }
  time_us_ = meas_package.timestamp_;

  // done initializing, no need to predict or update
  if(px == 0 && py==0){
        return false;
  }else{
        return true;
        }
}


void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  if (!is_initialized_){
	if(Initialization(meas_package)){
                        is_initialized_ = true;
              }
         return;
  }
 
 
 float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
 time_us_ = meas_package.timestamp_; 
 Prediction(dt); 


  //perform update step
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
	PredictRadarMeasurement();
	//cout << "radar update" << endl;
	UpdateRadar(meas_package);

  }else{
	//cout << "lidar update" << endl;
	UpdateLidar(meas_package);
	}
  Normalization(x_(3));

 return;
}
/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  
  //create augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug(n_x_) = 0;
  x_aug(n_x_+1) = 0;
   
  
  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(n_x_,n_x_) = std_a_*std_a_;
  P_aug(n_x_+1,n_x_+1) = std_yawdd_*std_yawdd_;
  
  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();
  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
 
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }

 

  //create matrix with predicted sigma points as columns
  //MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
  
   //double delta_t = 0.1; //time diff in sec

 //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //angle normalization
    Normalization(yaw_p);
   
    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
  

  //create vector for weights
  VectorXd weights = VectorXd::Zero(2*n_aug_+1);

  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  //create covariance matrix for prediction
   MatrixXd P = MatrixXd(n_x_, n_x_);

  // set weights
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights(i) = weight;
  }
  //predicted state mean
  x.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
     x = x + weights(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    //angle normalization
    Normalization(x_diff(3));
    P = P + weights(i) * x_diff * x_diff.transpose() ;
  }

 x_ =x;
 P_ =P;
 Normalization(x_(3));
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  VectorXd z(2);
  z << meas_package.raw_measurements_;
  MatrixXd R(2, 2);
  R << std_laspx_*std_laspx_, 0,
		0, std_laspy_ * std_laspy_;

  MatrixXd H(2, 5);
  H<< 1, 0, 0, 0,0,
	0, 1, 0, 0,0;
  VectorXd z_pred = H * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H.transpose();
  MatrixXd S = H * P_ * Ht + R;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H) * P_;

  //update NIS
  NIS_laser_ = ComputeNIS(z_pred, S, z);

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */

void UKF::PredictRadarMeasurement() {


  //set vector for weights
  VectorXd weights = VectorXd(2*n_aug_+1);
   double weight_0 = lambda_/(lambda_+n_aug_);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  
    double weight = 0.5/(n_aug_+lambda_);
    weights(i) = weight;
  }
  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig_(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig_(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig_(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }
  
  z_pred_.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred_ = z_pred_ + weights(i) * Zsig_.col(i);
  }
 
  //measurement covariance matrix S
  //MatrixXd S_ = MatrixXd(n_z_,n_z_);
  S_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;

    //angle normalization
    Normalization(z_diff(1));

    S_ = S_ + weights(i) * z_diff * z_diff.transpose();
  }
  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z_,n_z_);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S_ = S_ + R;


}
void UKF::UpdateRadar(MeasurementPackage meas_package) {

  //set vector for weights
  VectorXd weights = VectorXd(2*n_aug_+1);
   double weight_0 = lambda_/(lambda_+n_aug_);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights(i) = weight;
  }


  //create example vector for incoming radar measurement
  const VectorXd &z = meas_package.raw_measurements_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);


  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;
    //angle normalization
    Normalization(z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    Normalization(x_diff(3));

    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
  }
  
  //Kalman gain K;
  MatrixXd Kn = Tc * S_.inverse();
  //residual
  VectorXd z_diff = z - z_pred_;
  //angle normalization
  Normalization(z_diff(1));


  NIS_radar_ = ComputeNIS(z_pred_, S_, z);
  //update state mean and covariance matrix
  x_ = x_ + Kn * z_diff;
  P_ = P_ - Kn*S_*Kn.transpose();

}
void UKF::Normalization(double &angle) {
        while (angle> M_PI) {
                angle-=2.*M_PI;
        }
        while (angle<-M_PI) {
                angle+=2.*M_PI;
        }

}
double UKF::ComputeNIS(const VectorXd &z_pred, const MatrixXd &S, const VectorXd &z) {
        double nis = 0;
        VectorXd y = z - z_pred;
        MatrixXd Si = S.inverse();
        int mea_size = z.size();

        Eigen::Map<MatrixXd> yt_matrix(y.data(), 1,mea_size);
        Eigen::Map<MatrixXd> y_matrix(y.data(), mea_size,1);

        MatrixXd temp = yt_matrix * Si * y_matrix;
        nis = temp(0,0);
        return nis;

}

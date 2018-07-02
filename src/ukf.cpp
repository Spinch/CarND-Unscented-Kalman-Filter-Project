
/** @file ukf.cpp
 *  @author Andrey N. Zabegaev <speench@gmail.com>
 */

#include <math.h>
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"
#include "ukf.h"

#include <iostream>
#define TRACE std::cout << __LINE__ << std::endl;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF()
{
    // if this is false, laser measurements will be ignored (except during init)
    _use_laser = true;

    // if this is false, radar measurements will be ignored (except during init)
    _use_radar = true;

    _std_a = 0.4;  // Process noise standard deviation longitudinal acceleration in m/s^2
    _std_yawdd = 1.0;  // Process noise standard deviation yaw acceleration in rad/s^2

    //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
    _std_laspx = 0.15;  // Laser measurement noise standard deviation position1 in m
    _std_laspy = 0.15;  // Laser measurement noise standard deviation position2 in m
    _std_radr = 0.3;  // Radar measurement noise standard deviation radius in m
    _std_radphi = 0.03;  // Radar measurement noise standard deviation angle in rad
    _std_radrd = 0.3;  // Radar measurement noise standard deviation radius change in m/s
    //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
    
    _is_initialized = false;
    _n_x = 5;
    _n_aug = 7;
    _lambda = 3 - _n_aug;
    
    _x = Eigen::VectorXd(_n_x);
    _P = Eigen::MatrixXd(_n_x, _n_x);
    
    _weights = Eigen::VectorXd(2 * _n_aug + 1);
    
    _x_aug = Eigen::VectorXd(_n_aug);
    for (unsigned int i=_n_x; i<_n_aug; ++i)
	_x_aug(i) = 0.;
    _P_aug = Eigen::MatrixXd(_n_aug, _n_aug);
    _P_aug.fill(0);
    _P_aug(5,5) = _std_a*_std_a;
    _P_aug(6,6) = _std_yawdd*_std_yawdd;
    
    _Xsig_pred = Eigen::MatrixXd(_n_x, 2 * _n_aug + 1);
    
    _nisData[0].limit = 5.99;
    _nisData[1].limit = 7.81;
}

UKF::~UKF()
{
}

void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{
    
    // IInitialization
    if (!_is_initialized)
    {
	this->Init(meas_package);
	_previous_timestamp = meas_package.timestamp_;
	_is_initialized = true;
	return;
    }
    
    // Prediction
    float dt = (meas_package.timestamp_ - _previous_timestamp) / 1000000.0;	//dt - expressed in seconds
    _previous_timestamp = meas_package.timestamp_;
    
    this->Prediction(dt);
    
   // Update
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
	if (!_use_radar)
	    return;
	
	// Radar updates
	this->UpdateRadar(meas_package);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
	if (!_use_laser)
	    return;

	// Laser updates
	this->UpdateLidar(meas_package);
    }
    
    this->PrintNIS();
    
    return;
}

Eigen::VectorXd UKF::x()
{
    return _x;
}

void UKF::Init(MeasurementPackage meas_package)
{
    // initial state vector
    if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
	_x(0) = meas_package.raw_measurements_(0);
	_x(1) = meas_package.raw_measurements_(1);
	_x(2) = 0;
	_x(3) = 0;
	_x(4) = 0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
	_x(0) = meas_package.raw_measurements_(0)*cos(meas_package.raw_measurements_(1));
	_x(1) = meas_package.raw_measurements_(0)*sin(meas_package.raw_measurements_(1));
	_x(2) = 0;
	_x(3) = 0;
	_x(4) = 0;
    }
    
    // initial covariance matrix
    if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
	_P << _std_laspx*_std_laspx, 0, 0, 0, 0,
	      0, _std_laspy*_std_laspy, 0, 0, 0,
	      0, 0, 10, 0, 0,
	      0, 0, 0, 10, 0,
	      0, 0, 0, 0, 10;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
	double a = _std_radr+_std_radphi;
	a = a*a;
	_P << a, 0, 0, 0, 0,
		    0, a, 0, 0, 0,
		    0, 0, 10, 0, 0,
		    0, 0, 0, 10, 0,
		    0, 0, 0, 0, 10;
    }
    
    
    //generate waights
    _weights.fill(1./(2*(_lambda+_n_aug)));
    _weights(0) *= 2*_lambda;

    _HLidar = Eigen::MatrixXd(2, _n_x);
    _HLidar << 1, 0, 0, 0, 0,
		       0, 1, 0, 0, 0;
    
    _RLidar = Eigen::MatrixXd(2, 2);
    _RLidar << _std_laspx*_std_laspx, 0,
			0, _std_laspy*_std_laspy;
    _ILidar = Eigen::MatrixXd::Identity(_n_x, _n_x);
    
    return;
}

void UKF::Prediction(double delta_t)
{
    //Generate augmented sigma points
    
    _x_aug.head(_n_x) = _x;
    _P_aug.topLeftCorner(_n_x, _n_x) = _P;
    
    Eigen::MatrixXd Xsig_aug = Eigen::MatrixXd(_n_aug, 2 * _n_aug + 1);
    Eigen::MatrixXd Psq = _P_aug.llt().matrixL();
    Xsig_aug.col(0) = _x_aug;
    for (unsigned int i=1; i<=_n_aug; ++i)
    {
	Eigen::VectorXd Psql = Psq.col(i-1)*sqrt(_lambda+_n_aug);
	Xsig_aug.col(i) = _x_aug + Psql;
	Xsig_aug.col(i+_n_aug) = _x_aug - Psql;
    }
        
    //Apply proccess function to sigma points
    for (unsigned int i=0; i<2*_n_aug+1; ++i)
    {
	for (unsigned int j=0; j<5; ++j)
	    _Xsig_pred(j,i) = Xsig_aug(j,i);
      
	double d0 = Xsig_aug(4,i)*delta_t;
	double d2 = Xsig_aug(3,i)+d0;
	double dt2 = delta_t*delta_t/2;
	
	if (fabs(Xsig_aug(4,i)) > 0.0001)
	{
	    double d1 = Xsig_aug(2,i)/Xsig_aug(4,i);
	    _Xsig_pred(0,i) += d1*(sin(d2)-sin(Xsig_aug(3,i)));
	    _Xsig_pred(1,i) += d1*(-cos(d2)+cos(Xsig_aug(3,i)));
	    _Xsig_pred(3,i) += d0;
	}
	else
	{
	    float d0 = Xsig_aug(2,i)*delta_t;
	    _Xsig_pred(0,i) += d0*cos(Xsig_aug(3,i));
	    _Xsig_pred(1,i) += d0*sin(Xsig_aug(3,i));
	}
	
	_Xsig_pred(0,i) += dt2*cos(Xsig_aug(3,i))*Xsig_aug(5,i);
	_Xsig_pred(1,i) += dt2*sin(Xsig_aug(3,i))*Xsig_aug(5,i);
	_Xsig_pred(2,i) += delta_t*Xsig_aug(5,i);
	_Xsig_pred(3,i) += dt2*Xsig_aug(6,i);
	_Xsig_pred(4,i) += delta_t*Xsig_aug(6,i);
    }
    
    // Recostruct mean and covariance matrix from predicted sigma points
  
    //predict state mean
    _x.fill(0.);
    for (unsigned int i=0; i<2 * _n_aug + 1; ++i)
      _x += _weights(i) * _Xsig_pred.col(i);
    
    //predict state covariance matrix
    _P.fill(0.);
    for (unsigned int i=0; i<2 * _n_aug + 1; ++i)
    {
	Eigen::VectorXd Xdiff = (_Xsig_pred.col(i)-_x);
	Tools::CutAngle(Xdiff(3));
// 	while (Xdiff(3)> M_PI) Xdiff(3)-=2.*M_PI;
// 	while (Xdiff(3)<-M_PI) Xdiff(3)+=2.*M_PI;
	
	_P += _weights(i)*Xdiff*Xdiff.transpose();
    }
        
    return;
}

void UKF::UpdateLidar(MeasurementPackage meas_package)
{
    Eigen::VectorXd y = meas_package.raw_measurements_ - _HLidar*_x;
    Eigen::MatrixXd S = _HLidar * _P * _HLidar.transpose() + _RLidar;
    Eigen::MatrixXd Sinv = S.inverse();
    Eigen::MatrixXd K = _P * _HLidar.transpose() * Sinv;
    
    //new estimate
    _x = _x + (K * y);
    _P = (_ILidar - K * _HLidar) * _P;
    
    double nis = y.transpose()*Sinv*y;
    this->CheckNIS(nis, MeasurementPackage::LASER);
    
    return;
}

void UKF::UpdateRadar(MeasurementPackage meas_package)
{
    //set measurement dimension, radar can measure r, phi, and r_dot
    int n_z = 3;
    
    Eigen::MatrixXd Zsig = Eigen::MatrixXd(n_z, 2 * _n_aug + 1);
    
    //mean predicted measurement
    Eigen::VectorXd z_pred = Eigen::VectorXd(n_z);
    
    //transform sigma points into measurement space
    for (unsigned int i=0; i<2*_n_aug+1; ++i)
    {
	float px = _Xsig_pred(0,i);
	float py = _Xsig_pred(1,i);
	float v = _Xsig_pred(2,i);
	float phi = _Xsig_pred(3,i);
	Zsig(0,i) = sqrt(px*px+py*py);
	Zsig(1,i) = atan2(py, px);
	Zsig(2,i) = v*(px*cos(phi) + py*sin(phi)) / Zsig(0,i);
    }

    //calculate mean predicted measurement
    z_pred.fill(0.);
    for (unsigned int i=0; i<2 * _n_aug + 1; ++i)
	z_pred += _weights(i)*Zsig.col(i);

    //calculate innovation covariance matrix S
    Eigen::MatrixXd S = Eigen::MatrixXd(n_z,n_z);
    S.fill(0.);
    for (unsigned int i=0; i<2 * _n_aug + 1; ++i)
    {
	Eigen::VectorXd Zdiff = (Zsig.col(i)-z_pred);
	Tools::CutAngle(Zdiff(1));
// 	while (Zdiff(1)> M_PI) Zdiff(1)-=2.*M_PI;
// 	while (Zdiff(1)<-M_PI) Zdiff(1)+=2.*M_PI;
	
	S += _weights(i)*Zdiff*Zdiff.transpose();
    }
    S(0,0) += _std_radr*_std_radr;
    S(1,1) += _std_radphi*_std_radphi;
    S(2,2) += _std_radrd*_std_radrd;
    
    Eigen::MatrixXd Sinv = S.inverse();
    
    //calculate cross correlation matrix
    Eigen::MatrixXd Tc = Eigen::MatrixXd(_n_x, n_z);
    Tc.fill(0.);
    for (unsigned int i=0; i<2 * _n_aug + 1; ++i)
    {
	Eigen::VectorXd Zdiff = (Zsig.col(i)-z_pred);
	Tools::CutAngle(Zdiff(1));
// 	while (Zdiff(1)> M_PI) Zdiff(1)-=2.*M_PI;
// 	while (Zdiff(1)<-M_PI) Zdiff(1)+=2.*M_PI;
	
	Eigen::VectorXd Xdiff = _Xsig_pred.col(i) - _x;
	Tools::CutAngle(Xdiff(3));
// 	while (Xdiff(3)> M_PI) Xdiff(3)-=2.*M_PI;
// 	while (Xdiff(3)<-M_PI) Xdiff(3)+=2.*M_PI;
	
	Tc += _weights(i)*Xdiff*Zdiff.transpose();
    }
    
    //calculate Kalman gain K;
    Eigen::MatrixXd K = Tc*Sinv;

    //update state mean and covariance matrix
    Eigen::VectorXd y = meas_package.raw_measurements_ - z_pred;
    _x += K*y;
    _P -= K*S*K.transpose();
    
    double nis = y.transpose()*Sinv*y;
    this->CheckNIS(nis, MeasurementPackage::RADAR);
    
    return;
}

void UKF::CheckNIS(double nis, MeasurementPackage::SensorType sensorType)
{
    unsigned int n;
    if (sensorType == MeasurementPackage::LASER)
	n = 0;
    else if (sensorType == MeasurementPackage::RADAR)
	n = 1;
    else
	return;
    
    if (nis > _nisData[n].limit)
	++_nisData[n].overLimit;
    ++_nisData[n].total;
    
    return;
}

void UKF::PrintNIS()
{
    std::cout << "NIS:";
    for (unsigned int i=0; i<2; ++i)
    {
	std::cout << " " << _nisData[i].overLimit/double(_nisData[i].total);
    }
    std::cout << " " << _nisData[0].total << std::endl;
}

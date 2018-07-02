
/** @file ukf.h
 *  @author Andrey N. Zabegaev <speench@gmail.com>
 */

#ifndef _ukf_h_
#define _ukf_h_

#include "Eigen/Dense"
#include "measurement_package.h"

class UKF
{
public:

    /** @brief Constructor
     */
    UKF();
    
    /** @brief Destructor
     */
    virtual ~UKF();
    
    /** @brief Process measurement
     *  @param[in] meas_package The latest measurement data of either radar or laser
     */
    void ProcessMeasurement(MeasurementPackage meas_package);
    
    /** @brief Get current state vector
     *  @return current state vector
     */
    Eigen::VectorXd x();
    
protected:
    
    struct nisData
    {
	double limit = 0;
	unsigned int overLimit = 0;
	unsigned int total = 0;
    };
    
    /** @brief Init filter
     *  @param[in] meas_package first measurement
     */
    void Init(MeasurementPackage meas_package);
    
    /** @brief Predicts sigma points, the state, and the state covariance matrix
     *  @param[in] delta_t Time between k and k+1 in s
     */
    void Prediction(double delta_t);
    
    /** @brief Updates the state and the state covariance matrix using a laser measurement
     *  @param[in] meas_package The measurement at k+1
     */
    void UpdateLidar(MeasurementPackage meas_package);
    
    /** @brief Updates the state and the state covariance matrix using a radar measurement
     *  @param[in] meas_package The measurement at k+1
     */
    void UpdateRadar(MeasurementPackage meas_package);
    
    /** @brief NIS calculation procedure
     *  @param[in] nis nis for last measurement
     *  @param[in] sensorType last measurement type
     */
    void CheckNIS(double nis, MeasurementPackage::SensorType sensorType);
    
    /** @brief NIS output
     */
    void PrintNIS();
    
    bool					_is_initialized;				//!< initially set to false, set to true in first call of ProcessMeasurement
    bool					_use_laser;					//!< if this is false, laser measurements will be ignored (except for init)
    bool					_use_radar;					//!< if this is false, radar measurements will be ignored (except for init)
    Eigen::VectorXd			_x;							//!< state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
    Eigen::VectorXd			_x_aug;						//!< augmented state vector
    Eigen::MatrixXd			_P;							//!< state covariance matrix
    Eigen::MatrixXd			_P_aug;						//!< augmented covariance matrix
    Eigen::MatrixXd			_Xsig_pred;					//!< predicted sigma points matrix
    Eigen::VectorXd			_weights;					//!< Weights of sigma points
    Eigen::MatrixXd			_HLidar;						//!< H matrix for lidar measurement
    Eigen::MatrixXd			_RLidar;						//!< H matrix for lidar measurement
    Eigen::MatrixXd			_ILidar;						//!< H matrix for lidar measurement
    long long				_time_us;					//!< time when the state is true, in us
    long long				_previous_timestamp;			//!< time of previouse update
    double					_std_a;						//!< Process noise standard deviation longitudinal acceleration in m/s^2
    double					_std_yawdd;					//!< Process noise standard deviation yaw acceleration in rad/s^2
    double					_std_laspx;					//!< Laser measurement noise standard deviation position1 in m
    double					_std_laspy;					//!< Laser measurement noise standard deviation position2 in m
    double					_std_radr;					//!< Radar measurement noise standard deviation radius in m
    double					_std_radphi;					//!< Radar measurement noise standard deviation angle in rad
    double					_std_radrd;					//!< Radar measurement noise standard deviation radius change in m/s
    int					_n_x;						//!< State dimension
    int					_n_aug;						//!< Augmented state dimension
    double					_lambda;						//!< Sigma point spreading parameter
    nisData				_nisData[2];					//!< data fro NIS calculation, for Laser and Radar measurements
};

#endif /* _ukf_h_ */

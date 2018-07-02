
/** @file tools.cpp
 *  @author Andrey N. Zabegaev <speench@gmail.com>
 */

#include "tools.h"

Tools::Tools()
{
}

Tools::~Tools()
{
}

Eigen::VectorXd Tools::CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, const std::vector<Eigen::VectorXd> &ground_truth)
{
    Eigen::VectorXd rmse(ground_truth[0].size());
    rmse.fill(0.);
    
    if ( (estimations.size() == 0) || (estimations.size() != ground_truth.size()) )
	return rmse;
    
    //accumulate squared residuals
    for(int i=0; i < estimations.size(); ++i)
    {
	Eigen::VectorXd d = estimations[i] - ground_truth[i];
	d = d.array()*d.array();
	rmse += d;
    }
    
    //calculate the mean
    rmse = rmse / estimations.size();
    
    //calculate the squared root
    rmse = rmse.array().sqrt();
    
    //return the result
    return rmse;
}

double & Tools::CutAngle(double &angle)
{
    if ((angle < -M_PI) || (angle > M_PI))
    {
	double a = angle + M_PI;
	int b = floor(a / (2*M_PI));
	double c = a - b*2*M_PI - M_PI;
	angle = c;
    }
    
    return angle;
}

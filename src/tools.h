
/** @file tools.h
 *  @author Andrey N. Zabegaev <speench@gmail.com>
 */

#ifndef _tools_h_
#define _tools_h_

#include <vector>
#include "Eigen/Dense"

class Tools
{
public:

    /** @brief Constructor
     */
    Tools();
    
    /** @brief Destructor
     */
    virtual ~Tools();

    /** @breif A helper method to calculate RMSE
     */
    Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, const std::vector<Eigen::VectorXd> &ground_truth);
    
    /** @brief cut value to -pi,pi range
     *  @param[in,out] angle value to cut
     *  @return cutted value
     */
    static double &CutAngle(double &angle);
};

#endif /* _tools_h_ */

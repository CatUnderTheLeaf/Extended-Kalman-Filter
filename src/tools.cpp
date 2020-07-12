#include "tools.h"
#include "math.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  if (estimations.size()==0 and estimations.size()!=ground_truth.size()) {
      std::cout << "Error";
  }
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  else {    
     VectorXd err(4);
     err << 0,0,0,0;

      for (int i=0; i < estimations.size(); ++i) {
        err = estimations[i] - ground_truth[i];
        err = err.array()*err.array();
        rmse = rmse + err;        
      }
        rmse = rmse / estimations.size();
        rmse = rmse.array().sqrt();
  }
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj(3,4);
  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);


  // check division by zero
  if (px==0 and py==0) {
      std::cout << "Error, division by zero!";
  }

  
  // compute the Jacobian matrix
  else {
      double sum_sqr = px*px + py*py;
      float min_vx = vx*py - vy*px;
      float min_vy = vy*px - vx*py;
      Hj << px/pow(sum_sqr, 0.5), py/pow(sum_sqr, 0.5), 0, 0,
        -py/sum_sqr, px/sum_sqr, 0, 0,
        py*min_vx/pow(sum_sqr, 1.5), px*(min_vy)/pow(sum_sqr, 1.5), px/pow(sum_sqr, 0.5), py/pow(sum_sqr, 0.5);
  }
  
  return Hj;
}
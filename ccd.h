// Software License Agreement (BSD License)
// 
//   Copyright (c) 2011, Shulei Zhu <schuleichu@gmail.com>
//   All rights reserved.
// 
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided that the following conditions
//   are met:
// 
//    * Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//    * Redistributions in binary form must reproduce the above
//      copyright notice, this list of conditions and the following
//      disclaimer in the documentation and/or other materials provided
//      with the distribution.
//    * Neither the name of Shulei Zhu nor the names of its
//      contributors may be used to endorse or promote products derived
//      from this software without specific prior written permission.
// 
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
//   FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
//   COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
//   BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
//   CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//   LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
//   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//   POSSIBILITY OF SUCH DAMAGE.
// 
// 
// ccd.h --- 
// File            : ccd.h
// Created: Sa Jun 18 14:06:36 2011 (+0200)
// Author: Shulei Zhu
// 
// Code:

#include <Eigen/Core>

/* 
 * #pragma warning (disable:981)        
 * #pragma warning (disable:383)
 * #pragma warning (disable:15)
 */
struct CCDParams
{
 CCDParams(): gamma_1(0.5), gamma_2(5), gamma_3(7), gamma_4(5),alpha(1.2), beta(0.06), kappa(0.5),c(0.25), h(20), delta_h(1),resolution(500), degree(4), phi_dim(3)
  {
  }
  CCDParams(double p1,
            double p2,
            double p3,
            double p4,
            double p5,
            double p6,
            double p7,
            double p8,
            int p9,
            int p10,
            int p11,
            int p12,
            int p13
            )
  {
    gamma_1 = p1;
    gamma_2 = p2;
    gamma_3 = p3;
    gamma_4 = p4;
    alpha = p5;
    beta = p6;
    kappa = p7;
    c = p8;
    h = p9;
    delta_h = p10;
    resolution = p11;
    degree = p12;
    phi_dim = p13;
  }

  ~CCDParams()
  {
  }
  double gamma_1;
  double gamma_2;
  double gamma_3;
  double gamma_4;
  double alpha;
  double beta;
  double kappa;
  double c;
  int h;
  int delta_h;
  int resolution;
  int degree;
  int phi_dim;
};

struct pointCCD
{
	        double X;
			double Y;
			double Z;
            double x;
            double y;
            double nx;
            double ny;
			int xu;
			int xv;
};

class CCD
{
public:
  //cv::Mat image, canvas, tpl;
  std::vector<cv::Point3d> pts;
  std::vector<cv::Point3d> displacements;
  std::vector< std::vector<double> > errors;
  
  std::vector<pointCCD> pointsccdmin;
  std::vector<pointCCD> pointsccd2;

  /* 
   * CCD()
   * {
   *   Phi = cv::Mat::zeros(params_.phi_dim,1, CV_64F);
   *   Sigma_Phi = cv::Mat::zeros(params_.phi_dim,params_.phi_dim, CV_64F);
   *   delta_Phi = cv::Mat::zeros(params_.phi_dim,1, CV_64F);
   * }
   */
  void read_params( const std::string& filename);
  void init_mat();
  void init(Eigen::Matrix3f &_rgbIntrinsicMatrix, std::vector<pointCCD> &pointsccd);
  void updatePoints(std::vector<pointCCD> &pointsccd);
  
  void local_statistics(std::vector<pointCCD> &pointsccd, cv::Mat &image);
  void local_statistics_all(std::vector<pointCCD> &pointsccd, cv::Mat &image);
  void refine_parameters(std::vector<pointCCD> &pointsccd, cv::Mat &image);
  void run_ccd();
  double resolution(){return params_.resolution;}
  double degree(){return params_.degree;}
  /* inline void init_pts(int init_method); */
  void init_cov(int degree);
  ~CCD(){clear();}
private:
  void clear();
  /* 
   * void contour_sift();
   * void contour_manually();
   */
  /* void on_mouse( int event, int x, int y, int flags, void* param ); */
  CCDParams params_;
  cv::Mat vic;
  cv::Mat mean_vic;
  cv::Mat cov_vic;
  cv::Mat nv;
  cv::Mat Phi;
  cv::Mat Sigma_Phi;
  cv::Mat delta_Phi;
  cv::Mat bs_old;
  cv::Mat nabla_E;
  cv::Mat hessian_E;
  Eigen::Matrix3f rgbIntrinsicMatrix;
};

/* 
 * inline void CCD::init_pts(int init_method)
 * {
 *   if(init_method == 1)
 *     contour_manually();
 *   else if(init_method == 2)
 *     contour_sift();
 *   if((int)pts.size() > params_.degree)
 *   {
 *     for (int i = 0; i < params_.degree; ++i)
 *       pts.push_back(pts[i]);
 *   }
 * }
 */

/* void on_mouse(int event, int x, int y, int flags, void* param ); */
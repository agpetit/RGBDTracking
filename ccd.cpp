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
// ccd.cpp --- 
// File            : ccd.cpp
// Created: Sa Jun 18 14:03:47 2011 (+0200)
// Author: Shulei Zhu

// Code:

#include <opencv/cv.h>
#include <opencv2/core.hpp>
#include <opencv2/features2d.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

#include <vector>
#include <iostream>
#include <algorithm>
#include "ccd.h"
cv::Mat canvas_tmp;

inline double logistic(double x)
{
  return 1.0/(1.0+exp(-x));
}

inline double probit(double x)
{
  return 0.5*(1+1/sqrt(2)*erf(x));
}

inline cv::Scalar random_color(CvRNG* rng)
{
  int color = cvRandInt(rng);
  return CV_RGB(color&255, (color>>8)&255, (color>>16)&255);
}

void CCD::read_params( const std::string& filename)
{
  cv::FileStorage fs(filename, cv::FileStorage::READ);
  params_.gamma_1 = double(fs["gamma_1"]);      
  params_.gamma_2 = double(fs["gamma_2"]);      
  params_.gamma_3 = double(fs["gamma_3"]);      
  params_.gamma_4 = double(fs["gamma_4"]);      
  params_.alpha   = double(fs["alpha"]);
  params_.beta   = double(fs["beta"]);
  params_.kappa   = double(fs["kappa"]);        
  params_.c       = double(fs["c"]);            
  params_.h       = 40;//int(fs["h"]);            
  params_.delta_h = int(fs["delta_h"]);      
  params_.resolution = int(fs["resolution"]);   
  params_.degree  = int(fs["degree"]);
  params_.phi_dim  = int(fs["phi_dim"]);
  // std::cerr<< params_.gamma_1<< " ";
  // std::cerr<< params_.gamma_2<< " ";
  // std::cerr<< params_.gamma_3<< " ";
  // std::cerr<< params_.gamma_4<< " ";
  // std::cerr<< params_.alpha<< " ";
  // std::cerr<< params_.kappa<< " ";
  // std::cerr<< params_.c<< " ";
  // std::cerr<< params_.delta_h<< " ";
  // std::cerr<< params_.resolution<< " ";
  // std::cerr<< params_.degree<< " ";
  // std::cerr<< params_.phi_dim<< std::endl;
}

void CCD::init_mat()
{
  Phi = cv::Mat::zeros(params_.phi_dim,1, CV_64F);
  Sigma_Phi = cv::Mat::zeros(params_.phi_dim,params_.phi_dim, CV_64F);
  delta_Phi = cv::Mat::zeros(params_.phi_dim,1, CV_64F);
}

/*void CCD::init_cov(int degree)
{
  int n_dim = (int)pts.size() - 3;
  cv::Mat W = cv::Mat::zeros(2*n_dim, params_.phi_dim, CV_64F);
  cv::Mat U = cv::Mat::zeros(2*n_dim, 2*n_dim, CV_64F);
  for (int i = 0; i < n_dim; ++i)
  {
    double *W_ptr = W.ptr<double>(i);
    W_ptr[0] = 1;
    W_ptr[1] = 0;
    W_ptr[2] = pts[i].x;
    W_ptr[3] = 0;
    W_ptr[4] = 0;
    W_ptr[5] = pts[i].y;
    if(params_.phi_dim == 8)
    {
      W_ptr[6] = pts[i].z;
      W_ptr[7] = 0;
    }
    W_ptr = W.ptr<double>(i+params_.resolution);
    W_ptr[0] = 0;
    W_ptr[1] = 1;
    W_ptr[2] = 0;
    W_ptr[3] = pts[i].y;
    W_ptr[4] = pts[i].x;
    W_ptr[5] = 0;
    if(params_.phi_dim == 8)
    {
      W_ptr[6] = 0;
      W_ptr[7] = pts[i].z;
    }
  }

  cv::Mat tmp_mat = cv::Mat::zeros(n_dim, n_dim, CV_64F);
  int interval = params_.resolution/(pts.size() - degree);

  for (int i = 0; i < n_dim; ++i)
  {
    double *basic_mat_ptr = bs.basic_mat_.ptr<double>(i*interval);
    for (int m = 0; m < n_dim; ++m)
    {
      double *tmp_mat_ptr = tmp_mat.ptr<double>(m);
      for (int n = 0; n < n_dim; ++n)
        tmp_mat_ptr[n] += basic_mat_ptr[m]*basic_mat_ptr[n];
    }
  }
  
  for (int i = 0; i < n_dim; ++i)
  {
    double *tmp_mat_ptr = tmp_mat.ptr<double>(i);
    double *U_ptr = U.ptr<double>(i);
    for (int j = 0; j < n_dim; ++j)
    {
      U_ptr[j] = tmp_mat_ptr[j]/n_dim;
      U_ptr = U.ptr<double>(i+n_dim);
      U_ptr[j+n_dim] = tmp_mat_ptr[j]/n_dim;
    }
  }

  cv::Mat tmp_cov;
  cv::gemm(W, U, params_.beta,  cv::Mat(), 0, tmp_cov, cv::GEMM_1_T);
  cv::gemm(tmp_cov, W, 1,  cv::Mat(), 0, Sigma_Phi, 0);
}*/


void CCD::clear()
{
  vic.release();
  mean_vic.release();
  cov_vic.release();
  nv.release();
  Phi.release();
  Sigma_Phi.release();
  delta_Phi.release();
  //bs_old.release();
  nabla_E.release();
  hessian_E.release();
  //image.release();
  //canvas.release();
  //if(!tpl.empty()) tpl.release();
}

void CCD::local_statistics(std::vector<pointCCD> &pointsccd, cv::Mat &image)
{
	
  cv::Mat_<cv::Vec3b>& img = (cv::Mat_<cv::Vec3b>&)image;
  
  //cv::imwrite("image.png", image);
  
  // std::cout << params_.gamma_1 << " " << params_.gamma_2 << " " << params_.gamma_3 << " " << params_.gamma_4 << std::endl;
  double sigma = params_.h/(params_.alpha*params_.gamma_3);
  // sigma_hat = gamma_3 * sigma
    
  //  double sigma_hat = max(h/sqrt(2*gamma_2), gamma_4);
  double sigma_hat = params_.gamma_3*sigma + params_.gamma_4;

  // to save the normalized parameters of vic[i,8]
  // dimension: resolution x 2
  // the first column save the normalized coefficient outside the curve
  // the second column store the one inside the curve
  cv::Mat normalized_param = cv::Mat::zeros(params_.resolution, 2, CV_64F);
  
  vic = cv::Mat::zeros(params_.resolution, 20*floor(params_.h/params_.delta_h), CV_64F);

  // temporary points used to store those points in the
  // normal direction as well as negative normal direction
  cv::Point3d tmp1, tmp2;

  // store the distance from a point in normal(negative norml) direction
  // to the point on the curve
  cv::Point3d tmp_dis1, tmp_dis2;

  CvRNG rng;
  cv::Scalar color = random_color(&rng);
  
  for(int i=0; i < params_.resolution;i++)
  {
    // cv::circle(canvas, cv::Point2d(bs[i].x, bs[i].y), 1,color, 1);
    double *nv_ptr = nv.ptr<double>(i);
    // normal vector (n_x, n_y)
    // tagent vector (nv.at<double>(i,1), -n_x)
	
    nv_ptr[0] = pointsccd[i].nx;//-bs.dt(i).y/cvSqrt(bs.dt(i).x*bs.dt(i).x + bs.dt(i).y*bs.dt(i).y);
    nv_ptr[1] = pointsccd[i].ny;//bs.dt(i).x/cvSqrt(bs.dt(i).x*bs.dt(i).x + bs.dt(i).y*bs.dt(i).y);
    // dimension = 4
    // 0 - x value
    // 1 - y value
    // 2 - normal vector x
    // 3 - normal vector y
    //bs_old = cv::Mat::zeros(params_.resolution, 4, CV_64F);
    //double *bs_old_ptr = bs_old.ptr<double>(i);

    // save the old value of bspline
    /*bs_old_ptr[0] = bs[i].x;
    bs_old_ptr[1] = bs[i].y;

    // save the old normal vector of bspline
    bs_old_ptr[2] = nv_ptr[0];
    bs_old_ptr[3] = nv_ptr[1];*/
    
     //std::cout << " normal " << i << " " << params_.h << " " << nv_ptr[1] << std::endl;
    int k = 0;
    double alpha = 0.5;
    double *vic_ptr = vic.ptr<double>(i);
    for (int j = params_.delta_h; j <= params_.h; j += params_.delta_h, k++)
    {
      ///////////////////////////////////////////////////////////////////////////////////////////
      // calculate in the direction +n: (n_x, n_y)
      /////////////////////////////////////////////////////////////////////////////////////////

      // x_{k,l}
      tmp1.x = round(pointsccd[i].xu + j*nv_ptr[0]);

      // y_{k,l}
      tmp1.y = round(pointsccd[i].xv + j*nv_ptr[1]);
	  
	  //std::cout << " tmp x " << tmp1.x << " " << tmp1.y << std::endl;

      // cv::circle(canvas_tmp, cv::Point2d(tmp1.x, tmp1.y), 1, CV_RGB(255,0,0), 1);

      // distance between x_{k,l} and x_{k,0} in the normal direction
      // appoximately it is l*h, l = {1,2,3,.....}
      tmp_dis1.x = (tmp1.x - pointsccd[i].xu)*nv_ptr[0] + (tmp1.y - pointsccd[i].xv)*nv_ptr[1];

      // distance between y_{k,l} and y_{k,0} along the curve
      // it approximates 0
      tmp_dis1.y = (tmp1.x - pointsccd[i].xu)*nv_ptr[1] - (tmp1.y - pointsccd[i].xv)*nv_ptr[0];
      
      vic_ptr[10*k + 0] = tmp1.y;
      vic_ptr[10*k + 1] = tmp1.x;
      vic_ptr[10*k + 2] = tmp_dis1.x;
      vic_ptr[10*k + 3] = tmp_dis1.y;

      // fuzzy assignment a(d_{k,l}) = 1/2*(erf(d_{kl})/\sqrt(2)*sigma) + 1/2
      // vic_ptr[10*k + 4] = 0.5*(erf((tmp_dis1.x)/(sqrt(2)*sigma)) + 1);
      vic_ptr[10*k + 4] = logistic(tmp_dis1.x/(sqrt(2)*sigma));

      // wp1 = (a_{d,l} - gamm_1) /(1-gamma_1)
      double wp1 = (vic_ptr[10*k + 4] - params_.gamma_1)/(1-params_.gamma_1);

      // wp1^4, why? if a_{d,l} \approx 0.5, do not count the point
      vic_ptr[10*k + 5] = wp1*wp1*wp1*wp1;

      // wp1 = (1-a_{d,l} - gamm_1) /(1-gamma_1)
      // double wp2 = (1-vic_ptr[10*k + 4] - gamma_1)/(1-gamma_1);
      double wp2 = (1-vic_ptr[10*k + 4] - 0.25);
      vic_ptr[10*k + 6] = -64*wp2*wp2*wp2*wp2 + 0.25;

      // W_p(d_p, simga_p) = c*max[0, exp(-d_p^2/2*sigma_p'^2) - exp(-gamma_2))]
      vic_ptr[10*k + 7] = std::max((exp(-0.5*tmp_dis1.x*tmp_dis1.x/(sigma_hat*sigma_hat)) - exp(-params_.gamma_2)), 0.0);

      // W' = 0.5*exp(-|d_v= - d_p=|/alpha)/alpha
      vic_ptr[ 10*k + 8] = 0.5*exp(-abs(tmp_dis1.y)/alpha)/alpha;
      
      // the derivative of col_5: 1/(sqrt(2*PI)*sigma)*exp{-d_{k,l}^2/(2*sigma*sigma)}
      vic_ptr[ 10*k + 9] = exp(-tmp_dis1.x*tmp_dis1.x/(2*sigma*sigma))/(sqrt(2*CV_PI)*sigma);

      // m1_debug[0] += img(tmp1.y, tmp1.x)[0];
      // m1_debug[1] += img(tmp1.y, tmp1.x)[1];
      // m1_debug[2] += img(tmp1.y, tmp1.x)[2];
      // calculate the normalization parameter c 
      normalized_param.at<double>(i, 0) += vic_ptr[ 10*k + 7];

        
#ifdef DEBUG
      if(i == 0)
        std::cout << "tmp1 " << tmp1.x  << " " << tmp1.y << std::endl;
#endif
      
      // cv::circle(img1, tmp1, 1, CV_RGB(0, 255, 255), 1, 8 , 0);

      // for (int m = 0; m < 3; ++m)
      // {
      //   vic_pixels.at<Vec3b>(i, k)[m] = img(tmp1.y, tmp1.x)[m];          
      // }

      ///////////////////////////////////////////////////////////////////////////////////////////
      // calculate in the direction -n: (-n_x, -n_y)
      /////////////////////////////////////////////////////////////////////////////////////////      
      tmp2.x = round(pointsccd[i].xu - j*nv_ptr[0]);//round(bs[i].x - j*nv_ptr[0]);
      tmp2.y = round(pointsccd[i].xv - j*nv_ptr[1]);//round(bs[i].y - j*nv_ptr[1]);
      // cv::circle(canvas_tmp, cv::Point2d(tmp2.x, tmp2.y), 1, CV_RGB(255,0,0), 1);
#ifdef DEBUG
      if(i == 0)
        std::cout << "tmp2 " << tmp2.x  << " " << tmp2.y << std::endl;
#endif

      // start compute the size in the direction of -(n_x, n_y)
      tmp_dis2.x = (tmp2.x - pointsccd[i].xu)*nv_ptr[0] + (tmp2.y - pointsccd[i].xv)*nv_ptr[1];
      tmp_dis2.y = (tmp2.x - pointsccd[i].xu)*nv_ptr[1] - (tmp2.y - pointsccd[i].xv)*nv_ptr[0];
      int negative_normal = k + (int)floor(params_.h/params_.delta_h);
      vic_ptr[10*negative_normal + 0] = tmp2.y;
      vic_ptr[10*negative_normal + 1] = tmp2.x;
      vic_ptr[10*negative_normal + 2] = tmp_dis2.x;
      vic_ptr[10*negative_normal + 3] = tmp_dis2.y;
      // vic_ptr[10*negative_normal + 4] = 0.5*(erf(tmp_dis2.x/(cvSqrt(2)*sigma)) + 1);
      vic_ptr[10*negative_normal + 4] = logistic(tmp_dis2.x/(sqrt(2)*sigma));
      // vic_ptr[10*negative_normal + 4] = 0.5;
      wp1 = (vic_ptr[10*negative_normal + 4] - 0.25);
      vic_ptr[10*negative_normal + 5] = -64*wp1*wp1*wp1*wp1 + 0.25;
      wp2 = (1 - vic_ptr[10*negative_normal + 4] - params_.gamma_1)/(1-params_.gamma_1);
      vic_ptr[10*negative_normal + 6] = wp2*wp2*wp2*wp2;
      vic_ptr[10*negative_normal + 7] = std::max((exp(-0.5*tmp_dis2.x*tmp_dis2.x/(sigma_hat*sigma_hat)) - exp(-params_.gamma_2)), 0.0);
      vic_ptr[ 10*negative_normal + 8] = 0.5*exp(-abs(tmp_dis2.x)/alpha)/alpha;
      vic_ptr[ 10*negative_normal + 9] = exp(-tmp_dis2.x*tmp_dis2.x/(2*sigma*sigma))/(sqrt(2*CV_PI)*sigma);
      normalized_param.at<double>(i, 1) += vic_ptr[ 10*negative_normal + 7];
    }
  }
  

#ifdef DEBUG
  printf("%-5s  %-5s  %-5s  %-5s  %-5s  %-5s  %-5s  %-5s  %-5s  %-5s\n",
         "x", "y", "dist_x", "dist_y", "a", "w1^4", "w2^4", "prox", "edf", "erf'"
         );
  for (int  i = 0; i < 20*floor(params_.h/params_.delta_h); ++i)
  {
    // std::cout << vic.at<double>(0,i) << "    ";
    printf("%-5f   ", vic.at<double>(0,i));
    if((i+1)%10 == 0)
      std::cout << std::endl;
  }
#endif

  for (int i = 0; i < params_.resolution; ++i)
  {
    int k = 0;
    // w1 = \sum wp_1, w2 = \sum wp_2
    double w1 =0.0 , w2 = 0.0;

    // store mean value near the curve
    std::vector<double> m1(3,0.0), m2(3,0.0);
    
    // store the second mean value near the curve
    std::vector<double> m1_o2(9,0.0), m2_o2(9,0.0);

    ////////////////////////////////////////////////////////////////////////
    // compute local statistics
    // ////////////////////////////////////////////////////////////////////
    // start search the points in the +n direction as well as -n direction
    double wp1 = 0.0, wp2 = 0.0;

    double *vic_ptr = vic.ptr<double>(i);
    double *mean_vic_ptr = mean_vic.ptr<double>(i);
    double *cov_vic_ptr = cov_vic.ptr<double>(i);
    for (int j = params_.delta_h; j <= params_.h; j += params_.delta_h, k++)
    {
      wp1 = 0.0, wp2 = 0.0;
      int negative_normal = k + (int)floor(params_.h/params_.delta_h);
      
      // wp1 = w(a_{k,l})*w(d_{k,l})*w(d)
      wp1 = vic_ptr[ 10*k+ 5]*vic_ptr[ 10*k+ 7]/normalized_param.at<double>(i,0);

      // wp2 = w(a_{k,l})*w(d_{k,l})*w(d)
      wp2 = vic_ptr[ 10*k+ 6]*vic_ptr[ 10*k+ 7]/normalized_param.at<double>(i,1);;

      //w1 = \sum{wp1}
      w1 += wp1;

      //w2 = \sum{wp2}
      w2 += wp2;

      // compute the mean value in the vicinity of a point
      // m_{ks} = I{k}^{s} = \sum_{l} w_{kls}{I_{kl}} : s = 1 or 2

      m1[0] += wp1*img(vic_ptr[ 10*k + 0 ], vic_ptr[ 10*k + 1 ])[0];
      m1[1] += wp1*img(vic_ptr[ 10*k + 0 ], vic_ptr[ 10*k + 1 ])[1];
      m1[2] += wp1*img(vic_ptr[ 10*k + 0 ], vic_ptr[ 10*k + 1 ])[2];
      m2[0] += wp2*img(vic_ptr[ 10*k + 0 ], vic_ptr[ 10*k + 1 ])[0];
      m2[1] += wp2*img(vic_ptr[ 10*k + 0 ], vic_ptr[ 10*k + 1 ])[1];
      m2[2] += wp2*img(vic_ptr[ 10*k + 0 ], vic_ptr[ 10*k + 1 ])[2];

      // compute second order local statistics
      // m_{k,s} = \sum_{l} w_{kls} I_{kl}*I_{kl}^T

      for (int m = 0; m < 3; ++m)
      {
        for (int n =0; n < 3; ++n)
        {
          m1_o2[m*3+n] += wp1*img(vic_ptr[ 10*k + 0 ], vic_ptr[ 10*k + 1 ])[m]
                          *img(vic_ptr[ 10*k + 0 ], vic_ptr[ 10*k + 1 ])[n];
          m2_o2[m*3+n] += wp2*img(vic_ptr[ 10*k + 0 ], vic_ptr[ 10*k + 1 ])[m]
                          *img(vic_ptr[ 10*k + 0 ], vic_ptr[ 10*k + 1 ])[n];
        }
      }
        
      wp1 = vic_ptr[ 10*negative_normal+ 5]*vic_ptr[ 10*negative_normal+ 7]/normalized_param.at<double>(i,0);
      wp2 = vic_ptr[ 10*negative_normal+ 6]*vic_ptr[ 10*negative_normal+ 7]/normalized_param.at<double>(i,1);
      
      w1 += wp1;
      w2 += wp2;
      
      m1[0] += wp1*img(vic_ptr[10*negative_normal+0], vic_ptr[10*negative_normal+1])[0];
      m1[1] += wp1*img(vic_ptr[10*negative_normal+0], vic_ptr[10*negative_normal+1])[1];
      m1[2] += wp1*img(vic_ptr[10*negative_normal+0], vic_ptr[10*negative_normal+1])[2];
      m2[0] += wp2*img(vic_ptr[10*negative_normal+0], vic_ptr[10*negative_normal+1])[0];
      m2[1] += wp2*img(vic_ptr[10*negative_normal+0], vic_ptr[10*negative_normal+1])[1];
      m2[2] += wp2*img(vic_ptr[10*negative_normal+0], vic_ptr[10*negative_normal+1])[2];
      
      for (int m = 0; m < 3; ++m)
      {
        for (int n =0; n < 3; ++n)
        {
          m1_o2[m*3+n] += wp1*img(vic_ptr[ 10*negative_normal + 0 ], vic_ptr[ 10*negative_normal + 1 ])[m]
                          *img(vic_ptr[ 10*negative_normal + 0 ], vic_ptr[ 10*negative_normal + 1 ])[n];
          m2_o2[m*3+n] += wp2*img(vic_ptr[ 10*negative_normal + 0 ], vic_ptr[ 10*negative_normal + 1 ])[m]
                          *img(vic_ptr[ 10*negative_normal + 0 ], vic_ptr[ 10*negative_normal + 1 ])[n];
        }
      }
    }
    // std::cout << "w1: " << "              w2:" << std::endl;
    // std::cout << "w1 == " << w1 << "  w2== " << w2 << std::endl;
    mean_vic_ptr[0] = m1[0]/w1;
    mean_vic_ptr[1] = m1[1]/w1;
    mean_vic_ptr[2] = m1[2]/w1;
    mean_vic_ptr[3] = m2[0]/w2;
    mean_vic_ptr[4] = m2[1]/w2;
    mean_vic_ptr[5] = m2[2]/w2;
    
    for (int m = 0; m < 3; ++m)
    {
      for (int n = 0 ; n < 3; ++n)
      {
        cov_vic_ptr[ m*3+n] = m1_o2[m*3+n]/w1 -m1[m]*m1[n]/(w1*w1);
        cov_vic_ptr[ 9+m*3+n] = m2_o2[m*3+n]/w2 -m2[m]*m2[n]/(w2*w2);
        if(m == n)
        {
          cov_vic_ptr[ m*3+n] += params_.kappa;
          cov_vic_ptr[ 9+m*3+n] += params_.kappa;
        }
      }
    }
  }
  normalized_param.release();
}


void CCD::local_statistics_all(std::vector<pointCCD> &pointsccd, cv::Mat &image)
{
    //int k = 0;
    //double alpha = 0.5;
	//cv::imwrite("imageccd.png", image);
	
  cv::Mat_<cv::Vec3b>& img = (cv::Mat_<cv::Vec3b>&)image;
  cv::Mat tmp_cov = cv::Mat::zeros(3,3,CV_64F);
  cv::Mat tmp_cov_inv = cv::Mat::zeros(3,3,CV_64F);
  cv::Mat tmp_jacobian = cv::Mat::zeros(params_.phi_dim,3,CV_64F);
  cv::Mat tmp_pixel_diff = cv::Mat::zeros(3, 1, CV_64F);
	
	errors.resize(params_.h);
	pointsccdmin.resize(params_.resolution);
	pointsccd2.resize(params_.resolution);
	std::vector<pointCCD> pointsccd1;
	pointsccd1.resize(params_.resolution);

	//std::cout << " " << pointsccd2.size() << std::endl;
	
for (int j = 0; j < params_.h/2; j++/*= params_.delta_h*/)
    {
		
		//std::cout << " ok params " << " " << (int)2*floor(params_.h/params_.delta_h) << std::endl;
	for(int i0=0; i0 < params_.resolution;i0++)
	{    		
		
		//std::cout << " ok res 1 " << i0 << " " << params_.resolution << std::endl;

		    // cv::circle(canvas, cv::Point2d(bs[i].x, bs[i].y), 1,color, 1);
    double *nv_ptr = nv.ptr<double>(i0);
    // normal vector (n_x, n_y)
    // tagent vector (nv.at<double>(i,1), -n_x)
	
    nv_ptr[0] = pointsccd[i0].nx;//-bs.dt(i).y/cvSqrt(bs.dt(i).x*bs.dt(i).x + bs.dt(i).y*bs.dt(i).y);
    nv_ptr[1] = pointsccd[i0].ny;//bs.dt(i).x/cvSqrt(bs.dt(i).x*bs.dt(i).x + bs.dt(i).y*bs.dt(i).y);
	
	pointsccdmin[i0].nx = nv_ptr[0] ;
	pointsccdmin[i0].ny = nv_ptr[1] ;
	
	pointsccd2[i0].nx = nv_ptr[0] ;
	pointsccd2[i0].ny = nv_ptr[1] ;
	
	pointsccd1[i0].nx = nv_ptr[0] ;
	pointsccd1[i0].ny = nv_ptr[1] ;
		
    // dimension = 4
    // 0 - x value
    // 1 - y value
    // 2 - normal vector x
    // 3 - normal vector y
    //bs_old = cv::Mat::zeros(params_.resolution, 4, CV_64F);
    //double *bs_old_ptr = bs_old.ptr<double>(i);

    // save the old value of bspline
    
     //std::cout << " normal " << i << " " << params_.h << " " << nv_ptr[1] << std::endl;
   // int k = 0;
    //double alpha = 0.5;
    double *vic_ptr = vic.ptr<double>(i0);
      ///////////////////////////////////////////////////////////////////////////////////////////
      // calculate in the direction +n: (n_x, n_y)
      /////////////////////////////////////////////////////////////////////////////////////////
	 // j= 0;
	  	  // x_{k,l}
      pointsccd2[i0].xu = round(pointsccd[i0].xu - j*nv_ptr[0]);
      // y_{k,l}
      pointsccd2[i0].xv = round(pointsccd[i0].xv - j*nv_ptr[1]);
      // x_{k,l}
      pointsccd1[i0].xu = round(pointsccd[i0].xu + j*nv_ptr[0]);
      // y_{k,l}
      pointsccd1[i0].xv = round(pointsccd[i0].xv + j*nv_ptr[1]);
	  
	  //std::cout << " ok res 1 " << j << " " << params_.resolution << std::endl;
	}
	
	//std::cout << " ok local stat " << std::endl;
	local_statistics(pointsccd1, image);
	
	errors[j].resize(params_.resolution);
	
  for (int i1 = 0; i1 < params_.resolution; i1++)
  {
	  	//std::cout << " ok sort " << i1 << std::endl;

    double *vic_ptr = vic.ptr<double>(i1);
    double *nv_ptr = nv.ptr<double>(i1);
    double *mean_vic_ptr = mean_vic.ptr<double>(i1);
    double *cov_vic_ptr = cov_vic.ptr<double>(i1);
    double normal_points_number = floor(params_.h/params_.delta_h);
	
		  	//std::cout << " ok sort 1 " << std::endl;
	
	nabla_E = cv::Mat::zeros(params_.phi_dim,1, CV_64F);
	hessian_E = cv::Mat::zeros(params_.phi_dim,params_.phi_dim, CV_64F);
	
		  	//std::cout << " force 20 " << std::endl;
	
    for (int kk = 0; kk < 2*normal_points_number; ++kk)
    {
      tmp_cov = cv::Mat::zeros(3,3,CV_64F);
      tmp_cov_inv = cv::Mat::zeros(3,3,CV_64F);
	        
      for (int m = 0; m < 3; ++m)
      {
        double *tmp_cov_ptr = tmp_cov.ptr<double>(m);
        for (int n = 0; n < 3; ++n)
        {
          tmp_cov_ptr[n] = vic_ptr[10*kk+4] * cov_vic_ptr[m*3+n] +(1-vic_ptr[10*kk+4])* cov_vic_ptr[m*3+n+9];
        }
      }
      
      tmp_cov_inv = tmp_cov.inv(cv::DECOMP_SVD);
 
      tmp_pixel_diff = cv::Mat::zeros(3, 1, CV_64F);

      //compute the difference between I_{kl} and \hat{I_{kl}}
      for (int m = 0; m < 3; ++m)
      {
	//if (j==0&&m==1) std::cout << " mean " <<(int)img(vic_ptr[10*k+0], vic_ptr[10*k+1])[m] << " " << vic_ptr[10*k+4] << " " <<  mean_vic_ptr[m] << " " << mean_vic_ptr[m+3] << std::endl;
        tmp_pixel_diff.at<double>(m,0) = img(vic_ptr[10*kk+0], vic_ptr[10*kk+1])[m]- vic_ptr[10*kk+4] * mean_vic_ptr[m]- (1-vic_ptr[10*kk+4])* mean_vic_ptr[m+3];
	  }

      //\nabla{E_2} = \sum J * \Sigma_{kl}^{-1} * (I_{kl} - \hat{I_{kl}})
      nabla_E += /*tmp_cov_inv**/tmp_pixel_diff;
    }
	//std::cout << " yy " << (int)2*floor(params_.h/params_.delta_h)+1 << " " << j<< std::endl;
	//std::cout << " ok sort 2 " << nabla_E.at<double>(0,0) << " " << nabla_E.at<double>(1,0) << " " << nabla_E.at<double>(2,0) << std::endl;
  
	errors[j][i1] = (nabla_E.at<double>(0,0)*nabla_E.at<double>(0,0) + nabla_E.at<double>(1,0)*nabla_E.at<double>(1,0) + nabla_E.at<double>(2,0)*nabla_E.at<double>(2,0));	
    //if (i == 30) std::cout << " xx " << nabla_E << std::endl;

  }
 // std::cout << " local stat00 " << std::endl;

  	errors[(params_.h/(2*params_.delta_h)) + j].resize(params_.resolution);
	
		//double timeccd0 = (double)getTickCount();
	
  	local_statistics(pointsccd2, image);
	
		//double timeccd1 = ((double)getTickCount() - timeccd0)/getTickFrequency();
		
//std::cout << " local stat " << std::endl;
  
  for (int i1 = 0; i1 < params_.resolution; ++i1)
  {
	  	//std::cout << " ok sort " << std::endl;

    double *vic_ptr = vic.ptr<double>(i1);
    double *nv_ptr = nv.ptr<double>(i1);
    double *mean_vic_ptr = mean_vic.ptr<double>(i1);
    double *cov_vic_ptr = cov_vic.ptr<double>(i1);
    double normal_points_number = floor(params_.h/params_.delta_h);
	
		  	//std::cout << " ok sort 1 " << std::endl;
	
	nabla_E = cv::Mat::zeros(params_.phi_dim,1, CV_64F);
	hessian_E = cv::Mat::zeros(params_.phi_dim,params_.phi_dim, CV_64F);
	
		  	//std::cout << " force 20 " << std::endl;
	
    for (int k0 = 0; k0 < 2*normal_points_number; k0++)
    {
      tmp_cov = cv::Mat::zeros(3,3,CV_64F);
      tmp_cov_inv = cv::Mat::zeros(3,3,CV_64F);
	        
      for (int m = 0; m < 3; ++m)
      {
        double *tmp_cov_ptr = tmp_cov.ptr<double>(m);
        for (int n = 0; n < 3; n++)
        {
          tmp_cov_ptr[n] = vic_ptr[10*k0+4] * cov_vic_ptr[m*3+n] +(1-vic_ptr[10*k0+4])* cov_vic_ptr[m*3+n+9];
        }
      }
      
      tmp_cov_inv = tmp_cov.inv(cv::DECOMP_SVD);
 
      tmp_pixel_diff = cv::Mat::zeros(3, 1, CV_64F);

      //compute the difference between I_{kl} and \hat{I_{kl}}
      for (int m = 0; m < 3; ++m)
      {
	//std::cout << " yy " << (int)img(vic_ptr[10*k+0], vic_ptr[10*k+1])[m] << " " << vic_ptr[10*k+4] << " " << mean_vic_ptr[m] << std::endl;
        tmp_pixel_diff.at<double>(m,0) = img(vic_ptr[10*k0+0], vic_ptr[10*k0+1])[m]- vic_ptr[10*k0+4] * mean_vic_ptr[m]- (1-vic_ptr[10*k0+4])* mean_vic_ptr[m+3];
	  }

      //\nabla{E_2} = \sum J * \Sigma_{kl}^{-1} * (I_{kl} - \hat{I_{kl}})
      nabla_E += /*tmp_cov_inv**/tmp_pixel_diff;
    }
	//std::cout << " yy " << (int)2*floor(params_.h/params_.delta_h)+1 << " " << j<< std::endl;
	//std::cout << " ok sort 2 " << nabla_E.at<double>(0,0) << " " << nabla_E.at<double>(1,0) << " " << nabla_E.at<double>(2,0) << std::endl;
  
	errors[(params_.h/(2*params_.delta_h))+j][i1] = (nabla_E.at<double>(0,0)*nabla_E.at<double>(0,0) + nabla_E.at<double>(1,0)*nabla_E.at<double>(1,0) + nabla_E.at<double>(2,0)*nabla_E.at<double>(2,0));	
    //if (i1 == 30) std::cout << " yy " << errors[(params_.h/params_.delta_h)+j][i1] << std::endl;
  }
  }
  
  std::vector<int> indices_min_error;
  indices_min_error.resize(params_.resolution);
  
  double error_sum= 0;
  
	for(int ii=0; ii < params_.resolution;ii++)
    {
		
		double min = 1000000;
		int ind_min = 0;
		
		error_sum+=errors[0][ii];
		for (int j = 0; j < params_.h/2; j++/*= params_.delta_h*/)
		{
		if (ii == 30) std::cout << " xx1 " << errors[j][ii] << std::endl;
          if (errors[j][ii] < min)
		  {
			  min = errors[j][ii];
			  ind_min = j;
		  }
		  
		  //std::cout << " errors 1 " << errors[j][i] << std::endl;

		}
		
		for (int j = 0; j < params_.h/2; j ++)
		{
    if (ii == 30) std::cout << " yy1 " << errors[(params_.h/(2*params_.delta_h))+j][ii] << std::endl;

		  if (errors[(params_.h/(2*params_.delta_h)) + j][ii] < min)
		  {
			  min = errors[(params_.h/(2*params_.delta_h)) + j][ii];
			  ind_min = -j;
		  }
		  
		  		 // std::cout << " errors 2 " << errors[(params_.h/params_.delta_h) + j][i] << " min " << min << " ind_min " << ind_min << std::endl;
		  
		  //std::cout << " errors " << errors[(params_.h/params_.delta_h) +j][i] << std::endl;
		}
		//ind_min = 10;
	  pointsccdmin[ii].xu = round(pointsccd[ii].xu + ind_min*pointsccd[ii].nx);

      // y_{k,l}
      pointsccdmin[ii].xv = round(pointsccd[ii].xv + ind_min*pointsccd[ii].ny);
	  
	  
	  
	if (ii == 30) std::cout << " errors 2 " << ind_min << " min " << pointsccdmin[ii].xu << " ind_min " << pointsccdmin[ii].xv<< std::endl;
	  indices_min_error[ii] = ind_min;
	  
	  //std::cout << " ind min " << indices_min_error[i] << std::endl; 
	}
	
		  std::cout << " erorr sum " << error_sum << std::endl; 
	

}

void CCD::refine_parameters(std::vector<pointCCD> &pointsccd, cv::Mat &image)
{
  cv::Mat_<cv::Vec3b>& img = (cv::Mat_<cv::Vec3b>&)image;
  cv::Mat tmp_cov = cv::Mat::zeros(3,3,CV_64F);
  cv::Mat tmp_cov_inv = cv::Mat::zeros(3,3,CV_64F);
  cv::Mat tmp_jacobian = cv::Mat::zeros(params_.phi_dim,3,CV_64F);
  cv::Mat tmp_pixel_diff = cv::Mat::zeros(3, 1, CV_64F);
    
  displacements.resize(params_.resolution);
  
  // std::cout << "dimension: " << Sigma_Phi.cols << " " << Sigma_Phi.rows << std::endl;
  
  std::cout << " params_.resolution " << params_.resolution << std::endl;

  for (int i = 0; i < params_.resolution; ++i)
  {
    double *vic_ptr = vic.ptr<double>(i);
    double *nv_ptr = nv.ptr<double>(i);
    double *mean_vic_ptr = mean_vic.ptr<double>(i);
    double *cov_vic_ptr = cov_vic.ptr<double>(i);
    double normal_points_number = floor(params_.h/params_.delta_h);
	
	nabla_E = cv::Mat::zeros(params_.phi_dim,1, CV_64F);
	hessian_E = cv::Mat::zeros(params_.phi_dim,params_.phi_dim, CV_64F);
	
		  	//std::cout << " force 20 " << std::endl;
	
    for (int j = 0; j < 2*normal_points_number; ++j)
    {
      tmp_cov = cv::Mat::zeros(3,3,CV_64F);
      tmp_cov_inv = cv::Mat::zeros(3,3,CV_64F);
      
      for (int m = 0; m < 3; ++m)
      {
        double *tmp_cov_ptr = tmp_cov.ptr<double>(m);
        for (int n = 0; n < 3; ++n)
        {
          tmp_cov_ptr[n] = vic_ptr[10*j+4] * cov_vic_ptr[m*3+n] +(1-vic_ptr[10*j+4])* cov_vic_ptr[m*3+n+9];
        }
      }
      
      tmp_cov_inv = tmp_cov.inv(cv::DECOMP_SVD);
 
      tmp_pixel_diff = cv::Mat::zeros(3, 1, CV_64F);

      //compute the difference between I_{kl} and \hat{I_{kl}}
      for (int m = 0; m < 3; ++m)
      {
        tmp_pixel_diff.at<double>(m,0) = img(vic_ptr[10*j+0], vic_ptr[10*j+1])[m]- vic_ptr[10*j+4] * mean_vic_ptr[m]- (1-vic_ptr[10*j+4])* mean_vic_ptr[m+3];
      }
      
      //compute jacobian matrix        
      tmp_jacobian = cv::Mat::zeros(params_.phi_dim,3,CV_64F);
 
      for (int n = 0; n < 3; ++n)
      {
		tmp_jacobian.at<double>(0,n) = rgbIntrinsicMatrix(0,0)*vic_ptr[10*j + 9]*(mean_vic_ptr[n] - mean_vic_ptr[n+3])*(nv_ptr[0]/pointsccd[i].Z);
		tmp_jacobian.at<double>(1,n) = rgbIntrinsicMatrix(0,0)*vic_ptr[10*j + 9]*(mean_vic_ptr[n] - mean_vic_ptr[n+3])*(nv_ptr[1]/pointsccd[i].Z);
		tmp_jacobian.at<double>(2,n) = rgbIntrinsicMatrix(0,0)*vic_ptr[10*j + 9]*(mean_vic_ptr[n] - mean_vic_ptr[n+3])*(-(nv_ptr[0]*pointsccd[i].x+nv_ptr[1]*pointsccd[i].y)/pointsccd[i].Z);
		//std::cout << " jaco " <<pointsccd[i].x << " " << pointsccd[i].y << " " << pointsccd[i].Z  << std::endl;
         //std::cout <<" mean vic " << mean_vic_ptr[n] << " " << mean_vic_ptr[n+3] << std::endl;
      }
      
      //\nabla{E_2} = \sum J * \Sigma_{kl}^{-1} * (I_{kl} - \hat{I_{kl}})
      nabla_E += tmp_jacobian*tmp_cov_inv*tmp_pixel_diff;
      //Hessian{E_2} = \sum J * \Sigma_{kl}^{-1} * J
      hessian_E += tmp_jacobian*tmp_cov_inv*tmp_jacobian.t();
    }

    cv::Mat hinvnab = hessian_E.inv(cv::DECOMP_SVD)*nabla_E;
	displacements[i].x = 0.001*nabla_E.at<double>(0,0);
	displacements[i].y = 0.001*nabla_E.at<double>(1,0);
	displacements[i].z = 0.001*nabla_E.at<double>(2,0);
		
			//std::cout << " force 0 " << displacements[i].x << " " << displacements[i].y << " " << displacements[i].z << std::endl;
			
			//std::cout << " force 0 " << displacements[i].x << " "<< hinvnab.at<double>(0,0) << std::endl;

  }
  
  //cv::Mat Sigma_Phi_inv = Sigma_Phi.inv(cv::DECOMP_SVD);
  //hessian_E += Sigma_Phi_inv;
  //nabla_E += 2*Sigma_Phi_inv*Phi;
  
    	  	  	//std::cout << " ok stats 2 " << std::endl;


  cv::Mat hessian_E_inv;//= hessian_E.inv(cv::DECOMP_SVD);
  //delta_Phi = hessian_E_inv*nabla_E;
  // #ifdef DEBUG
  // std::cout << delta_Phi.at<double>(0,0) << " "
  //           << delta_Phi.at<double>(1,0) << " " 
  //           << delta_Phi.at<double>(2,0) << " " 
  //           << delta_Phi.at<double>(3,0) << " " 
  //           << delta_Phi.at<double>(4,0) << " " 
  //           << delta_Phi.at<double>(5,0) << " ";
  // if(params_.phi_dim == 8)
  // {
  //   std::cout<< delta_Phi.at<double>(6,0) << " " 
  //            << delta_Phi.at<double>(7,0) << " ";
  // }
  // std::cout<< std::endl;
  // #endif
  // cv::norm(delta_Phi);
  
  //std::cout << " ok " << displacements[0].x << std::endl;
  
  /*Phi -= delta_Phi;
  Sigma_Phi = params_.c*Sigma_Phi + 2*(1-params_.c)*hessian_E_inv;

  Sigma_Phi_inv.release();
  hessian_E_inv.release();*/
  tmp_cov.release();
  tmp_cov_inv.release();
  tmp_jacobian.release();
  tmp_pixel_diff.release();
}

void CCD::init(Eigen::Matrix3f &_rgbIntrinsicMatrix, std::vector<pointCCD> &pointsccd)
{
  // store the control points trasformed in the shape space
  int iter = 0;
  double tol = 0.0;
  double tol_old = 0.0;
  bool convergence = false;
  double norm = 0.0;
  
  	rgbIntrinsicMatrix = _rgbIntrinsicMatrix;
	/*Phi = cv::Mat::zeros(ccdParameters.phi_dim,1, CV_64F);
	Sigma_Phi = cv::Mat::zeros(ccdParameters.phi_dim,ccdParameters.phi_dim, CV_64F);
	sigmaF.resize(6,6);
	sigmaP.resize(6,6);
	delta_Phi = cv::Mat::zeros(ccdParameters.phi_dim,1, CV_64F);
    tol = 0.0;
    tol_old = 0.0;
    convergence = false;
    norm = 0.0;
	resolution = pointsCCD[scaleLevel].size();
    nv = cv::Mat::zeros(resolution, 2, CV_64F);
    mean_vic = cv::Mat::zeros(resolution, 6, CV_64F);
    cov_vic = cv::Mat::zeros(resolution, 18, CV_64F);

    nv_prev = cv::Mat::zeros(resolution, 2, CV_64F);
    mean_vic_prev = cv::Mat::zeros(resolution, 6, CV_64F);
    cov_vic_prev = cv::Mat::zeros(resolution, 18, CV_64F);
    vic_prev = cv::Mat::zeros(resolution, 20*floor(ccdParameters.h/ccdParameters.delta_h), CV_64F);

    nabla_E = cv::Mat::zeros(ccdParameters.phi_dim,1, CV_64F);
    hessian_E = cv::Mat::zeros(ccdParameters.phi_dim,ccdParameters.phi_dim, CV_64F);
    hessian_E_partial = cv::Mat::zeros(ccdParameters.phi_dim,ccdParameters.phi_dim, CV_64F);
    vic = cv::Mat::zeros(resolution, 20*floor(ccdParameters.h/ccdParameters.delta_h), CV_64F);
    p_old = cv::Mat::zeros(resolution, 4, CV_64F);
	//int normal_points_number = floor(ccdParameters.h/ccdParameters.delta_h);
	nerror_ccd = 2*normal_points_number*3*resolution;*/
	
	params_.resolution = pointsccd.size();
	
	nv = cv::Mat::zeros(params_.resolution, 2, CV_64F);
    mean_vic = cv::Mat::zeros(params_.resolution, 6, CV_64F);
    cov_vic = cv::Mat::zeros(params_.resolution, 18, CV_64F);
    nabla_E = cv::Mat::zeros(params_.phi_dim,1, CV_64F);
    hessian_E = cv::Mat::zeros(params_.phi_dim,params_.phi_dim, CV_64F);

  //do{
    
	// update model parameters
    // for (int i = 0; i < 6; ++i)
    //   Phi.at<double>(i,0) = Phi.at<double>(i,0) - delta_Phi.at<double>(i,0);

    // update the control points in terms of the change of
    // model parameters    
    /*for (size_t i = 0; i < pts.size(); ++i)
    {
      // C = W*\Phi + C_0
      //           1  0  x_0  0  0  y_0
      //     C = [                       ][\phi_0 \phi_1 \phi_2 \phi_3 \phi_4 \phi_5 ]^T + C_0
      //           0  1  0   y_0 x_0  0
      //
      pts[i].x = Phi.at<double>(0,0) + (1+Phi.at<double>(2,0))*pts[i].x + Phi.at<double>(5,0)*pts[i].y;
      if(params_.phi_dim == 8)
        pts[i].x +=  Phi.at<double>(6,0)*pts[i].z;
      pts[i].y = Phi.at<double>(1,0) + (1+Phi.at<double>(3,0))*pts[i].y + Phi.at<double>(4,0)*pts[i].x;
      if(params_.phi_dim == 8)
        pts[i].y += Phi.at<double>(7,0)*pts[i].z;
      pts[i].z = pts[i].z;
    }*/

    // create a new B-spline curve
    //BSpline bs(params_.degree , params_.resolution, pts);
    
	//image.copyTo(canvas_tmp);

    /*for (int i = 0; i < params_.resolution; ++i)
    {
        int j = (i+1)%params_.resolution;
        cv::line(canvas_tmp, cv::Point2d(bs[i].x, bs[i].y),cv::Point2d(bs[j].x, bs[j].y),CV_RGB(255, 0, 0 ),2,8,0);
        // cv::line(canvas, cv::Point2d(bs[i].x, bs[i].y),cv::Point2d(bs[j].x, bs[j].y),CV_RGB( 0, 0, 255 ),2,8,0);
    }*/

    // converge condition
    // tol = \int (r - r_f)*n
    tol = 0.0;
    /*if(iter > 0)
    {
      for (int k = 0; k < params_.resolution; ++k)
      {
        tol += pow((bs[k].x - bs_old.at<double>(k, 0))*bs_old.at<double>(k, 2) +
                   (bs[k].y - bs_old.at<double>(k, 1))*bs_old.at<double>(k, 3), 2);
      }
      // if(iter > 1)
      // params_.h = std::max(params_.h/log(tol_old/tol),10.0);
      tol_old = tol;
    }*/
    /*local_statistics();
    
    refine_parameters();*/

    /*norm =  0.0;
    for (int i = 0; i < params_.phi_dim; ++i){
      if(i == 0 || i == 1)
        norm += delta_Phi.at<double>(i, 0)*delta_Phi.at<double>(i, 0)/10000;
      else
        norm += delta_Phi.at<double>(i, 0)*delta_Phi.at<double>(i, 0);
    }
    norm = cv::sqrt(norm);*/
    //std::cerr << "iter: " << iter << "   tol: " << tol  << " norm: " << cv::norm(delta_Phi, cv::NORM_L2)  << " norm_tmp:" << norm << " params.h: " << params_.h << std::endl;
    // if(iter == 19)
    //   for (int i = 0; i < params_.resolution; ++i){
    //     int j = (i+1)%params_.resolution;
    //     cv::line(canvas, cv::Point2d(bs[i].x, bs[i].y),cv::Point2d(bs[j].x, bs[j].y),CV_RGB( 0, 0, 255 ),2,8,0);        
    //   }

    /*std::stringstream name;
    name << iter;
    cv::imwrite(name.str() + ".png", canvas_tmp);
    canvas_tmp.release();*/
    // cv::imwrite(name.str() + ".png", canvas);

    // cv::imshow("CCD", canvas);    
    // cv::waitKey(200);

    /*if(iter >=1 && tol < 1 && params_.h > 10)
    {
      params_.h = params_.h/sqrt(2);
    }*/
    
    //if(iter >= 1 && tol < 0.001)
    {
      /*for (int i = 0; i < params_.resolution; ++i)
      {
        // std::cout << bs[i].x << " " << bs[i].y << " " <<bs[i].z << std::endl;
        int j = (i+1)%params_.resolution;
        // cv::line(canvas_tmp, cv::Point2d(bs[i].x, bs[i].y),cv::Point2d(bs[j].x, bs[j].y),CV_RGB( 255, 0, 0 ),2,8,0);
        cv::line(canvas, cv::Point2d(bs[i].x, bs[i].y),cv::Point2d(bs[j].x, bs[j].y),CV_RGB( 255, 0, 0 ),2,8,0);
        // cv::circle(canvas_tmp, cv::Point2d(bs[i].x, bs[i].y), 2 ,CV_RGB(255,0,0), 2); 
        // cv::circle(canvas, cv::Point2d(bs[i].x, bs[i].y), 1 ,CV_RGB(0,255,0), 1); 
      }*/
      
	  //convergence = true;
      //init_cov(params_.degree);
    }
    //iter += 1;
    // bs.release();
  //}while(!convergence);
}
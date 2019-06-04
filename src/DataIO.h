/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 MGH, INRIA, USTL, UJF, CNRS                    *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this library; if not, write to the Free Software Foundation,     *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
*******************************************************************************
*                               SOFA :: Modules                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/

#ifndef SOFA_RGBDTRACKING_DATAIO_H
#define SOFA_RGBDTRACKING_DATAIO_H

#include <pcl/common/common_headers.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/search/impl/search.hpp>



#include <RGBDTracking/config.h>
#include <sofa/core/core.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/objectmodel/Event.h>
#include <sofa/simulation/AnimateBeginEvent.h>
#include <sofa/simulation/AnimateEndEvent.h>

#include <sofa/core/behavior/BaseForceField.h>
#include <sofa/core/behavior/ForceField.h>
#include <sofa/defaulttype/VecTypes.h>

#include <opencv/cv.h>
#include <opencv2/core.hpp>
#include <opencv2/features2d.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

#include <set>


#ifdef WIN32
    #include <process.h>
#else
    #include <pthread.h>
#endif

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <string>


#include <sys/times.h>

#include <boost/thread.hpp>

//using namespace std;
//using namespace cv;

namespace sofa {

namespace rgbdtracking {

//using helper::vector;
//using namespace sofa::defaulttype;

template<class DataTypes>
class DataIOInternalData {
public:
};

/*!
 * \brief The RGBDFileSystemIO class
 * util class used by data io
 * to store rgbd images in files
 */
class RGBDFileSystemIO {
public :
    RGBDFileSystemIO (std::string outpath, int iter_p)
        : depth_path     (outpath + "/depth1%06d.png")
        , depth_path2    (outpath + "/depth2%06d.png")
        , depthfile_path (outpath + "/depthfile%06d.txt")

        , image_path  (outpath + "/img1%06d.png")
        , segimg_path (outpath + "/imgseg1%06d.png")

        , rtt_path            (outpath + "/rtt%06d.png")
        , rttstress_path      (outpath + "/rttstress%06d.png") //deprecated
        , rttstressplast_path (outpath + "/rttstressplast%06d.png") //deprected

        , klt_path (outpath +  + "/imgklt%06d.png")
        , iter (iter_p)
    {}

    /// reads rgb image
    cv::Mat read_image () {
        return cv::imread(process_filename(image_path, iter));
    }
    cv::Mat read_segimage () {
        return cv::imread(process_filename(segimg_path, iter));
    }
    cv::Mat read_depths () {
        return cv::imread(process_filename(segimg_path, iter));
    }
    /// write rgbd
    void write_image (const cv::Mat &I) {
        cv::imwrite(process_filename(image_path, iter), I) ;
    }
    void write_segimage (const cv::Mat &I) {
        cv::imwrite(process_filename(segimg_path, iter), I) ;
    }
    void write_depth (const cv::Mat &I) {
        cv::imwrite(process_filename(depth_path, iter), I) ;
    }
    void write_depth_bis (const cv::Mat &I) {
        cv::imwrite(process_filename(depth_path2, iter), I) ;
    }
    void write_klt (const cv::Mat & I) {
        cv::imwrite(process_filename(klt_path, iter), I) ;
    }
    void write_rtt (const cv::Mat & I) {
        cv::imwrite(process_filename(rtt_path, iter), I) ;
    }
    void write_rtt_stress (const cv::Mat & I) {
        cv::imwrite(process_filename(rttstress_path, iter), I) ;
    }
    void write_rtt_stress_plast (const cv::Mat & I) {
        cv::imwrite(process_filename(rttstressplast_path, iter), I) ;
    }


    /// depth_files IO
    cv::Mat read_depth_file() {
        //create the file stream
        std::string filename(process_filename(depthfile_path, iter));
        std::ifstream file(filename.c_str(), std::ios::in | std::ios::binary );
        if (!file) {
            std::cerr
                << "(RGBDFileSystem) Error: "
                << "can't open the file stream to read depth data"
            << std::endl ;
            return cv::Mat::zeros (1,1, CV_32F) ;
        }

        //declare image parameters
        int matWidth = 0,
            matHeight = 0,
            type = 0 ;
        //read type and size of the matrix first
        file.read((char*) &type, sizeof(type));
        file.read((char*) &matWidth, sizeof(matWidth));
        file.read((char*) &matHeight, sizeof(matHeight));

        //change Mat type
        cv::Mat I = cv::Mat::zeros(matHeight, matWidth, type);

        //std::cout << " width " << matWidth << " " << matHeight <<  std::endl;

        //write data depending on the image's type
        switch (type) {
            // FLOAT ONE CHANNEL
            case CV_32F:
                std::cout << "Reading CV_32F image" << std::endl;
                read_image_from_stream<float>(I, matHeight, matWidth, file);
                std::cout << "Read CV_32F image" << std::endl;
            break;
            // DOUBLE ONE CHANNEL
            case CV_64F:
                std::cout << "Reading CV_64F image" << std::endl;
                read_image_from_stream<double>(I, matHeight, matWidth, file);
            break;

            // FLOAT THREE CHANNELS
            case CV_32FC3:
                std::cout << "Reading CV_32FC3 image" << std::endl;
                read_image_from_stream<defaulttype::Vec3f>(I, matHeight, matWidth, file);
            break;

            // DOUBLE THREE CHANNELS
            case CV_64FC3:
                std::cout << "Reading CV_64FC3 image" << std::endl;
                read_image_from_stream<defaulttype::Vec3d>(I, matHeight, matWidth, file);
            break;
            default:
                std::cerr
                    << "(RGBDFileSystem) Error: wrong Mat type" << std::endl
                    << "must be CV_32F, CV_64F, CV_32FC3 or CV_64FC3"
                << std::endl;
            break;
        }
        //close file
        file.close();

        return I;
    }
    bool write_depth_file (const cv::Mat &I) {
        std::string filename(process_filename(depthfile_path, iter));
        //create the file stream
        std::ofstream file(filename.c_str(), std::ios::out | std::ios::binary );
        if (!file) {
            std::cerr
                << "(RGBDFileSystem) Error: "
                << "can't open the file stream to read depth data"
            << std::endl ;
            return false ;
        }

        //load the matrix param
        int
            matWidth = I.size().width,
            matHeight = I.size().height,
            type = I.type();

        //write type and size of the matrix first
        file.write((const char*) &type, sizeof(type));
        file.write((const char*) &matWidth, sizeof(matWidth));
        file.write((const char*) &matHeight, sizeof(matHeight));

        //write data depending on the image's type
        switch (type) {
            // FLOAT ONE CHANNEL
            case CV_32F:
                std::cout << "Writing CV_32F image" << std::endl;
                write_image_from_stream<float>(I, matHeight, matWidth, file);
            break;
            // DOUBLE ONE CHANNEL
            case CV_64F:
                std::cout << "Writing CV_64F image" << std::endl;
                write_image_from_stream<double>(I, matHeight, matWidth, file);
            break;

            // FLOAT THREE CHANNELS
            case CV_32FC3:
                std::cout << "Writing CV_32FC3 image" << std::endl;
                write_image_from_stream<defaulttype::Vec3f>(I, matHeight, matWidth, file);
            break;
            // DOUBLE THREE CHANNELS
            case CV_64FC3:
                std::cout << "Writing CV_64FC3 image" << std::endl;
                write_image_from_stream<defaulttype::Vec3d>(I, matHeight, matWidth, file);
            break;
            default:
                std::cerr
                    << "(RGBDFileSystem) Error: wrong Mat type" << std::endl
                    << "must be CV_32F, CV_64F, CV_32FC3 or CV_64FC3"
                << std::endl;
            break;
        }

        //close file
        file.close();

        return 0;
    }


private :
    template <typename Type>
    void write_image_from_stream(const cv::Mat &I, int matHeight, int matWidth, std::ofstream & file) {
        for (int i=0; i < matWidth*matHeight; ++i) {
            Type value = I.at<Type>(i);
            file.write((const char*) &value, sizeof(value));
        }
    }

    template<typename Type>
    void read_image_from_stream(cv::Mat & I, int matHeight, int matWidth, std::ifstream & file) {
        for (int i=0; i < matWidth*matHeight; ++i) {
            Type value ;
            file.read((char*) &value, sizeof(value));
            I.at<Type>(i) = value;
        }
    }

    inline std::string process_filename (std::string path, int index) {
        char buf[FILENAME_MAX];
        sprintf(buf, path.c_str(), index);
        return std::string (buf) ;
    }

    std::string depth_path ;
    std::string depthfile_path ;
    std::string depth_path2 ;

    std::string image_path ;
    std::string segimg_path ;

    std::string rtt_path ;
    std::string rttstress_path ;
    std::string rttstressplast_path ;

    std::string klt_path ;

    int iter ;
} ;


template<class DataTypes>
class DataIO : public virtual sofa::core::objectmodel::BaseObject {
public:
    SOFA_CLASS(SOFA_TEMPLATE(DataIO,DataTypes), sofa::core::objectmodel::BaseObject);

    // Basic typedefs
    typedef sofa::core::objectmodel::BaseObject Inherit;
    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecReal VecReal;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef Data<typename DataTypes::VecCoord> DataVecCoord;
    typedef Data<typename DataTypes::VecDeriv> DataVecDeriv;

    // Paths
    Data<std::string> inputPath;
    Data<std::string> outputPath;

    // Number of iterations
    Data<int> nimages;
    Data<int> startimage;
    Data<int> niterations;

    int npasses;

    bool pcl;
    bool disp;

//    Data<bool> useGroundTruth;
//    Data<bool> useKLTPoints;
    Data<bool> newImages; // is output

    int ntargetcontours;
	
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr target;
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr targetP;
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr targetGt;

    int ind;
    cv::Mat depth;
    cv::Mat color ;
    cv::Mat color_1 ;

    int hght,wdth;
	
//    Data< VecCoord > targetPositions;
//    Data< VecCoord > targetGtPositions;

    std::vector<cv::Mat*> listimg;
    std::vector<cv::Mat*> listimgklt;
    std::vector<cv::Mat*> listimgseg;
    std::vector<cv::Mat*> listdepth;
    std::vector<cv::Mat*> listrtt;
    std::vector<cv::Mat*> listrttstress;
    std::vector<cv::Mat*> listrttstressplast;

//    std::vector<std::vector<defaulttype::Vec3d>*> listpcd;
//    std::vector<std::vector<bool>*> listvisible;

//    std::vector<std::vector<double>*> listvm;
//    std::vector<std::vector<double>*> listps;
//    std::vector<std::vector<double>*> listes;
//    std::vector<std::vector<double>*> listts;
//    std::vector<std::vector<double>*> listesnode;
//    std::vector<std::vector<double>*> listpsnode;
//    std::vector<std::vector<double>*> listtsnode;
	
//    std::vector<double> vonmisesstressGt;
//    std::vector<double> elasticstrainsGt;
//    std::vector<double> plasticstrainsGt;
//    std::vector<double> totalstrainsGt;
//    std::vector<double> elasticstrainsnodeGt;
//    std::vector<double> plasticstrainsnodeGt;
//    std::vector<double> totalstrainsnodeGt;

    cv::Mat* imgl;
//    cv::Mat* imgklt;
    cv::Mat* imglsg;
    cv::Mat* depthl;
    cv::Mat* rtt;
	
    int iter_im;
    cv::Mat rtd;
	
    DataIO();
    virtual ~DataIO();

    void init();
    void handleEvent(sofa::core::objectmodel::Event *event);
	
    void readImages();
	
    void writeImages();
		
};


#if defined(SOFA_EXTERN_TEMPLATE) && !defined(DataIO_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_RGBDTRACKING_API DataIO<defaulttype::Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_RGBDTRACKING_API DataIO<defaulttype::Vec3fTypes>;
#endif
#endif


} // rgbdtracking

} // namespace sofa

#endif

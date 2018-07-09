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

#define SOFA_RGBDTRACKING_CLOSESTPOINT_INL

#include <limits>
#include <iterator>
#include <sofa/helper/gl/Color.h>
#include <sofa/simulation/Simulation.h>
#include <iostream>
#include <map>

#ifdef USING_OMP_PRAGMAS
    #include <omp.h>
#endif

#include <algorithm> 
#include "ClosestPoint.h"

using std::cerr;
using std::endl;

namespace sofa
{

namespace core
{

namespace objectmodel
{

using namespace sofa::defaulttype;
using namespace helper;

template <class DataTypes>
ClosestPoint<DataTypes>::ClosestPoint()
    : Inherit()
    , blendingFactor(initData(&blendingFactor,(Real)1,"blendingFactor","blending between projection (=0) and attraction (=1) forces."))
    , outlierThreshold(initData(&outlierThreshold,(Real)7,"outlierThreshold","suppress outliers when distance > (meandistance + threshold*stddev)."))
    , rejectBorders(initData(&rejectBorders,false,"rejectBorders","ignore border vertices."))
    , useContour(initData(&useContour,false,"useContour","Emphasize forces close to the target contours"))
    , useVisible(initData(&useVisible,true,"useVisible","Use the vertices of the viisible surface of the source mesh"))
    , useDistContourNormal(initData(&useDistContourNormal,false,"useVisible","Use the vertices of the visible surface of the source mesh"))
{
    iter_im = 0;
    Vector4 camParam = cameraIntrinsicParameters.getValue();
	
    rgbIntrinsicMatrix(0,0) = camParam[0];
    rgbIntrinsicMatrix(1,1) = camParam[1];
    rgbIntrinsicMatrix(0,2) = camParam[2];
    rgbIntrinsicMatrix(1,2) = camParam[3];
	
}

template <class DataTypes>
ClosestPoint<DataTypes>::~ClosestPoint()
{
}

template <class DataTypes>
void ClosestPoint<DataTypes>::init()
{

    this->Inherit::init();

    // Get source normals
    //if(!sourceNormals.getValue().size()) serr<<"normals of the source model not found"<<sendl;		
}

template<class DataTypes>
void ClosestPoint<DataTypes>::detectBorder(vector<bool> &border,const helper::vector< tri > &triangles)
{
    unsigned int nbp=border.size();
    unsigned int nbt=triangles.size();
    for(unsigned int i=0;i<nbp;i++) border[i]=false;

    if(!nbt) return;
    vector<vector< unsigned int> > ngbTriangles((int)nbp);
    for(unsigned int i=0;i<nbt;i++) for(unsigned int j=0;j<3;j++)	ngbTriangles[triangles[i][j]].push_back(i);
    for(unsigned int i=0;i<nbp;i++) if(ngbTriangles[i].size()==0) border[i]=true;
    for(unsigned int i=0;i<nbt;i++)
        for(unsigned int j=0;j<3;j++)
        {
            unsigned int id1=triangles[i][j],id2=triangles[i][(j==2)?0:j+1];
            if(!border[id1] || !border[id2]) {
                bool bd=true;
                for(unsigned int i1=0;i1<ngbTriangles[id1].size() && bd;i1++)
                    for(unsigned int i2=0;i2<ngbTriangles[id2].size() && bd;i2++)
                        if(ngbTriangles[id1][i1]!=i)
                            if(ngbTriangles[id1][i1]==ngbTriangles[id2][i2])
                                bd=false;
                if(bd) border[id1]=border[id2]=true;
            }
        }
}

template<class DataTypes>
void ClosestPoint<DataTypes>::initSource()
{
    // build k-d tree
    const VecCoord&  p = sourcePositions.getValue();
    sourceKdTree.build(p);
	
    // detect border
   /* if(sourceBorder.size()!=p.size())
	{ sourceBorder.resize(p.size()); 
	//detectBorder(sourceBorder,sourceTriangles.getValue()); 
    }*/
}

template<class DataTypes>
void ClosestPoint<DataTypes>::initSourceSurface()
{
    // build k-d tree
    const VecCoord&  p = sourcePositions.getValue();
	
	ReadAccessor< Data< VecCoord > > pSurface(sourceSurfacePositions);
	ReadAccessor< Data< VecCoord > > nSurface(sourceSurfaceNormals);
	VecCoord nSurfaceM;
		
	sourceSurface.resize(p.size());
	Vector3 pos;
	Vector3 col;
	nSurfaceM.resize(p.size());

	for (unsigned int i=0; i<p.size(); i++)
	{
			sourceSurface[i] = false;
			for (unsigned int j=0; j<pSurface.size(); j++)
			{
				Real dist=(p[i]-pSurface[j]).norm();
				if ((float)dist < 0.001)
				{
					sourceSurface[i] = true;
					nSurfaceM[i] = nSurface[j];
				}
			}
	}
			
	sourceSurfaceNormalsM.setValue(nSurfaceM);
	
}

template<class DataTypes>
void ClosestPoint<DataTypes>::updateSourceSurface()
{
    // build k-d tree
    const VecCoord&  p = sourcePositions.getValue();
    const VecCoord&  pSurface = sourceSurfacePositions.getValue();
		
    sourceSurface.resize(p.size());
    Vector3 pos;
    Vector3 col;

        for (unsigned int i=0; i<p.size(); i++)
	{
            sourceSurface[i] = false;
                for (unsigned int j=0; j<pSurface.size(); j++)
                {
                    Real dist=(p[i]-pSurface[j]).norm();
                        if (dist < 10e-3f)
                        {
                            sourceSurface[i] = true;
                        }
                }
	}
	
}

template<class DataTypes>
void ClosestPoint<DataTypes>::initSourceVisible()
{
    // build k-d tree
	
    const VecCoord&  p = sourceVisiblePositions.getValue();
	
    sourceKdTree.build(p);

    // detect border
    /*if(sourceBorder.size()!=p.size())
    {
        sourceBorder.resize(p.size());

	//detectBorder(sourceBorder,sourceTriangles.getValue()); 
    }*/
}

template<class DataTypes>
void ClosestPoint<DataTypes>::initTarget()
{
    const VecCoord&  p = targetPositions.getValue();
	
    targetKdTree.build(p);

    // updatebbox
    for(unsigned int i=0;i<p.size();++i)    targetBbox.include(p[i]);

    // detect border
    //if(targetBorder.size()!=p.size()) { targetBorder.resize(p.size()); detectBorder(targetBorder,targetTriangles.getValue()); }

}

template<class DataTypes>
void ClosestPoint<DataTypes>::initTargetContour()
{

    const VecCoord&  p = targetContourPositions.getValue();
    targetContourKdTree.build(p);

    // updatebbox
    //for(unsigned int i=0;i<p.size();++i)    targetBbox.include(p[i]);

    // detect border
    //if(targetBorder.size()!=p.size()) { targetBorder.resize(p.size()); detectBorder(targetBorder,targetTriangles.getValue()); }

}

template<class DataTypes>
void ClosestPoint<DataTypes>::updateClosestPoints()
{
    VecCoord x;

    std::cout<<" source size " <<std::endl;

	if (!useVisible.getValue() || timer <= 2)
            x = sourcePositions.getValue();
        else
            x = sourceVisiblePositions.getValue();
	
    const VecCoord&  tp = targetPositions.getValue();
    unsigned int nbs=x.size(), nbt=tp.size();

    distanceSet emptyset;

    std::cout<<" source size "<< nbs <<std::endl;
	
        if(nbs!=closestSource.size()) {if (!useVisible.getValue() || timer <= 2) initSource(); else initSourceVisible();  closestSource.resize(nbs);	closestSource.fill(emptyset); cacheDist.resize(nbs); cacheDist.fill((Real)0.); cacheDist2.resize(nbs); cacheDist2.fill((Real)0.); previousX.assign(x.begin(),x.end());}
	
	/*if(nbtc!=closestSourceContour.size()) {initSource();  closestSourceContour.resize(nbtc);	
	closestSourceContour.fill(emptyset); 
	cacheDist.resize(nbtc); 
	cacheDist.fill((Real)0.); 
	cacheDist2.resize(nbtc); 
	cacheDist2.fill((Real)0.); 
        previousX.assign(x.begin(),x.end());}*/

        if(nbt!=closestTarget.size()) {initTarget();  closestTarget.resize(nbt);	closestTarget.fill(emptyset);}
	
        if(blendingFactor.getValue()<1)
        {

        //unsigned int count=0;
#ifdef USING_OMP_PRAGMAS
#pragma omp parallel for
#endif
            for(int i=0;i<(int)nbs;i++)
            {
            /*Real dx=(previousX[i]-x[i]).norm();
            //  closest point caching [cf. Simon96 thesis]
            if(dx>=cacheDist[i] || closestSource[i].size()==0)
            {
                targetKdTree.getNClosest(closestSource[i],x[i],this->cacheSize.getValue() );
                typename distanceSet::iterator it0=closestSource[i].begin(), it1=it0; it1++;
                typename distanceSet::reverse_iterator itn=closestSource[i].rbegin();
                cacheDist[i] =((itn->first)-(it0->first))*(Real)0.5;
                cacheDist2[i]=((it1->first)-(it0->first))*(Real)0.5;
                previousX[i]=x[i];
            }
            else if(dx>=cacheDist2[i]) // in the cache -> update N-1 distances 
            {
                targetKdTree.updateCachedDistances(closestSource[i],x[i]);
                //count++;
            }*/
            targetKdTree.getNClosest(closestSource[i],x[i],targetPositions.getValue(),1);
            //std::cout << " ok kdtree " << x[i][0] << " " << x[i][1] << std::endl;
            }

        //std::cout<<(Real)count*(Real)100./(Real)nbs<<" % cached"<<std::endl;
        }
		
        // closest source points from target points
        if(blendingFactor.getValue()>0)
        {
            if (!useVisible.getValue() || timer <= 2) initSource();
            else initSourceVisible();

        /*#ifdef USING_OMP_PRAGMAS
            #pragma omp parallel for
        #endif*/
            for(int i=0;i<(int)nbt;i++)
            {
                //if(!targetBackground[i])
                {
                    sourceKdTree.getNClosest(closestTarget[i],tp[i],sourcePositions.getValue(),1);
                }
            }
        }
    this->sourceIgnored.resize(nbs); sourceIgnored.fill(false);
    this->targetIgnored.resize(nbt); targetIgnored.fill(false);

    // prune outliers
        if(outlierThreshold.getValue()!=0 && timer > 5)
        {
        Real mean=0,stdev=0,count=0;
            for(unsigned int i=0;i<nbs;i++) if(closestSource[i].size() )
            {
                count++; stdev+=(closestSource[i].begin()->first)*(closestSource[i].begin()->first);
		mean+=(Real)(closestSource[i].begin()->first);
		//std::cout << " distances source " << nbs << " " << (double)(closestSource[i].begin()->first) << std::endl;
            }
            for(unsigned int i=0;i<nbt;i++) if(closestTarget[i].size())
            {
                count++;
                stdev+=(closestTarget[i].begin()->first)*(closestTarget[i].begin()->first);
                mean+=(Real)(closestTarget[i].begin()->first);
            }
		
        mean=mean/count; 
        stdev=(Real)sqrt(stdev/count-mean*mean);
        mean+=stdev*outlierThreshold.getValue();
        //std::cout << " distances " << count << " " << stdev << " " << mean << " " << outlierThreshold.getValue() << std::endl;
        //mean*=mean;
        for(unsigned int i=0;i<nbs;i++)
        {
            if(closestSource[i].size()) if(closestSource[i].begin()->first>mean )
            sourceIgnored[i]=true;
        }
        for(unsigned int i=0;i<nbt;i++) if(closestTarget[i].size()) if(closestTarget[i].begin()->first>mean )
        {
            targetIgnored[i]=true;
        }
		
        for(unsigned int i=0;i<nbs;i++)
            if(closestSource[i].size())
                for(unsigned int j=0;j<nbt;j++)
                    if(j == closestSource[i].begin()->second && i == closestTarget[j].begin()->second)
                    {
					//sourceIgnored[i]=true;
					//targetIgnored[j] = true;
                    }
			
        for(unsigned int i=0;i<nbs;i++)
            if(closestSource[i].size())// && sourceBorder[i])
                for(unsigned int j=0;j<nbt;j++)
                {
                    if(i == closestTarget[j].begin()->second)
                    {
                        if(closestSource[i].begin()->first < closestTarget[j].begin()->first)//&& i == closestTarget[j].begin()->second)
                        {
                            //sourceIgnored[i]=true;
                            //targetIgnored[j] = true;
                        }
                            //else targetIgnored[j] = true;
					
                    }
                }

        }

        if(rejectBorders.getValue())
        {
            for(unsigned int i=0;i<nbs;i++) if(closestSource[i].size()) if(targetBorder[closestSource[i].begin()->second]) sourceIgnored[i]=true;
            for(unsigned int i=0;i<nbt;i++) if(closestTarget[i].size()) if(sourceBorder[closestTarget[i].begin()->second]) targetIgnored[i]=true;
        }
    /*if(normalThreshold.getValue()>(Real)-1. && sourceNormals.getValue().size()!=0 && targetNormals.getValue().size()!=0) {
        ReadAccessor< Data< VecCoord > > sn(sourceNormals);
        ReadAccessor< Data< VecCoord > > tn(targetNormals);
        for(unsigned int i=0;i<nbs;i++) if(closestSource[i].size()) if(dot(sn[i],tn[closestSource[i].begin()->second])<normalThreshold.getValue()) sourceIgnored[i]=true;
        for(unsigned int i=0;i<nbt;i++) if(closestTarget[i].size()) if(dot(tn[i],sn[closestTarget[i].begin()->second])<normalThreshold.getValue()) targetIgnored[i]=true;
    }*/
}

template<class DataTypes>
void ClosestPoint<DataTypes>::updateClosestPointsContours()
{
	
    VecCoord x,x0;
    std::cout << " source size 2 " <<targetPositions.getValue().size() << " " << sourceVisiblePositions.getValue().size() << std::endl;
    if (!useVisible.getValue() || timer <= 2)
    x = sourcePositions.getValue();
    else  x = sourceVisiblePositions.getValue();

    x0 = sourcePositions.getValue();
	
    const VecCoord& tp = targetPositions.getValue();
    const VecCoord& xcp = sourceContourPositions.getValue();
    const VecCoord& tcp = targetContourPositions.getValue();

    unsigned int nbs=x.size(), nbt=tp.size(), nbtc = tcp.size(), nbsc = xcp.size(), nbs0=x0.size();

    distanceSet emptyset;
	
        if(nbs!=closestSource.size()) {if (!useVisible.getValue() || timer <= 2) initSource(); else initSourceVisible();  closestSource.resize(nbs);	closestSource.fill(emptyset); cacheDist.resize(nbs); cacheDist.fill((Real)0.); cacheDist2.resize(nbs); cacheDist2.fill((Real)0.); previousX.assign(x.begin(),x.end());}

        if(nbt!=closestTarget.size()) {initTarget();  initTargetContour(); closestTarget.resize(nbt);	closestTarget.fill(emptyset);}

    indicesTarget.resize(0);
		
        // closest target points from source points
        if(blendingFactor.getValue()<1)
        {
		
            /*double distmean = 0;
            vector<double> distm;
            distm.resize(0);
            double stddist = 0;*/

		for(int i=0;i<(int)nbt;i++)
                {
                    //int id = indicesVisible[i];
                    if(targetBorder[i])// && t%niterations.getValue() == 0)
                    {
                        double distmin = 10;
                        double dist;
                        int kmin;
                            for (int k = 0; k < x0.size(); k++)
                            {
                                if (sourceBorder[k])
                                {
                                    dist = (tp[i][0] - x0[k][0])*(tp[i][0] - x0[k][0]) + (tp[i][1] - x0[k][1])*(tp[i][1] - x0[k][1]) + (tp[i][2] - x0[k][2])*(tp[i][2] - x0[k][2]);
					if (dist < distmin)
					{
						distmin = dist;
						kmin = k;
					}
                                }
                            }
                        //distmean += distmin;
                        //distm.push_back(distmin);
                        indicesTarget.push_back(kmin);
                    }
        }
		
		/*distmean /= (double)distm.size();
		
        for(int i=0;i<(int)distm.size();i++)
        {
            stddist += (distm[i] - distmean);
        }*/

    //unsigned int count=0;
#ifdef USING_OMP_PRAGMAS
        #pragma omp parallel for
#endif
        for(int i=0;i<(int)nbs;i++)
        {	
            //if(sourceVisible[i])
            targetKdTree.getNClosest(closestSource[i],x[i],targetPositions.getValue(),1);
        }
    //std::cout<<(Real)count*(Real)100./(Real)nbs<<" % cached"<<std::endl;
    }	
		
    indices.resize(0);
		
#ifdef USING_OMP_PRAGMAS
        #pragma omp parallel for
#endif
        //for(int i=0;i<(int)nbs;i++)
		
    int kc = 0;
        for(int i=0;i<(int)nbs0;i++)
        {
            if(sourceBorder[i])// && t%niterations.getValue() == 0)
            {

                double distmin = 1000;
                double distmin1;
                double dist, dist1,dist2;
                int kmin2,kmin1;
				
                    for (int k = 0; k < tcp.size(); k++)
                    {
                        dist = (x0[i][0] - tcp[k][0])*(x0[i][0] - tcp[k][0]) + (x0[i][1] - tcp[k][1])*(x0[i][1] - tcp[k][1]) + (x0[i][2] - tcp[k][2])*(x0[i][2] - tcp[k][2]);

                        //if (dist < distmin)
                            if (dist < distmin)
                            {
                                distmin = dist;
                                kmin1 = k;
                            }
                    }
                double x_u_1 = ((x0[i][0])*rgbIntrinsicMatrix(0,0)/x0[i][2] + rgbIntrinsicMatrix(0,2)) - ((tcp[kmin1][0])*rgbIntrinsicMatrix(0,0)/tcp[kmin1][2] + rgbIntrinsicMatrix(0,2));
                double x_v_1 = ((x0[i][1])*rgbIntrinsicMatrix(1,1)/x0[i][2] + rgbIntrinsicMatrix(1,2)) - ((tcp[kmin1][1])*rgbIntrinsicMatrix(1,1)/tcp[kmin1][2] + rgbIntrinsicMatrix(1,2));
                double distmin0 = distmin;
                distmin = 1000;
                    for (int k = 0; k < tcp.size(); k++)
                    {
                        dist = (x0[i][0] - tcp[k][0])*(x0[i][0] - tcp[k][0]) + (x0[i][1] - tcp[k][1])*(x0[i][1] - tcp[k][1]) + (x0[i][2] - tcp[k][2])*(x0[i][2] - tcp[k][2]);
                        double x_u_2 = ((x0[i][0])*rgbIntrinsicMatrix(0,0)/x0[i][2] + rgbIntrinsicMatrix(0,2)) - ((tcp[k][0])*rgbIntrinsicMatrix(0,0)/tcp[k][2] + rgbIntrinsicMatrix(0,2));
                        double x_v_2 = ((x0[i][1])*rgbIntrinsicMatrix(1,1)/x0[i][2] + rgbIntrinsicMatrix(1,2)) - ((tcp[k][1])*rgbIntrinsicMatrix(1,1)/tcp[k][2] + rgbIntrinsicMatrix(1,2));

                        dist2 = abs(normalsContour[kc].y*x_u_2 - normalsContour[kc].x*x_v_2);
                        dist1 = x_u_2*x_u_1 + x_v_2*x_v_1;

                            //if (dist < distmin)
                            if (dist2 < distmin && sqrt(dist) < 0.10 && dist1 > 0 && sqrt(dist)/sqrt(distmin0)< 5)
                            {
                                distmin = dist2;
                                kmin2 = k;
                                distmin1 = dist1;
                            }
                    }
				
                        if (useDistContourNormal.getValue())
                            indices.push_back(kmin2);
                        else indices.push_back(kmin1);
                    kc++;
			
                }
            }
        std::cout << " indices size " << indices.size() << " tcp size " << tcp.size() << " xcp size " << xcp.size() << std::endl;

		
        // closest source points from target points
        if(blendingFactor.getValue()>0)
        {
            if (!useVisible.getValue() || timer <= 2)
                initSource();
            else initSourceVisible();
    #ifdef USING_OMP_PRAGMAS
            #pragma omp parallel for
    #endif
            for(int i=0;i<(int)nbt;i++)
                {
                    //if(!targetBackground[i])
                        sourceKdTree.getNClosest(closestTarget[i],tp[i],sourcePositions.getValue(),1);

                }
        }
	
    this->sourceIgnored.resize(nbs); sourceIgnored.fill(false);
    this->targetIgnored.resize(nbt); targetIgnored.fill(false);

    // prune outliers
    if(outlierThreshold.getValue()!=0) {
        Real mean=0,stdev=0,count=0;
        for(unsigned int i=0;i<nbs;i++) if(closestSource[i].size()) {count++; stdev+=closestSource[i].begin()->first; 
		mean+=(Real)(closestSource[i].begin()->first); 
		}
        for(unsigned int i=0;i<nbt;i++) if(closestTarget[i].size()) {count++; stdev+=closestTarget[i].begin()->first; mean+=(Real)(closestTarget[i].begin()->first); 
		//std::cout << " distances " << (double)(closestTarget[i].begin()->first) << std::endl;
		}
        mean=mean/count; stdev=(Real)sqrt(stdev/count-mean*mean);
        mean+=stdev*outlierThreshold.getValue();
        mean*=mean;
		int kkk=0;
        for(unsigned int i=0;i<nbs;i++) {
			if(closestSource[i].size()) 
			if(closestSource[i].begin()->first>mean ) 
				sourceIgnored[i]=true;
				
				if(sourceBorder[i])
				{
					
					double dists = (x[i][0] - tcp[indices[kkk]][0])*(x[i][0] - tcp[indices[kkk]][0]) + (x[i][1] - tcp[indices[kkk]][1])*(x[i][1] - tcp[indices[kkk]][1]) + (x[i][2] - tcp[indices[kkk]][2])*(x[i][2] - tcp[indices[kkk]][2]);

					if (sqrt(dists) > mean)
					sourceIgnored[i]=true;
					
					kkk++;

				}
				
		}
        for(unsigned int i=0;i<nbt;i++) if(closestTarget[i].size()) if(closestTarget[i].begin()->first>mean ) targetIgnored[i]=true;
		
		for(unsigned int i=0;i<nbs;i++) 
			if(closestSource[i].size()) 
				for(unsigned int j=0;j<nbt;j++)  
				if(j == closestSource[i].begin()->second && i == closestTarget[j].begin()->second) 
				{
					//sourceIgnored[i]=true;
					//targetIgnored[j] = true;
				}

    }
    if(rejectBorders.getValue()) {
        for(unsigned int i=0;i<nbs;i++) if(closestSource[i].size()) if(targetBorder[closestSource[i].begin()->second]) sourceIgnored[i]=true;
        for(unsigned int i=0;i<nbt;i++) if(closestTarget[i].size()) if(sourceBorder[closestTarget[i].begin()->second]) targetIgnored[i]=true;
    }
}

template<class DataTypes>
void ClosestPoint<DataTypes>::updateClosestPointsContoursNormals()
{
    const VecCoord& x = sourcePositions.getValue();
    const VecCoord& tp = targetPositions.getValue();
    const VecCoord& xcp = sourceContourPositions.getValue();
    const VecCoord& tcp = targetContourPositions.getValue();
    const VecCoord& ssn = sourceSurfaceNormalsM.getValue();

    unsigned int nbs=x.size(), nbt=tp.size(), nbtc = tcp.size(), nssn = ssn.size();

    distanceSet emptyset;
    if(nbs!=closestSource.size()) {initSource();  closestSource.resize(nbs);	closestSource.fill(emptyset); cacheDist.resize(nbs); cacheDist.fill((Real)0.); cacheDist2.resize(nbs); cacheDist2.fill((Real)0.); previousX.assign(x.begin(),x.end());}

	if(nbt!=closestTarget.size()) {initTarget();  /*initTargetContour();*/ closestTarget.resize(nbt);	closestTarget.fill(emptyset);}
					//std::cout << " tcp size () " << tcp.size() << std::endl;

    //if(nbt!=closestTarget.size()) {extractTargetPCD() ; closestTarget.resize(nbt);	closestTarget.fill(emptyset);}


    //if(nbs==0 || nbt==0) return;
	
	indices.resize(0);

    // closest target points from source points
    if(blendingFactor.getValue()<1) {

    //unsigned int count=0;
#ifdef USING_OMP_PRAGMAS
        #pragma omp parallel for
#endif
        for(int i=0;i<(int)nbs;i++)
        {
				if(sourceBorder[i])
				{

				//targetContourKdTree.getNClosest(closestSource[i],x[i],1);
				double distmin = 1000;
				double dist,dist0;
				double sign;
				int kmin = -1;
				//std::cout << " tcp size () " << tcp.size() << std::endl;
				for (int k = 0; k < tcp.size(); k++)
				{
					//dist = (ssn[i][0] - tcp[k][0])*(ssn[i][0] - tcp[k][0]) + (ssn[i][1] - tcp[k][1])*(ssn[i][1] - tcp[k][1]) +(ssn[i][2] - tcp[k][2])*(ssn[i][2] - tcp[k][2]);
					double norm = sqrt(ssn[i][0]*ssn[i][0] + ssn[i][1]*ssn[i][1] + ssn[i][2]*ssn[i][2]);
					double dist_ = sqrt((ssn[i][0]*(tcp[k][1]-x[i][1])-ssn[i][1]*(tcp[k][0]-x[i][0]))*(ssn[i][0]*(tcp[k][1]-x[i][1])-ssn[i][1]*(tcp[k][0]-x[i][0]))+(ssn[i][2]*(tcp[k][0]-x[i][0])-ssn[i][0]*(tcp[k][2]-x[i][2]))*(ssn[i][2]*(tcp[k][0]-x[i][0])-ssn[i][0]*(tcp[k][2]-x[i][2]))+(ssn[i][1]*(tcp[k][2]-x[i][2])-ssn[i][2]*(tcp[k][1]-x[i][1]))*(ssn[i][1]*(tcp[k][2]-x[i][2])-ssn[i][2]*(tcp[k][1]-x[i][1])));
					dist = dist_/norm;
					sign = (x[i][0] - tcp[k][0])*(ssn[k][0]) + (x[i][1] - tcp[k][1])*(ssn[k][1]) + (x[i][2] - tcp[k][2])*(ssn[k][2]);
					dist0 = (x[i][0] - tcp[k][0])*(x[i][0] - tcp[k][0]) + (x[i][1] - tcp[k][1])*(x[i][1] - tcp[k][1]) + (x[i][2] - tcp[k][2])*(x[i][2] - tcp[k][2]);
					//dist = (x[i][0] - tcp[k][0])*(x[i][0] - tcp[k][0]) + (x[i][1] - tcp[k][1])*(x[i][1] - tcp[k][1]) + (x[i][2] - tcp[k][2])*(x[i][2] - tcp[k][2]);
					//std::cout << " tcp size () " << sqrt(dist0) << std::endl;
					if (sqrt(dist0) < distmin)// && sqrt(dist0) < 0.20)
					{
						distmin = sqrt(dist0);
						kmin = k;
					}
				}
				indices.push_back(kmin);								
				unsigned int id=closestSource[i].begin()->second;
				int id1 = indices[i];
				//targetKdTree.getNClosest(closestSource[i],x[i],1);
			}
			
			targetKdTree.getNClosest(closestSource[i],x[i],targetPositions.getValue(),1);
			
        }
    //std::cout<<(Real)count*(Real)100./(Real)nbs<<" % cached"<<std::endl;
    }

    // closest source points from target points
    if(blendingFactor.getValue()>0)
    {
        initSource();
#ifdef USING_OMP_PRAGMAS
        #pragma omp parallel for
#endif
        for(int i=0;i<(int)nbt;i++)
		{
			{
            sourceKdTree.getNClosest(closestTarget[i],tp[i],sourcePositions.getValue(),1);
			}
			
		}
    }

    this->sourceIgnored.resize(nbs); sourceIgnored.fill(false);
    this->targetIgnored.resize(nbt); targetIgnored.fill(false);

    // prune outliers
    if(outlierThreshold.getValue()!=0) {
        Real mean=0,stdev=0,count=0;
        for(unsigned int i=0;i<nbs;i++) if(closestSource[i].size()) {count++; stdev+=closestSource[i].begin()->first; 
		mean+=(Real)(closestSource[i].begin()->first); 
		}
        for(unsigned int i=0;i<nbt;i++) if(closestTarget[i].size()) {count++; stdev+=closestTarget[i].begin()->first; mean+=(Real)(closestTarget[i].begin()->first); 
		//std::cout << " distances " << (double)(closestTarget[i].begin()->first) << std::endl;
		}
        mean=mean/count; stdev=(Real)sqrt(stdev/count-mean*mean);
        mean+=stdev*outlierThreshold.getValue();
        mean*=mean;
        for(unsigned int i=0;i<nbs;i++) if(closestSource[i].size()) if(closestSource[i].begin()->first>mean ) sourceIgnored[i]=true;
        for(unsigned int i=0;i<nbt;i++) if(closestTarget[i].size()) if(closestTarget[i].begin()->first>mean ) targetIgnored[i]=true;
		
		for(unsigned int i=0;i<nbs;i++) 
			if(closestSource[i].size()) 
				for(unsigned int j=0;j<nbt;j++)  
				if(j == closestSource[i].begin()->second && i == closestTarget[j].begin()->second) 
				{
					//sourceIgnored[i]=true;
					//targetIgnored[j] = true;
				}

    }
    if(rejectBorders.getValue()) {
        for(unsigned int i=0;i<nbs;i++) if(closestSource[i].size()) if(targetBorder[closestSource[i].begin()->second]) sourceIgnored[i]=true;
        for(unsigned int i=0;i<nbt;i++) if(closestTarget[i].size()) if(sourceBorder[closestTarget[i].begin()->second]) targetIgnored[i]=true;
    }
    /*if(normalThreshold.getValue()>(Real)-1. && sourceNormals.getValue().size()!=0 && targetNormals.getValue().size()!=0) {
        ReadAccessor< Data< VecCoord > > sn(sourceNormals);
        ReadAccessor< Data< VecCoord > > tn(targetNormals);
        for(unsigned int i=0;i<nbs;i++) if(closestSource[i].size()) if(dot(sn[i],tn[closestSource[i].begin()->second])<normalThreshold.getValue()) sourceIgnored[i]=true;
        for(unsigned int i=0;i<nbt;i++) if(closestTarget[i].size()) if(dot(tn[i],sn[closestTarget[i].begin()->second])<normalThreshold.getValue()) targetIgnored[i]=true;
    }*/
}
            

}
}
} // namespace sofa




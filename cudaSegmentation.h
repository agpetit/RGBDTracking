/*
 * cudaSegmentation.h
 *
 *  Created on: Apr 24, 2014
 *      Author: apetit
 */

/*
 * Copyright 1993-2014 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */

#ifndef _CUDASEGMENTATION_H_
#define _CUDASEGMENTATION_H_

#include <cuda_runtime.h>
#include <nppi.h>

class cudaSegmentation
{

    public:
	cudaSegmentation(const uchar4 *image, int image_pitch, unsigned char *trimap, int trimap_pitch, int width, int height);
        ~cudaSegmentation();

        void computeSegmentationFromTrimap();
        void computeSegmentation();
        void updateSegmentation();
        void updateSegmentationCrop();

        void updateImage(const uchar4 *image);
	void updateImageCrop(const uchar4 *crop_image, int _crop_image_pitch, int crop_width, int crop_height);
        void updateTrimap(uchar *trimap);
        void updateTrimapCrop(uchar *trimap, uchar *crop_trimap, int _crop_trimap_pitch);


        void setNeighborhood(int n)																												
        {
            m_neighborhood = n;
            computeSegmentationFromTrimap();
        };

        const unsigned char *getAlpha() const
        {
            return d_alpha[current_alpha];
        }
        int getAlphaPitch() const
        {
            return (int)alpha_pitch;
        };

        const unsigned char *getAlphaCrop() const
        {
            return d_crop_alpha[current_alpha];
        }
        int getAlphaPitchCrop() const
        {
            return (int)crop_alpha_pitch;
        };

    private:
        void createSmallImage(int max_dim);
        void createSmallTrimap();

        uchar4 *d_image;
        size_t image_pitch;
        float edge_strength;

        uchar4 *d_small_image;
        size_t small_pitch;
        NppiSize small_size;

        uchar4 *d_crop_image;
        size_t crop_image_pitch;

        const unsigned char *d_trimap;
        int trimap_pitch;

        unsigned char *d_small_trimap[2];
        size_t small_trimap_pitch[2];
        int small_trimap_idx;

        unsigned char *d_crop_trimap;
        size_t crop_trimap_pitch;

        NppiSize size;
        NppiSize crop_size;

        Npp32s *d_terminals;
        Npp32s *d_left_transposed, *d_right_transposed;
        Npp32s *d_top, *d_bottom, *d_topleft, *d_topright, *d_bottomleft, *d_bottomright;
        size_t pitch, transposed_pitch;
        int m_neighborhood;

        Npp32s *d_crop_terminals;
        Npp32s *d_crop_left_transposed, *d_crop_right_transposed;
        Npp32s *d_crop_top, *d_crop_bottom, *d_crop_topleft, *d_crop_topright, *d_crop_bottomleft, *d_crop_bottomright;
        size_t crop_pitch, crop_transposed_pitch;

        unsigned char *d_alpha[2];
        size_t alpha_pitch;
        int current_alpha;

        unsigned char *d_crop_alpha[2];
        size_t crop_alpha_pitch;
        int crop_current_alpha;

        Npp8u *d_scratch_mem;
        NppiGraphcutState *pState;
        
	Npp8u *d_crop_scratch_mem;
        NppiGraphcutState *crop_pState;

        float *d_gmm;
        size_t gmm_pitch;

        int *d_histogram;

        int gmms;
        int blocks;
        int crop_blocks;

        cudaEvent_t start, stop;
};

#endif //_CUDASEGMENTATION_H_

/*=========================================================================
 
 Program:   Insight Segmentation & Registration Toolkit
 Module:    $RCSfile: itkCreateRGBProjectionsFromDWI.h,v $
 Language:  C++
 Date:      $Date: 2006/03/27 17:01:10 $
 Version:   $Revision: 1.15 $
 
 Copyright (c) Insight Software Consortium. All rights reserved.
 See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.
 
 This software is distributed WITHOUT ANY WARRANTY; without even 
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
 PURPOSE.  See the above copyright notices for more information.
 
 =========================================================================*/
#ifndef __itkCreateRGBProjectionsFromDWI_h
#define __itkCreateRGBProjectionsFromDWI_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkVector.h"
#include "itkArray.h"
#include <vector>

namespace itk
{
    
    
    template <class TInputImage, class TOutputImage>
    class ITK_EXPORT CreateRGBProjectionsFromDWI : public ImageToImageFilter< TInputImage, TOutputImage >
    {
    public:
        /** Standard class typedefs. */
        typedef CreateRGBProjectionsFromDWI           Self;
        /** Convenient typedefs for simplifying declarations. */
        typedef TInputImage                           InputImageType;
        typedef typename InputImageType::Pointer      InputImagePointer;
        typedef typename InputImageType::ConstPointer InputImageConstPointer;
        typedef TOutputImage                          OutputImageType;
        typedef typename OutputImageType::Pointer     OutputImagePointer;
        
        /** Standard class typedefs. */
        typedef ImageToImageFilter< InputImageType, OutputImageType> Superclass;
        typedef SmartPointer<Self>                                   Pointer;
        typedef SmartPointer<const Self>                             ConstPointer;
        
        /** Method for creation through the object factory. */
        itkNewMacro(Self);
        
        /** Run-time type information (and related methods). */
        itkTypeMacro( CreateRGBProjectionsFromDWI, ImageToImageFilter );
        
        /** Image typedef support. */
        typedef typename InputImageType::PixelType           InputPixelType;
        typedef typename OutputImageType::PixelType          OutputPixelType;
        typedef typename InputImageType::RegionType          InputImageRegionType;
        typedef typename InputImageType::SizeType            InputImageSizeType;
        typedef typename InputImageType::IndexType           InputImageIndexType;
        typedef typename OutputImageType::RegionType         OutputImageRegionType;
        
        /** Indicator type (which channels are actual gradients -as opposed to baselines-)?: */
        typedef itk::Array<unsigned int> IndicatorType;
        void SetDWI( const IndicatorType idx )
        {
            m_DWI = idx;
        }
        
        /** Type to store the weights of each channel */
        typedef std::vector<float> WeightsType;
        
        /** Types to store the gradients table */
        typedef itk::Vector<double, 3>    GradientType;
        typedef std::vector<GradientType> GradientListType;
        /** Set the gradients list: */
        void SetGradientsTable( const GradientListType grads )
        {
            m_GradientsTable = grads;
        }
        
        /** Method to set the dimension to operate along:*/
        void SetProjectionCoordinate( const unsigned int d )
        {
            m_ProjectionCoordinate = d;
        }
        
        /** Method to compute the amount of noise in the projection
         as a function of the amount of nosie in the original DWI
         volume
         // NOTE: if the weights are w_i, we estimate the variance is
         // sum_i w_i^2 sigma^2, with sigma the original std of noise in
         // the complex domain of the x-space. This is not strictly valid
         // since we have Rician (not Gaussian) noise and the std is a 
         // function of the mean. This approximation is only valid for
         // high SNR.
         //
         // HOWEVER: in the large SNR limit (A>>sigma), the noise is 
         // almost Gaussian and sigma_Rician \simeq sigam_Gaussian. In
         // the low SNR case (A<<sigma), the noise is almost Rayleigh
         // and sigma_Rician \simeq sqrt((4-pi)/2)*sigma_Gaussian
         // < Sigma_Gaussian. I.e: for low SNR we over-estimate the
         // noise, hence the filtering is more agressive. For high
         // SNR the estimation is accurate, hence we preserve the
         // details. This makes a lot of sense, by the way...
         */
        float GetCorrectionFactor( void ) const
        {
            return m_CorrectionFactor;
        }
        
    protected:
        CreateRGBProjectionsFromDWI();
        virtual ~CreateRGBProjectionsFromDWI() {}
#if ITK_VERSION_MAJOR < 4
        void ThreadedGenerateData( const OutputImageRegionType & outputRegionForThread, int threadId );
        
#else
        void ThreadedGenerateData( const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId ) ITK_OVERRIDE;
        
#endif
        void BeforeThreadedGenerateData( void ) ITK_OVERRIDE;
    private:
        CreateRGBProjectionsFromDWI(const Self&); // purposely not implemented
        void operator=(const Self&);              // purposely not implemented
        /** The actual list of channels corresponding to gradients*/
        IndicatorType    m_DWI;
        /** The table with the gradients: */
        GradientListType m_GradientsTable;
        /** The list of weights: */
        WeightsType      m_Weights;
        /** The dimension to operate along: */
        unsigned int     m_ProjectionCoordinate;
        /** The resulting std of noise: */
        float            m_CorrectionFactor;
    };
    
    
    
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkCreateRGBProjectionsFromDWI.txx"
#endif

#endif

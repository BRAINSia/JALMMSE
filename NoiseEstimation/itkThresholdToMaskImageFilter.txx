/*=========================================================================
 
 Program:   Insight Segmentation & Registration Toolkit
 Module:    $RCSfile: itkThresholdToMaskImageFilter.txx,v $
 Language:  C++
 Date:      $Date: 2005/05/4 14:28:51 $
 Version:   $Revision: 1.1
 =========================================================================*/
#ifndef _itkThresholdToMaskImageFilter_txx
#define _itkThresholdToMaskImageFilter_txx

#include "itkThresholdToMaskImageFilter.h"
#include "itkImageRegionIterator.h"

namespace itk
{
    
    /** Constructor */
    template <class TInputImage, class TOutputImage>
    ThresholdToMaskImageFilter<TInputImage, TOutputImage>::ThresholdToMaskImageFilter()
    {
        m_Threshold = itk::NumericTraits<InputPixelType>::Zero;
    }
    
   
    
    template <class TInputImage, class TOutputImage>
    void ThresholdToMaskImageFilter<TInputImage, TOutputImage>
    ::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread,
                           ThreadIdType itkNotUsed(threadId) )
    {
        // Input and output
        InputImageConstPointer input   =  this->GetInput();
        OutputImagePointer     output  =  this->GetOutput();
        // Iterators:
        ImageRegionConstIterator<InputImageType> bit( input, outputRegionForThread );
        ImageRegionIterator<OutputImageType>     it( output, outputRegionForThread );
        
        for( bit.GoToBegin(), it.GoToBegin(); !bit.IsAtEnd(); ++bit, ++it ){
            OutputPixelType val = static_cast<OutputPixelType>( bit.Get()>m_Threshold ? 1 : 0 );
            it.Set( val );
        }
    }
    
} // end namespace itk

#endif

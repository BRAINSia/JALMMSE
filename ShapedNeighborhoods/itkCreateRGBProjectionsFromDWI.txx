/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkitkCreateRGBProjectionsFromDWI.txx,v $
  Language:  C++
  Date:      $Date: 2006/01/11 19:43:31 $
  Version:   $Revision: 1.21 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkCreateRGBProjectionsFromDWI_txx
#define _itkCreateRGBProjectionsFromDWI_txx

#include "itkCreateRGBProjectionsFromDWI.h"
#include "itkImageRegionIterator.h"
#include "math.h"

namespace itk
{
template< class TInputImage, class TOutputImage >
CreateRGBProjectionsFromDWI< TInputImage, TOutputImage >
::CreateRGBProjectionsFromDWI()
{
    m_DWI.SetSize(0);
    m_GradientsTable.clear();
    m_Weights.clear();
    m_ProjectionCoordinate = 0;
    m_CorrectionFactor = 1;
}

template< class TInputImage, class TOutputImage >
void CreateRGBProjectionsFromDWI< TInputImage, TOutputImage >
::BeforeThreadedGenerateData( void )
{
    Superclass::BeforeThreadedGenerateData();
    // Perform sanity checks:
    if( !this->GetInput() )
        itkExceptionMacro( << "You must set the input before updating the filter" );
    if( m_GradientsTable.size()==0 )
        itkExceptionMacro( << "You must set the gradients table before running the filter" );
    if( m_GradientsTable.size() != m_DWI.GetSize() )
        itkExceptionMacro( << "The number of DWI channels must match the number of entries in the gradients table" );
    if( m_ProjectionCoordinate >= TInputImage::ImageDimension )
        itkExceptionMacro( << "The projection coordinate exceeds the maximum dimension of the image" );
    if( m_DWI.GetSize() > this->GetInput()->GetVectorLength() )
        itkExceptionMacro( << "The number of gradients is larger than the actual number of channels in the volume!!!" );

    /** In case all these things succeeded, compute the weights:*/
    double norm = itk::NumericTraits<double>::Zero; // For normalization purposes
    m_Weights.resize( m_DWI.GetSize() );
    for( unsigned int k=0; k<m_DWI.GetSize(); ++k ){
        m_Weights[k] = std::abs( m_GradientsTable[k][ m_ProjectionCoordinate ] );
        norm += m_Weights[k];
    }
    /** Normalize the weights so that they sum up to 1 */
    for( unsigned int k=0; k<m_DWI.GetSize(); ++k )
        m_Weights[k] /= norm;

    /** Update the correction factor to apply to the original std of noise in the
    DWI image to compute the std of noise in the projected image:
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
    m_CorrectionFactor = itk::NumericTraits<float>::Zero;
    for( unsigned int k=0; k<m_DWI.GetSize(); ++k )
        m_CorrectionFactor += ( m_Weights[k] * m_Weights[k] );
    m_CorrectionFactor = std::sqrt(m_CorrectionFactor);
    // Ready to go!
}


template< class TInputImage, class TOutputImage >
#if ITK_VERSION_MAJOR < 4
void CreateRGBProjectionsFromDWI< TInputImage, TOutputImage >
::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int itkNotUsed(threadId) )
#else
void CreateRGBProjectionsFromDWI< TInputImage, TOutputImage >
::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, ThreadIdType itkNotUsed(threadId) )
#endif
{
   // Input and output
   InputImageConstPointer   input   =  this->GetInput();
   OutputImagePointer       output  =  this->GetOutput();

   ImageRegionConstIterator<InputImageType>        bit = ImageRegionConstIterator<InputImageType>( input,  outputRegionForThread );
   ImageRegionIterator<OutputImageType>            it  = ImageRegionIterator<OutputImageType>(     output, outputRegionForThread );

   for( bit.GoToBegin(),it.GoToBegin(); !bit.IsAtEnd(); ++bit,++it )
   {
      const InputPixelType & inputPixel = bit.Get();
      double outputPixel = itk::NumericTraits<double>::Zero;
      for( unsigned int j=0; j<m_DWI.GetSize(); ++j )
        {
        outputPixel += ( m_Weights[j] ) * ( inputPixel[ m_DWI[j] ] );
        }
      it.Set( static_cast<OutputPixelType>(outputPixel) );
   }
}

} // end namespace itk

#endif

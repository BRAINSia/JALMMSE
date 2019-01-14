/*=========================================================================
 
 Program:   Insight Segmentation & Registration Toolkit
 Module:    $RCSfile: itkComputeRestrictedHistogram.txx,v $
 Language:  C++
 Date:      $Date: 2006/01/11 19:43:31 $
 Version:   $Revision: 1.14 $
 
 Copyright (c) Insight Software Consortium. All rights reserved.
 See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.
 
 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.
 
 =========================================================================*/
#ifndef _itkComputeRestrictedHistogram_txx
#define _itkComputeRestrictedHistogram_txx
#include "itkComputeRestrictedHistogram.h"

#include "itkImageRegionConstIterator.h"
#include "itkNumericTraits.h"

namespace itk
{
    
    template <class TInputImage, class TOutputImage>
    ComputeRestrictedHistogram<TInputImage, TOutputImage>
    ::ComputeRestrictedHistogram()
    {
        m_Ready = false;
        m_Min   = 0.0f;
        m_Max   = 5.0f;
        m_Bins  = 1024;
        m_Mode  = 0.0f;
    }
    
    template <class TInputImage, class TOutputImage>
    void ComputeRestrictedHistogram<TInputImage, TOutputImage>
    ::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, ThreadIdType threadId )
    {
        // Allocate input:
        typename InputImageType::ConstPointer input  = this->GetInput();
        // Iterators:
        itk::ImageRegionConstIterator<InputImageType> bit = itk::ImageRegionConstIterator<InputImageType>(
                                                                                                          input, outputRegionForThread );
        // The scale factor to compute histogram position:
        double scale = ( (double)m_Bins ) / ( m_Max - m_Min );
        // Initiallise max and min:
        for( bit.GoToBegin(); !bit.IsAtEnd(); ++bit ){
            for( unsigned int p=0; p<this->GetInput()->GetVectorLength(); ++p ){
                double current = bit.Get()[p];
                if(   ( current >= m_Min ) && ( current < m_Max)   ){
                    unsigned int pos = (unsigned int)(   ( current - m_Min ) * scale   );
                    m_THist[threadId][pos]++;
                }
            }
        }
    }
    
    template <class TInputImage, class TOutputImage>
    void ComputeRestrictedHistogram<TInputImage, TOutputImage>
    ::BeforeThreadedGenerateData()
    {
        this->Superclass::BeforeThreadedGenerateData();
        if( m_Bins<128 )
            m_Bins = 128;
#if ITK_VERSION_MAJOR >= 5
        unsigned int k = this->GetNumberOfWorkUnits();
#else
        unsigned int k = this->GetNumberOfThreads();
#endif
        m_THist.SetSize( k, m_Bins );
        m_THist.Fill( 0 );
    }
    
    template <class TInputImage, class TOutputImage>
    void ComputeRestrictedHistogram<TInputImage, TOutputImage>
    ::AfterThreadedGenerateData()
    {
#if ITK_VERSION_MAJOR >= 5
        unsigned int  k = this->GetNumberOfWorkUnits();
#else
        unsigned int  k = this->GetNumberOfThreads();
#endif
        unsigned long max    = 0;
        unsigned int  maxPos = 0;
        
        // Gather all the partial hisotgrams to compute the final histogram:
        m_Histogram.SetSize( m_Bins );
        for( unsigned int i = 0; i < m_Bins; ++i ){
            m_Histogram[i] = m_THist[0][i];
            for( unsigned int j = 1; j < k; ++j )
                m_Histogram[i] += m_THist[j][i];
        }
        // Smooth the histogram:
        unsigned int R = (unsigned int)((double)m_Bins/64);
        if( R<1 )
            R = 1;
        double* window = new double[2*R+1]; // This mimicks matlab's gausswin
        double norm   = itk::NumericTraits<double>::Zero;
        for( int r=-((int)R); r<=((int)R); ++r ){
            window[(unsigned int)(r+(int)R)] = std::exp( -3.125 * (double)(r*r) / (double)(R*R)  );
            norm += window[(unsigned int)(r+(int)R)];
        }
        for( unsigned int r=0; r<2*R+1; ++r )
            window[r] /= norm;
        for( unsigned int i=R; i < (unsigned int)((int)m_Bins-(int)R); ++i ){
            m_Histogram[i] *= window[R];
            for( unsigned int r=1; r<=R; ++r ){
                m_Histogram[i] += (m_Histogram[i+r]*window[R+r]);
                m_Histogram[i] += (m_Histogram[i-r]*window[R-r]);
            }
        }
        delete[] window;
        // Find the maximum of the histogram, i.e. the arg mode:
        for( unsigned int i = 0; i < m_Bins; ++i ){
            if( m_Histogram[i] > max ){
                max    = m_Histogram[i];
                maxPos = i;
            }
        }
        m_Mode  = m_Min + ( (double)maxPos + 0.5f ) * ( m_Max - m_Min ) / ( (double)m_Bins );
        m_Ready = true;
        this->Superclass::AfterThreadedGenerateData();
    }
    
} // end namespace itk

#endif

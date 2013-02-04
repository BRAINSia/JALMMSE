/*=========================================================================
 
 Program:   Insight Segmentation & Registration Toolkit
 Module:    $RCSfile: itkVectorLocalStdImageFilter.txx,v $
 Language:  C++
 Date:      $Date: 2005/05/4 14:28:51 $
 Version:   $Revision: 1.1
 =========================================================================*/
#ifndef _itkVectorLocalStdImageFilter_txx
#define _itkVectorLocalStdImageFilter_txx

#include "itkVectorLocalStdImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"

namespace itk
{
    
    /** Constructor */
    template <class TInputImage, class TOutputImage>
    VectorLocalStdImageFilter<TInputImage, TOutputImage>::VectorLocalStdImageFilter()
    {
        m_Mask = NULL;
        m_Radius.Fill(1);
        m_PerThreadSum.SetSize( this->GetNumberOfThreads() );
        m_PerThreadCount.SetSize( this->GetNumberOfThreads() );
        m_LocalStdMean = itk::NumericTraits<double>::Zero;
    }
    
    template <class TInputImage, class TOutputImage>
    void VectorLocalStdImageFilter<TInputImage, TOutputImage>
    ::BeforeThreadedGenerateData( void )
    {
        m_PerThreadSum.SetSize( this->GetNumberOfThreads() );
        m_PerThreadSum.Fill( itk::NumericTraits<double>::Zero );
        m_PerThreadCount.SetSize( this->GetNumberOfThreads() );
        m_PerThreadCount.Fill( itk::NumericTraits<unsigned long>::Zero );
        this->Superclass::BeforeThreadedGenerateData();
    }
    
    template <class TInputImage, class TOutputImage>
    void VectorLocalStdImageFilter<TInputImage, TOutputImage>
    ::AfterThreadedGenerateData( void )
    {
        double        cumsum = itk::NumericTraits<double>::Zero;
        unsigned long counts = 0;
        for( unsigned int k=0; k<(unsigned int)(this->GetNumberOfThreads()); ++k ){
            cumsum += m_PerThreadSum[k];
            counts += m_PerThreadCount[k];
        }
        cumsum         /= static_cast<double>(counts);
        m_LocalStdMean  = static_cast<double>(cumsum);
        this->Superclass::AfterThreadedGenerateData();
    }
    
    template <class TInputImage, class TOutputImage>
#if ITK_VERSION_MAJOR < 4
    void VectorLocalStdImageFilter<TInputImage, TOutputImage>
    ::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread,
                           int threadId )
#else
    void VectorLocalStdImageFilter<TInputImage, TOutputImage>
    ::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread,
                           ThreadIdType threadId )
    
#endif
    {
        // Input and output
        InputImageConstPointer input   =  this->GetInput();
        OutputImagePointer     output  =  this->GetOutput();
        // Iterators:
        itk::ConstNeighborhoodIterator<InputImageType> bit;
        itk::ImageRegionIterator<OutputImageType>      it;
        itk::ImageRegionConstIterator<MaskImageType>   mit;
        
        itk::ZeroFluxNeumannBoundaryCondition<InputImageType> nbc;
        
#if ITK_VERSION_MAJOR < 4
#else
        typedef typename itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>::FaceListType FaceListType;
        typedef typename FaceListType::iterator                                                                 FaceIteratorType;
        typedef typename itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>               FaceCalculatorType;
        
        FaceCalculatorType faceCalculator;
        FaceListType       faceList = faceCalculator( input, outputRegionForThread, m_Radius );
        FaceIteratorType   faceIterator;
#endif
        double        cumsum = itk::NumericTraits<double>::Zero;
        unsigned long counts = 0;
        
#if ITK_VERSION_MAJOR < 4
        /** In ITKv3 the corners of the outputRegionForThread are processed
         twice!!! :-(. This is the reason why I need this clause to ensure
         exactly the same behavior in ITKv3 and ITKv4 */
        for( unsigned int k=0; k<1; ++k ){
            OutputImageRegionType currentRegion = outputRegionForThread;
#else
        for( faceIterator=faceList.begin(); faceIterator!=faceList.end(); ++faceIterator ){
            OutputImageRegionType currentRegion = *faceIterator;
#endif
            bit = itk::ConstNeighborhoodIterator<InputImageType>( m_Radius, input,  currentRegion );
            it  = itk::ImageRegionIterator<OutputImageType>(                output, currentRegion );
            if( m_Mask ){
                mit = itk::ImageRegionConstIterator<MaskImageType>( m_Mask, currentRegion );
                mit.GoToBegin();
            }
            bit.OverrideBoundaryCondition( &nbc );
            
            const unsigned long N = bit.Size();
            
            for( bit.GoToBegin(), it.GoToBegin(); !bit.IsAtEnd(); ++bit, ++it ){
                
                OutputPixelType avg( this->GetInput()->GetVectorLength() );
                avg.Fill( itk::NumericTraits<OutputScalarType>::Zero );
                OutputPixelType sqavg( this->GetInput()->GetVectorLength() );
                sqavg.Fill( itk::NumericTraits<OutputScalarType>::Zero );
                
                if( m_Mask ){
                    bool doIt = static_cast<bool>(mit.Get());
                    ++mit;
                    if( !doIt ){
                        it.Set( sqavg );
                        continue;
                    }
                }
                
                for( unsigned int n=0; n<N; ++n ){
                    InputPixelType  ipix = bit.GetPixel(n);
                    for( unsigned int g=0; g<this->GetInput()->GetVectorLength(); ++g ){
                        OutputScalarType value = static_cast<OutputScalarType>(ipix[g]);
                        avg[g]   += value;
                        sqavg[g] += (value*value);
                    }
                }
                
                for( unsigned int g=0; g<this->GetInput()->GetVectorLength(); ++g ){
                    double value = itk::NumericTraits<double>::Zero;
                    double mean  = itk::NumericTraits<double>::Zero;
                    if( N>1 ){
                        mean   = static_cast<double>(avg[g]) / static_cast<double>(N);
                        value  = ( static_cast<double>(sqavg[g]) - static_cast<double>((avg[g])*(avg[g]))/static_cast<double>(N) );
                        value /= static_cast<double>(N-1);
                        if( value>itk::NumericTraits<double>::Zero )
                            value = vcl_sqrt( value );
                        value = this->ComputeStdCorrection( mean, value );
                    }
                    sqavg[g]  = static_cast<OutputScalarType>(value);
                    cumsum   += value;
                    counts++;
                }
                
                it.Set( sqavg );
                
            }
        }
        
        m_PerThreadCount[static_cast<unsigned int>(threadId)] = counts;
        m_PerThreadSum[static_cast<unsigned int>(threadId)]   = cumsum;
        
        return;
    }
    
    template <class TInputImage, class TOutput>
    double VectorLocalStdImageFilter<TInputImage, TOutput>
    ::ComputeStdCorrection( const double inavg, const double instd ) const
    {
        if( inavg > 2*instd )
            return instd;
        else if( instd > 2*inavg )
            return ( 1.5264f * instd );
        else{
            double sigma = instd/inavg;
            sigma = 0.5f + 2.5528*(sigma - 0.5)/( 2 - 0.5 );
            return( sigma*inavg );
        }
    }
    
} // end namespace itk

#endif

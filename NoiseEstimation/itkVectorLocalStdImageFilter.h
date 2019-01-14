/*=========================================================================
 
 Program:   Insight Segmentation & Registration Toolkit
 Module:    $RCSfile: itkVectorLocalStdImageFilter.h,v $
 Language:  C++
 Date:      $Date: 2008/02/7 14:28:51 $
 Version:   $Revision: 0.0 $
 =========================================================================*/
#ifndef __itkVectorLocalStdImageFilter_h
#define __itkVectorLocalStdImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkArray.h"

namespace itk
{
    
    template <class TInputImage, class TOutputImage>
    class ITK_EXPORT VectorLocalStdImageFilter : public ImageToImageFilter<TInputImage,TOutputImage>
    {
    public:
        /** Convenient typedefs for simplifying declarations. */
        typedef TInputImage                            InputImageType;
        typedef typename InputImageType::Pointer       InputImagePointer;
        typedef typename InputImageType::ConstPointer  InputImageConstPointer;
        typedef TOutputImage                           OutputImageType;
        typedef typename OutputImageType::Pointer      OutputImagePointer;
        typedef typename OutputImageType::ConstPointer OutputImageConstPointer;
        
        /** Standard class typedefs. */
        typedef VectorLocalStdImageFilter                           Self;
        typedef ImageToImageFilter<InputImageType, OutputImageType> Superclass;
        typedef SmartPointer<Self>                                  Pointer;
        typedef SmartPointer<const Self>                            ConstPointer;
        
        /** Method for creation through the object factory. */
        itkNewMacro(Self);
        
        /** Run-time type information (and related methods). */
        itkTypeMacro( VectorLocalStdImageFilter, ImageToImageFilter );
        
        /** Image typedef support. */
        typedef typename InputImageType::PixelType   InputPixelType;
        typedef typename InputPixelType::ValueType   InputScalarType;
        typedef typename InputImageType::SizeType    InputSizeType;
        typedef typename OutputImageType::PixelType  OutputPixelType;
        typedef typename OutputPixelType::ValueType  OutputScalarType;
        typedef typename OutputImageType::RegionType OutputImageRegionType;
        
        typedef unsigned short                                    LabelType;
        typedef itk::Image<LabelType,TInputImage::ImageDimension> MaskImageType;
        typedef typename MaskImageType::Pointer                   MaskImagePointer;
        
        typedef itk::Array<double>                   FloatSumType;
        typedef itk::Array<unsigned long>            LongSumType;
        
        itkSetMacro( Radius, InputSizeType );
        itkGetMacro( Radius, InputSizeType );
        
        double GetLocalStdMean( void ) const
        {
            return m_LocalStdMean;
        }
        
        void SetMask( const MaskImagePointer mask )
        {
            m_Mask = mask;
        }
        
    protected:
        VectorLocalStdImageFilter();
        virtual ~VectorLocalStdImageFilter(){}
        void ThreadedGenerateData( const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId ) ITK_OVERRIDE;
        
        void BeforeThreadedGenerateData(void) ITK_OVERRIDE;
        void AfterThreadedGenerateData(void) ITK_OVERRIDE;
        double ComputeStdCorrection( const double, const double ) const;
        
    private:
        VectorLocalStdImageFilter(const Self &); // purposely not implemented
        void operator=(const Self &);            // purposely not implemented
        
        InputSizeType    m_Radius;
        
        FloatSumType     m_PerThreadSum;
        LongSumType      m_PerThreadCount;
        
        double           m_LocalStdMean;
        
        MaskImagePointer m_Mask;
    };
    
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVectorLocalStdImageFilter.txx"
#endif

#endif

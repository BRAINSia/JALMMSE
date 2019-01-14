/*=========================================================================
 
 Program:   Insight Segmentation & Registration Toolkit
 Module:    $RCSfile: itkLMMSEVectorImageFilter.h,v $
 Language:  C++
 Date:      $Date: 2008/02/7 14:28:51 $
 Version:   $Revision: 0.0 $
 =========================================================================*/
#ifndef __itkLMMSEVectorImageFilter_h
#define __itkLMMSEVectorImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImage.h"
#include <vector>
#include "itkVector.h"
#include "itkArray.h"
#include "itkArray2D.h"
#include "itkFixedArray.h"

// Shaped neighborhoods includes
#include "ShapedNeighborhoods/itkLSDerivatives.h"
#include "ShapedNeighborhoods/itkCreateRGBProjectionsFromDWI.h"

namespace itk
{
    /* *************************************************************************** */
    // We use tis piece of code to sort the gradients with respect to a given one
    // in an ascending order with respect to the angular distance
    typedef itk::FixedArray<double, 2> OrderType;
    // To use with the sort method of std::vector
    bool UNLM_gradientDistance_smaller( OrderType e1, OrderType e2 )
    {
        return e1[1] < e2[1];
    }
    /* *************************************************************************** */
    
    /** \class UNLMFilter
     *
     * DO NOT assume a particular image or pixel type, which is, the input image
     * may be a VectorImage as well as an Image obeject with vectorial pixel type.
     *
     * \sa Image
     */
    
    template <class TInputImage, class TOutputImage>
    class ITK_EXPORT LMMSEVectorImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
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
        typedef LMMSEVectorImageFilter                              Self;
        typedef ImageToImageFilter<InputImageType, OutputImageType> Superclass;
        typedef SmartPointer<Self>                                  Pointer;
        typedef SmartPointer<const Self>                            ConstPointer;
        
        /** Method for creation through the object factory. */
        itkNewMacro(Self);
        
        /** Run-time type information (and related methods). */
        itkTypeMacro( LMMSEVectorImageFilter, ImageToImageFilter );
        
        /** Image typedef support. */
        typedef typename InputImageType::PixelType   InputPixelType;
        typedef typename InputImageType::SizeType    InputImageSizeType;
        typedef typename InputImageType::RegionType  InputImageRegionType;
        typedef typename InputImageType::IndexType   InputImageIndexType;
        typedef typename OutputImageType::PixelType  OutputPixelType;
        typedef typename OutputPixelType::ValueType  ScalarType;
        typedef typename OutputImageType::RegionType OutputImageRegionType;
        
        typedef typename InputImageType::SizeType InputSizeType;
        
        /** ############################################################################################################# */
        // THIS TYPE IS USED TO INCLUDE A MASK IN THE COMPUTATIONS. THE VALUES OUTSIDE THIS MASK ARE NOT FILTERED TO
        // SAVE UNNECESSARY COMPUTATIONS:
        typedef itk::Image<unsigned short,TInputImage::ImageDimension> LabelImageType;
        typedef typename LabelImageType::Pointer                       LabelImagePointer;
        typedef itk::ImageRegionConstIterator<LabelImageType>          MaskIteratorType;
        /** ############################################################################################################# */
        
        
        /** ############################################################################################################# */
        // THESE TYPES ARE USED TO COMPUTE THE SHAPED NEIGHBORHOODS:
        typedef itk::Image<float,TInputImage::ImageDimension>                   FloatImageType;
        typedef itk::CreateRGBProjectionsFromDWI<InputImageType,FloatImageType> RGBProjectionType;
        typedef typename RGBProjectionType::Pointer                             RGBProjectionPointer;
        
        typedef itk::Image<LSGradientsL2,TInputImage::ImageDimension >   FeaturesMapType;
        typedef typename FeaturesMapType::Pointer                        FeaturesMapPointer;
        
        typedef itk::LSDerivativesL0<FloatImageType>                     L0Type;
        typedef typename L0Type::Pointer                                 L0Pointer;
        typedef itk::LSDerivativesL1<TInputImage::ImageDimension>        L1Type;
        typedef typename L1Type::Pointer                                 L1Pointer;
        typedef itk::LSDerivativesL2<TInputImage::ImageDimension>        L2Type;
        typedef typename L2Type::Pointer                                 L2Pointer;
        /** ############################################################################################################# */
        
        
        
        /** Set and get the radius of the neighborhood used to compute the statistics. */
        itkSetMacro(               Radius, InputSizeType );
        itkGetConstReferenceMacro( Radius, InputSizeType );
        
        /** Set and get the radius of the neighborhood used to compute the salient features. */
        itkSetMacro(               RadiusFeatures, InputSizeType );
        itkGetConstReferenceMacro( RadiusFeatures, InputSizeType );
        
        
        // Type of the gradients direction:
        typedef itk::Vector<double, 3> GradientType;
        // Type to store all the gradient directions
        typedef std::vector<GradientType> GradientListType;
        // Indicator type:
        typedef itk::Array<unsigned int> IndicatorType;
        // The type to store the closest gradient directions to a given one.
        typedef itk::Array2D<unsigned int> NeighboursIndType;
        /** Set and get the parameters */
        itkSetMacro( NDWI,           unsigned int );
        itkGetMacro( NDWI,           unsigned int );
        itkSetMacro( NBaselines,     unsigned int );
        itkGetMacro( NBaselines,     unsigned int );
        itkSetMacro( Sigma,          float        );
        itkGetMacro( Sigma,          float        );
        itkSetMacro( H,              float        );
        itkGetMacro( H,              float        );
        itkSetMacro( Neighbours,     unsigned int );
        itkGetMacro( Neighbours,     unsigned int );
        itkSetMacro( SetZeroBck,     bool         );
        itkGetMacro( SetZeroBck,     bool         );
        itkSetMacro( OnlyUNLM,       bool         );
        itkGetMacro( OnlyUNLM,       bool         );
        itkSetMacro( FilterOutliers, bool         );
        itkGetMacro( FilterOutliers, bool         );
        /** Add a new gradient direction: */
        void AddGradientDirection( GradientType grad )
        {
            m_GradientList.push_back( grad );
            return;
        }
        
        /** Set the vector with the DWI channels that are going to be used: */
        void SetDWI( IndicatorType ind )
        {
            m_DWI = ind;
        }
        
        void SetBaselines( IndicatorType ind )
        {
            m_Baselines = ind;
        }
        
        /** Set the vector of DWI channels using c-style vector. The user must set
         *  m_Channels before */
        void SetDWI( unsigned int* ind )
        {
            m_DWI.SetSize( m_NDWI );
            for( unsigned int k = 0; k < m_NDWI; ++k )
            {
                m_DWI[k] = ind[k];
            }
        }
        
        void SetBaselines( unsigned int* ind )
        {
            m_Baselines.SetSize( m_NBaselines );
            for( unsigned int k = 0; k < m_NBaselines; ++k )
            {
                m_Baselines[k] = ind[k];
            }
        }
        
        IndicatorType GetDWI(void)
        {
            return m_DWI;
        }
        IndicatorType GetBaselines(void)
        {
            return m_Baselines;
        }
        
        const GradientListType GetGradientList( void )
        {
            return m_GradientList;
        }
        
        /** In case we use a mask, use this method to fix it */
        void SetMask( const LabelImagePointer& mask ){
            m_Mask = mask;
        }
    protected:
        LMMSEVectorImageFilter();
        virtual ~LMMSEVectorImageFilter()
        {
        }
        void PrintSelf(std::ostream& os, Indent indent) const ITK_OVERRIDE;
        
#if ITK_VERSION_MAJOR < 4
        void ThreadedGenerateData( const OutputImageRegionType & outputRegionForThread, int threadId );
        
#else
        void ThreadedGenerateData( const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId ) ITK_OVERRIDE;
        
#endif
        void BeforeThreadedGenerateData( void ) ITK_OVERRIDE;
        
        virtual void GenerateInputRequestedRegion() throw (InvalidRequestedRegionError) ITK_OVERRIDE;
        
        bool ComputeInverseMatrix( const double *, const double *, double, double *, unsigned int ) const;
        
        void CMMInversion( const double *, const double *, double, double *, unsigned int, unsigned int) const;
        
        float ComputeTraceMO1( const InputImageSizeType& rcomp );
        
        float ComputeTraceMO0( const InputImageSizeType& rcomp );
        
    private:
        LMMSEVectorImageFilter(const Self &); // purposely not implemented
        void operator=(const Self &);         // purposely not implemented
        
        // The size of the nieghbourhood to compute the statistics:
        InputSizeType m_Radius;
        InputSizeType m_RadiusFeatures;
        // What should we do with negative values of the estimated square?
        bool m_UseAbsoluteValue;
        bool m_KeepValue;
        // The minimum number of voxels that we allow to compute local statistics and filter:
        unsigned int m_MinimumNumberOfUsedVoxelsFiltering;
        // The noise variance; this filter itself does not estimate this parameter, so
        // it should be supplied externally:
        float             m_Sigma;
        float             m_SigmaR; // Noise std in the R (x) projected channel
        float             m_SigmaG; // Noise std in the G (y) projected channel
        float             m_SigmaB; // Noise std in the B (z) projected channel
        float             m_H;
        unsigned int      m_NDWI;
        unsigned int      m_NBaselines;
        IndicatorType     m_DWI;
        IndicatorType     m_Baselines;
        GradientListType  m_GradientList;
        unsigned int      m_Neighbours;
        NeighboursIndType m_NeighboursInd;
        // This is used when a mask is applied to the DWI input:
        LabelImagePointer m_Mask;
        bool              m_SetZeroBck;
        // This is to tell the filter to NLM-filter (or not)
        // those locations with very large signal variabilities,
        // such as the background, the CSF or certain grey
        // matter locations
        bool              m_FilterOutliers;
        // This is to tell the filter not to apply the LMMSE
        // correction, turning out into a "simple" UNLM:
        bool              m_OnlyUNLM;
        
        // Where the features will be stored:
        FeaturesMapPointer   m_Featuresx;
        FeaturesMapPointer   m_Featuresy;
        FeaturesMapPointer   m_Featuresz;
    };
    
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLMMSEVectorImageFilter.txx"
#endif

#endif

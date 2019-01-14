/*=========================================================================
 
 Program:   Insight Segmentation & Registration Toolkit
 Module:    $RCSfile: itkLMMSEVectorImageFilter.txx,v $
 Language:  C++
 Date:      $Date: 2005/05/4 14:28:51 $
 Version:   $Revision: 1.1
 =========================================================================*/
#ifndef _itkLMMSEVectorImageFilter_txx
#define _itkLMMSEVectorImageFilter_txx
#include "itkLMMSEVectorImageFilter.h"

#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkMath.h"

/** &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */
// COMMENT THIS LINE TO AVOID THE DENUG CODE
//#define USE_DEBUG_CODE
/** &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */


/** &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */
// DEBUG:
#ifdef USE_DEBUG_CODE
#define DEBUG_FILE "/Users/atriveg/Downloads/WorkModes.nrrd"
#include "itkImageFileWriter.h"
#endif
/** &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */

namespace itk
{
    
    /** Constructor */
    template <class TInputImage, class TOutputImage>
    LMMSEVectorImageFilter<TInputImage, TOutputImage>::LMMSEVectorImageFilter()
    {
        m_Radius.Fill(1);
        m_RadiusFeatures.Fill(1);
        
        m_NDWI           = 0;
        m_NBaselines     = 0;
        m_DWI            = IndicatorType( 0 );
        m_Baselines      = IndicatorType( 0 );
        
        m_Sigma          = 20.0f;
        m_H              = 1.2;
        m_SigmaR         = 0;
        m_SigmaG         = 0;
        m_SigmaB         = 0;
        
        m_SetZeroBck     = false;
        m_OnlyUNLM       = false;
        m_FilterOutliers = false;
        
        m_GradientList   = GradientListType(0);
        m_Neighbours     = 1;    // By default, we use the gradient by gradient behaviour
        m_NeighboursInd  = NeighboursIndType(0, 0);
        
        m_Mask      = NULL;
        m_Featuresx = NULL;
        m_Featuresy = NULL;
        m_Featuresz = NULL;
/** &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */
// DEBUG:
#ifdef USE_DEBUG_CODE
        this->SetNumberOfThreads(1);
#endif
/** &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */
    }
    
    template <class TInputImage, class TOutputImage>
    void LMMSEVectorImageFilter<TInputImage, TOutputImage>
    ::BeforeThreadedGenerateData( void )
    {
        
        //======================================================================
        // BLOCK I: ORDER THE GRADIENTS AND COMPUTE THE NEIGHBORS
        if( m_Neighbours > m_NDWI )
        {
            m_Neighbours = m_NDWI;
        }
        // Find the closest neighbours to each gradient direction
        if( m_NDWI != m_DWI.GetSize() || m_NBaselines != m_Baselines.GetSize() ||
           (m_NDWI < 1 && m_NBaselines < 1) || m_GradientList.size() != m_NDWI || m_Neighbours < 1 || m_Neighbours > m_NDWI )
        {
            itkExceptionMacro( << "Bad initialisation of the filter!!! Check parameters, please" );
        }
        if( (m_NDWI+m_NBaselines) != this->GetInput()->GetVectorLength() )
            itkExceptionMacro( << "Bad initialisation of the filter!!! Check parameters, please" );
        m_NeighboursInd = NeighboursIndType( m_NDWI, m_Neighbours );
        
        // Vectors to compute the distance from each gradient direction to each other gradient direction; we need to sort to
        // find the closest
        // gradient directions to each of one.
        std::vector<OrderType> distances;
        OrderType              orderElement;
        for( unsigned int g = 0; g < m_NDWI; ++g )           // For each gradient direction
        {
            distances.clear();
            for( unsigned int k = 0; k < m_NDWI; ++k )       // Compare to each other gradient direction
            {
                orderElement[0] = (double)k;
                orderElement[1] = itk::NumericTraits<double>::Zero;
                for( unsigned int d = 0; d < TInputImage::ImageDimension; ++d ) // Sum of squared differences (euclidean norm)
                {
                    orderElement[1] += ( m_GradientList[g][d] * m_GradientList[k][d] );
                }
                if( orderElement[1] < -1.0f || orderElement[1] > 1.0f )
                {
                    orderElement[1] = 0.0f;
                }
                else
                {
                    orderElement[1] = ::acos( orderElement[1] );
                }
                if( 3.141592654f - orderElement[1] < orderElement[1] )
                {
                    orderElement[1] = 3.141592654f - orderElement[1];
                }
                distances.push_back( orderElement );
            }
            std::sort( distances.begin(), distances.end(), UNLM_gradientDistance_smaller );
            for( unsigned int k = 0; k < m_Neighbours; ++k )
            {
                m_NeighboursInd[g][k] = m_DWI[(unsigned int)(distances[k][0])];
            }
        }
        
        //======================================================================
        // BLOCK II: COMPUTE THE RGB PROJECTIONS FROM THE DWI CHANNELS:
        RGBProjectionPointer projx = RGBProjectionType::New();
        projx->SetInput( this->GetInput() );
        projx->SetProjectionCoordinate(0);
        projx->SetDWI( m_DWI );
        projx->SetGradientsTable( m_GradientList );
        RGBProjectionPointer projy = RGBProjectionType::New();
        projy->SetInput( this->GetInput() );
        projy->SetProjectionCoordinate(1);
        projy->SetDWI( m_DWI );
        projy->SetGradientsTable( m_GradientList );
        RGBProjectionPointer projz = RGBProjectionType::New();
        projz->SetInput( this->GetInput() );
        projz->SetProjectionCoordinate(2);
        projz->SetDWI( m_DWI );
        projz->SetGradientsTable( m_GradientList );
        //======================================================================
        // BLOCK III: COMPUTE THE SALIENT FEATURES RELATED TO THE PATCH DISTANCES:
        //  R (x) - Channel
        L0Pointer l0x = L0Type::New();
		L1Pointer l1x = L1Type::New();
		L2Pointer l2x = L2Type::New();
		l0x->SetRadius( this->GetRadiusFeatures()[0] );
		l0x->SetCoordinate( 0 );
		l1x->SetRadius( this->GetRadiusFeatures()[1] );
		l1x->SetCoordinate( 1 );
		l2x->SetRadius( this->GetRadiusFeatures()[2] );
		l2x->SetCoordinate( 2 );
        l0x->SetInput( projx->GetOutput() );
		l1x->SetInput( l0x->GetOutput() );
		l2x->SetInput( l1x->GetOutput() );
		//  G (y) - Channel
        L0Pointer l0y = L0Type::New();
		L1Pointer l1y = L1Type::New();
		L2Pointer l2y = L2Type::New();
		l0y->SetRadius( this->GetRadiusFeatures()[0] );
		l0y->SetCoordinate( 0 );
		l1y->SetRadius( this->GetRadiusFeatures()[1] );
		l1y->SetCoordinate( 1 );
		l2y->SetRadius( this->GetRadiusFeatures()[2] );
		l2y->SetCoordinate( 2 );
        l0y->SetInput( projy->GetOutput() );
		l1y->SetInput( l0y->GetOutput() );
		l2y->SetInput( l1y->GetOutput() );
        //  B (z) - Channel
        L0Pointer l0z = L0Type::New();
		L1Pointer l1z = L1Type::New();
		L2Pointer l2z = L2Type::New();
		l0z->SetRadius( this->GetRadiusFeatures()[0] );
		l0z->SetCoordinate( 0 );
		l1z->SetRadius( this->GetRadiusFeatures()[1] );
		l1z->SetCoordinate( 1 );
		l2z->SetRadius( this->GetRadiusFeatures()[2] );
		l2z->SetCoordinate( 2 );
        l0z->SetInput( projz->GetOutput() );
		l1z->SetInput( l0z->GetOutput() );
		l2z->SetInput( l1z->GetOutput() );
        
        //======================================================================
        // BLOCK III: Update the pipelines and keep the usable outputs:
        l2x->Update();
        l2y->Update();
        l2z->Update();
        m_Featuresx = l2x->GetOutput();
        m_Featuresy = l2y->GetOutput();
        m_Featuresz = l2z->GetOutput();
        // Compute the amount of residual noise in the projected images
        // as a function of the amount of noise in the original MRI
        // image and the weighting factors using in each channel.
        //
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
        m_SigmaR    = m_Sigma * projx->GetCorrectionFactor();
        m_SigmaG    = m_Sigma * projy->GetCorrectionFactor();
        m_SigmaB    = m_Sigma * projz->GetCorrectionFactor();
        
        return;
    }
    
    /** The requested input region is larger than the corresponding output, so we need to override this method: */
    template <class TInputImage, class TOutputImage>
    void LMMSEVectorImageFilter<TInputImage, TOutputImage>
    ::GenerateInputRequestedRegion()
    throw (InvalidRequestedRegionError)
    {
        // Call the superclass' implementation of this method
        Superclass::GenerateInputRequestedRegion();
        
        // Get pointers to the input and output
        InputImagePointer  inputPtr  = const_cast<TInputImage *>( this->GetInput() );
        OutputImagePointer outputPtr = this->GetOutput();
        
        if( !inputPtr || !outputPtr )
        {
            return;
        }
        
        // Get a copy of the input requested region (should equal the output
        // requested region)
        InputImageRegionType inputRequestedRegion = inputPtr->GetRequestedRegion();
        
        // Pad the input requested region by the operator radius
        inputRequestedRegion.PadByRadius( m_Radius );
        
        // Crop the input requested region at the input's largest possible region
        inputRequestedRegion.Crop(inputPtr->GetLargestPossibleRegion() );
        inputPtr->SetRequestedRegion( inputRequestedRegion );
        return;
    }
    
    
    /**
     This method completely differs from the former implementation. Instead of 
     having square (cubic) neighborhoods, the neighborhoods have an arbitrary
     shape adapted to the actual contents of the image: we use weigts computed
     in a similar fashion as in the non-local means approach to account only
     for those neighbors whose diffusion structure is similar enough to that
     of the voxel being studied.
     */
    template <class TInputImage, class TOutputImage>
    void LMMSEVectorImageFilter<TInputImage, TOutputImage>
    ::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread,
                           ThreadIdType itkNotUsed(threadId) )
    {
        //==================================================================================================================================
        // Iterators:
        ImageRegionConstIteratorWithIndex<FeaturesMapType> mitx;    // Iterator for the map of local featrues
        ImageRegionConstIteratorWithIndex<FeaturesMapType> mity;    // Iterator for the map of local featrues
        ImageRegionConstIteratorWithIndex<FeaturesMapType> mitz;    // Iterator for the map of local featrues
        ImageRegionIterator<OutputImageType>               it;      // Iterator for the output image
        ImageRegionConstIterator<InputImageType>           bit;     // Iterator for the output image
        ImageRegionConstIterator<InputImageType>           search;  // Search iterator
        ImageRegionConstIterator<FeaturesMapType>          msitx;   // Iterator for search in the map of local features
        ImageRegionConstIterator<FeaturesMapType>          msity;   // Iterator for search in the map of local features
        ImageRegionConstIterator<FeaturesMapType>          msitz;   // Iterator for search in the map of local features
        // Input and output
        InputImageConstPointer   input   =  this->GetInput();
        OutputImagePointer       output  =  this->GetOutput();
        
        //==================================================================================================================================        
        unsigned int numNeighbours = 1;
        InputImageSizeType baseSearchSize, searchSize;
        for( unsigned int d=0; d<TInputImage::ImageDimension; ++d ){
            // The number of voxels which are going to be accounted in the WA
            numNeighbours    *= ( 2*m_Radius[d] + 1 );
            baseSearchSize[d] = ( 2*m_Radius[d] + 1 );
        }
        InputImageRegionType searchRegion;
        
        //==================================================================================================================================
        float normNoisex = ( m_H * m_SigmaR * m_SigmaR ) * ComputeTraceMO1( this->GetRadiusFeatures() );
        normNoisex       = 1.0f/normNoisex;
        float normNoisey = ( m_H * m_SigmaG * m_SigmaG ) * ComputeTraceMO1( this->GetRadiusFeatures() );
        normNoisey       = 1.0f/normNoisey;
        float normNoisez = ( m_H * m_SigmaB * m_SigmaB ) * ComputeTraceMO1( this->GetRadiusFeatures() );
        normNoisez       = 1.0f/normNoisez;
        float lsnorm[TInputImage::ImageDimension];
        for( unsigned int k=0; k<TInputImage::ImageDimension; ++k ){
            lsnorm[k] = itk::NumericTraits<float>::Zero;
            //=====================================================================
            float* weight = new float[ m_RadiusFeatures[k] ];
            float  wsum   = itk::NumericTraits<float>::Zero;
            for( int j=0; j<((int)m_RadiusFeatures[k]); ++j ){
                weight[j]  = ::exp( -((float)(m_RadiusFeatures[k]-j)*(m_RadiusFeatures[k]-j))/2.0f );
                wsum      += 2.0f*weight[j];
            }
            wsum += weight[m_RadiusFeatures[k]-1];
            wsum  = 1.0f/wsum;
            //=====================================================================
            for( int j=-((int)m_RadiusFeatures[k]); j<0; ++j )
                lsnorm[k] += 2.0f * j*j * weight[j+m_RadiusFeatures[k]] * wsum;
            //=====================================================================
            delete[] weight;
            lsnorm[k]  = 1.0f/lsnorm[k];
        }
        // This constant is used to assess the theoretical self-similarity of
        // the central pixel to avoid over-weighting:
        double centerSelfSimilarity = std::exp( -itk::NumericTraits<double>::One / m_H );
        //==================================================================================================================================
        // CREATE THE ITERATORS:
        mitx = ImageRegionConstIteratorWithIndex<FeaturesMapType>( m_Featuresx, outputRegionForThread );
        mity = ImageRegionConstIteratorWithIndex<FeaturesMapType>( m_Featuresy, outputRegionForThread );
        mitz = ImageRegionConstIteratorWithIndex<FeaturesMapType>( m_Featuresz, outputRegionForThread );
        bit  = ImageRegionConstIterator<InputImageType>(           input,       outputRegionForThread );
        it   = ImageRegionIterator<OutputImageType>(               output,      outputRegionForThread );
        InputImageIndexType originR;
        InputImageSizeType  radiusR;
        radiusR = m_Radius;
        
        //==================================================================================================================================
        // ALLOCATE MEMORY FOR THE VECTORS OF MOMENTS TO BE COMPUTED:
        double* diff                  = new double[m_NDWI + m_NBaselines];
        double* dSecondAveragedMoment = new double[m_NDWI + m_NBaselines];
        double* dSquaredMagnitude     = new double[m_NDWI + m_NBaselines];
        double* dFiltered             = new double[m_NDWI + m_NBaselines];
        double* dFourthAveragedMoment = new double[m_NDWI + m_NBaselines];
        
        double* bSqMag = new double[m_NBaselines];
        double* bSqAvg = new double[m_NBaselines];
        double* bRes   = new double[m_NBaselines];
        double* dSqMag = new double[m_Neighbours];
        double* dSqAvg = new double[m_Neighbours];
        double* dRes   = new double[m_Neighbours];
        
        //==================================================================================================================================
        // PREPARE THE MASK IN CASE IT IS USED:
        MaskIteratorType maskIterator;
        if( m_Mask ){
            maskIterator = MaskIteratorType( m_Mask, outputRegionForThread );
            maskIterator.GoToBegin();
        }
        
/** &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */
// DEBUG:
#ifdef USE_DEBUG_CODE
        typedef itk::Image<unsigned char,3> DebugImageType;
        typedef DebugImageType::Pointer    DebugImagePointer;
        DebugImagePointer debugImage = DebugImageType::New();
        debugImage->SetOrigin( input->GetOrigin() );
        debugImage->SetSpacing( input->GetSpacing() );
        debugImage->SetDirection( input->GetDirection() );
        debugImage->SetRegions( input->GetLargestPossibleRegion() );
        debugImage->Allocate();
        debugImage->FillBuffer(0);
        ImageRegionIterator<DebugImageType> debugIt = ImageRegionIterator<DebugImageType>( debugImage, outputRegionForThread );
        debugIt.GoToBegin();
#endif
/** &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */
        //==================================================================================================================================
        // DO THE ACTUAL COMPUTATIONS:
        for( it.GoToBegin(),bit.GoToBegin(),mitx.GoToBegin(),mity.GoToBegin(),mitz.GoToBegin();
            !it.IsAtEnd();
            ++it,++bit,++mitx,++mity,++mitz ){
            
            //-------------------------------------------------------------------------------------------------------------
            // In case we use a mask, there is a chance we have to pass the input directly
            // to the output:
            if( m_Mask ){
                if( !maskIterator.Get() ){ // The value of the mask is 0
                    if( m_SetZeroBck ){
                        OutputPixelType outpx( this->GetInput()->GetVectorLength() );
                        outpx.Fill( itk::NumericTraits<ScalarType>::Zero );
                        it.Set( outpx );
                    }
                    else
                        it.Set( bit.Get() );
                    ++maskIterator;        // increment the iterator
/** &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */
// DEBUG:
#ifdef USE_DEBUG_CODE
                    debugIt.Set(0);
                    ++debugIt;
#endif
/** &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */
                    continue;              // and go for the next pixel
                }
                else
                    ++maskIterator;        // only increment the iterator;
            }
            
            //-------------------------------------------------------------------------------------------------------------
            // CREATE THE REGION TO SEARCH AND THE ITERATORS:
            searchSize = baseSearchSize;
            originR    = mitx.GetIndex() - radiusR;
            bool         needToComputeCenter = false;
            unsigned int midPosition         = numNeighbours/2;
            for( unsigned int d=0; d<TInputImage::ImageDimension; ++d ){
                if( originR[d]<0 ){
                    searchSize[d] += originR[d];
                    originR[d]     = 0;
                    needToComputeCenter = true;
                }
                if( originR[d]+searchSize[d] > input->GetLargestPossibleRegion().GetSize()[d] ){
                    searchSize[d]       = input->GetLargestPossibleRegion().GetSize()[d] - originR[d];
                    needToComputeCenter = true;
                }
            }
            
            //-------------------------------------------------------------------------------------------------------------
            // Compute the index corresponding to the original center:
            if( needToComputeCenter ){
                unsigned int aux = 1;
                for( unsigned int d=0; d<TInputImage::ImageDimension; ++d )
                    aux *= searchSize[d];
                midPosition = 0;
                if( aux>0 ){
                    for( int d=(int)TInputImage::ImageDimension-1; d>=0; --d ){
                        aux /= searchSize[d];
                        midPosition += ( mitx.GetIndex()[d] - originR[d] )*aux;
                    }
                }
            }
            
            //-------------------------------------------------------------------------------------------------------------
            // Initialize the search iterators:
            searchRegion.SetIndex( originR );
            searchRegion.SetSize( searchSize );
            search = ImageRegionConstIterator<InputImageType>(   input,       searchRegion  );
            msitx  = ImageRegionConstIterator<FeaturesMapType>(  m_Featuresx, searchRegion  );
            msity  = ImageRegionConstIterator<FeaturesMapType>(  m_Featuresy, searchRegion  );
            msitz  = ImageRegionConstIterator<FeaturesMapType>(  m_Featuresz, searchRegion  );
            
            //-------------------------------------------------------------------------------------------------------------
            // Initalize the vectors to compute the moments:
            for( unsigned int ch=0; ch<this->GetInput()->GetVectorLength(); ++ch ){
                dSecondAveragedMoment[ch] = itk::NumericTraits<double>::Zero;
                dFourthAveragedMoment[ch] = itk::NumericTraits<double>::Zero;
            }
            
            //-------------------------------------------------------------------------------------------------------------
            // FILTER THE PIXEL
            LSGradientsL2 centerx = mitx.Get();
            LSGradientsL2 centery = mity.Get();
            LSGradientsL2 centerz = mitz.Get();
            float norm     = itk::NumericTraits<float>::Zero;    // To normalize the weights to sum to 1
            float weight, weightx, weighty, weightz;
            unsigned int pos; // Auxiliar variable
            for( pos=0,search.GoToBegin(),msitx.GoToBegin(),msity.GoToBegin(),msitz.GoToBegin(); !search.IsAtEnd(); ++search,++msitx,++msity,++msitz,++pos ){
                // Compute the weight associated to the current voxel:
                if( pos!=midPosition ){
                    LSGradientsL2 valuex = msitx.Get();
                    LSGradientsL2 valuey = msity.Get();
                    LSGradientsL2 valuez = msitz.Get();
                    weightx  = (centerx.LLL-valuex.LLL)*(valuex.LLL-centerx.LLL);
                    weightx += (centerx.HLL-valuex.HLL)*(valuex.HLL-centerx.HLL)*lsnorm[0];
                    weightx += (centerx.LHL-valuex.LHL)*(valuex.LHL-centerx.LHL)*lsnorm[1];
                    weightx += (centerx.LLH-valuex.LLH)*(valuex.LLH-centerx.LLH)*lsnorm[2];
                    weightx *= normNoisex;
                    weighty  = (centery.LLL-valuey.LLL)*(valuey.LLL-centery.LLL);
                    weighty += (centery.HLL-valuey.HLL)*(valuey.HLL-centery.HLL)*lsnorm[0];
                    weighty += (centery.LHL-valuey.LHL)*(valuey.LHL-centery.LHL)*lsnorm[1];
                    weighty += (centery.LLH-valuey.LLH)*(valuey.LLH-centery.LLH)*lsnorm[2];
                    weighty *= normNoisey;
                    weightz  = (centerz.LLL-valuez.LLL)*(valuez.LLL-centerz.LLL);
                    weightz += (centerz.HLL-valuez.HLL)*(valuez.HLL-centerz.HLL)*lsnorm[0];
                    weightz += (centerz.LHL-valuez.LHL)*(valuez.LHL-centerz.LHL)*lsnorm[1];
                    weightz += (centerz.LLH-valuez.LLH)*(valuez.LLH-centerz.LLH)*lsnorm[2];
                    weightz *= normNoisez;
                    weight   = std::exp( (weightx+weighty+weightz)/3 );
                    norm    += weight;
                }
                else{
                    weight   = centerSelfSimilarity;
                    norm    += weight;
                    // In the center of the neighborhood we have to keep the
                    // non-filtered value too:
                    InputPixelType ipx = search.Get();
                    for( unsigned int ch=0; ch<this->GetInput()->GetVectorLength(); ++ch ){
                        double pix  = ipx[ch];
                        dSquaredMagnitude[ch] = pix*pix;
                    }
                }
                // Compute the actual moments:
                InputPixelType cipx = search.Get();
                for( unsigned int ch=0; ch<this->GetInput()->GetVectorLength(); ++ch ){
                    double pix  = cipx[ch];
                    pix        *= pix;
                    dSecondAveragedMoment[ch] += ( pix * weight );
                    pix        *= pix;
                    dFourthAveragedMoment[ch] += ( pix * weight );
                }
            }
            
            //-------------------------------------------------------------------------------------------------------------
            // Now we have searched all the neighborhood, we can normalize the 
            // sums to compute the actual moments; these are the moments of the
            // measurements M, so we need to correct them to compute the moments
            // of the original magnitude A^2:
            norm = itk::NumericTraits<float>::One / norm;
            for( unsigned int ch=0; ch<this->GetInput()->GetVectorLength(); ++ch ){
                dSecondAveragedMoment[ch] *= norm;
                dFourthAveragedMoment[ch] *= norm;
                diff[ch]                   = dSquaredMagnitude[ch] - dSecondAveragedMoment[ch];
                dSecondAveragedMoment[ch] -= 2*m_Sigma*m_Sigma;
                if( dSecondAveragedMoment[ch] < 100000 * std::numeric_limits<double>::epsilon() )
                    dSecondAveragedMoment[ch] = 100000 * std::numeric_limits<double>::epsilon();
                dFourthAveragedMoment[ch] -= 8*m_Sigma*m_Sigma*( dSecondAveragedMoment[ch] + m_Sigma*m_Sigma );
                if( dFourthAveragedMoment[ch] < 100000 * std::numeric_limits<double>::epsilon() )
                    dFourthAveragedMoment[ch] = 100000 * std::numeric_limits<double>::epsilon();
            }
            
            // -----------------------------------------------------------------------------------------------------------------------
            // Now, we have estimates of the moments of A. We have computed as well the difference M - E{M}, that has to be
            // filtered with the inverse of the covariance matrix C_M2M2.
            const unsigned int MAX_ALLOWED_VAR = 1000;
            const float        CFACT1          = 5.0f;
            OutputPixelType outPixel = bit.Get(); // Auxiliar output pixel
            //    -Normalization factor:
            unsigned int count = 0;
            double normal      = itk::NumericTraits<double>::Zero;
            for( unsigned int ch=0; ch<this->GetInput()->GetVectorLength(); ++ch ){
                double dsqMVar = dFourthAveragedMoment[ch] - dSecondAveragedMoment[ch]*dSecondAveragedMoment[ch];
                if( dsqMVar>0 ){
                    if( dSecondAveragedMoment[ch] > 100 * std::numeric_limits<double>::epsilon() ){
                        normal += ( dsqMVar / (dSecondAveragedMoment[ch]*dSecondAveragedMoment[ch]) );
                        count++;
                    }
                }
            }
            if( count>0 )
                normal /= count;
            
            // -----------------------------------------------------------------------------------------------------------------------
            // If the "OnlyUNLM" mode has been set, we just fix count=0 here,
            // so that we go straight to the "else" statement, i.e. we simply
            // keep the second order moment <A^2>0<M^2> - 2·sigma^2. This is 
            // equivalent to an unbiased non-local means whose widths are
            // computed from the RGB projections:
            if( m_OnlyUNLM )
                count = 0;
            
            //    - Background checking:
            if( count >= m_NBaselines ){
                //    - Variability checking:
                if( normal <= 100 * std::numeric_limits<double>::epsilon() ){
                    // The Variability is extremely low, so it is  likely that an
                    // homogeneous region is being filtered. In this case,
                    // ||C_A2A2|| << ||C_M2M2||, so we simply use the unbiased
                    // estimate of the second order  moment:
                    for( unsigned int ch=0; ch<this->GetInput()->GetVectorLength(); ++ch )
                        dFiltered[ch] = dSecondAveragedMoment[ch];
/** &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */
// DEBUG:
#ifdef USE_DEBUG_CODE
                    debugIt.Set( 50 );
#endif
/** &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */
                }
                else if( normal > MAX_ALLOWED_VAR )
                {
                    // The variability is too high, so C_M2M2 is close to singular
                    // and numerical problems could arise.
                    if( m_FilterOutliers ){
                        for( unsigned int ch=0; ch<this->GetInput()->GetVectorLength(); ++ch )
                            dFiltered[ch] = dSecondAveragedMoment[ch];
                    }
                    else{
                        for( unsigned int ch=0; ch<this->GetInput()->GetVectorLength(); ++ch )
                            dFiltered[ch] = dSquaredMagnitude[ch];
                    }
/** &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */
// DEBUG:
#ifdef USE_DEBUG_CODE
                    debugIt.Set( 100 );
#endif
/** &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */
                }
                else
                {
                    // This is the normal case, and should be the one present in the majority of the voxels of the image
                    // -----------------------------------------------------------------------------------------------------------------------
                    // First, filter the baseline images, all together:
                    double minSqAvg = itk::NumericTraits<double>::max();
                    for( unsigned int ch=0; ch<m_NBaselines; ++ch ){
                        bSqMag[ch] = diff[ m_Baselines[ch] ];
                        bSqAvg[ch] = dSecondAveragedMoment[ m_Baselines[ch] ];
                        if( bSqAvg[ch] < minSqAvg )
                            minSqAvg = bSqAvg[ch];
                    }
                    //    - Pre-whitening of the input:
                    if( minSqAvg > CFACT1*m_Sigma*m_Sigma ){
                        // In this case the power series expansion is convergent:
                        this->CMMInversion( bSqMag, bSqAvg, normal, bRes, 10, m_NBaselines );
                    }
                    else{
                        // The serie expansion is not convergent, and the linear
                        // correction is not stable; the aproximation is
                        // not accurate, but this corresponds mainly to background
                        // pixels, so it is not so important
                        this->ComputeInverseMatrix( bSqMag, bSqAvg, normal, bRes, m_NBaselines );
                    }
                    //    - Product with C_A2M2
                    //          Scalar product with the vector of second order moments:
                    double dp = itk::NumericTraits<double>::Zero;
                    for( unsigned int ch=0; ch<m_NBaselines; ++ch )
                        dp += bRes[ch] * bSqAvg[ch];
                    //    - Correction of the output value:
                    for( unsigned int ch=0; ch<m_NBaselines; ++ch )
                        dFiltered[m_Baselines[ch]] = (1.0f+normal*dp) * bSqAvg[ch];
                    // -----------------------------------------------------------------------------------------------------------------------
                    // Now, filter the gradient images
                    unsigned int top = m_NDWI;
                    if( m_Neighbours == m_NDWI )
                        top = 1;
/** &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */
// DEBUG:
#ifdef USE_DEBUG_CODE
                    unsigned char myMode = 150;
#endif
/** &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */
                    for( unsigned int g=0; g<top; ++g ){
                        minSqAvg = itk::NumericTraits<double>::max(); // Initialise maximum
                        // Generate the vector with the appropriate measures, i.e., the ones from the closest gradient directions
                        for( unsigned int ch=0; ch<m_Neighbours; ++ch ){
                            dSqMag[ch] = diff[ m_NeighboursInd[g][ch] ];
                            dSqAvg[ch] = dSecondAveragedMoment[ m_NeighboursInd[g][ch] ];
                            if( dSqAvg[ch] < minSqAvg )
                                minSqAvg = dSqAvg[ch];
                        }
                        //    - Pre-whitening of the input:
                        if( minSqAvg > CFACT1*m_Sigma*m_Sigma ){
                            // In this case, the series expansion is convergent, so we may
                            // perform the linear correction
                            this->CMMInversion( dSqMag, dSqAvg, normal, dRes, 10, m_Neighbours );
                        }
                        else{
                            // The series expansion is not convergent, and the linear correction is not stable; the aproximation is
                            // not accurate, but this corresponds mainly to background pixels, so it is not so important
                            this->ComputeInverseMatrix( dSqMag, dSqAvg, normal, dRes, m_Neighbours );
/** &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */
// DEBUG:
#ifdef USE_DEBUG_CODE
                            myMode = 200;
#endif
/** &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */
                        }
                        //    - Product with C_A2M2
                        //          Scalar product with the vector of second order moments:
                        dp = itk::NumericTraits<double>::Zero;
                        for( unsigned int ch=0; ch<m_Neighbours; ++ch )
                            dp += dRes[ch] * dSqAvg[ch];
                        if( m_Neighbours==m_NDWI ){
                            //    - Correction of the output value:
                            for( unsigned int ch=0; ch<m_Neighbours; ++ch )
                                dFiltered[m_NeighboursInd[g][ch]] = (1.0f + normal*dp) * dSqAvg[ch];
                        }
                        else
                            dFiltered[m_DWI[g]] = (1.0f + normal*dp) * dSqAvg[0];
                    }
/** &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */
// DEBUG:
#ifdef USE_DEBUG_CODE
                    debugIt.Set( myMode );
#endif
/** &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */
                }
                // Compute the square root of the output, and check if the result is physically consisitent:
                for( unsigned int ch=0; ch<this->GetInput()->GetVectorLength(); ++ch ){
                    if( dFiltered[ch] > 0 )
                        dFiltered[ch] = std::sqrt( dFiltered[ch] );
                    else
                        dFiltered[ch] = 0;
                }
            }
            else{ // In this case, the second order moment is too small; this is likely to occur in the background
                for( unsigned int ch=0; ch<this->GetInput()->GetVectorLength(); ++ch ){
                    dFiltered[ch] = dSecondAveragedMoment[ch];
                    if( dFiltered[ch]>0 )
                        dFiltered[ch] = std::sqrt(dFiltered[ch]);
                    else
                        dFiltered[ch] = itk::NumericTraits<double>::Zero;
                }
/** &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */
// DEBUG:
#ifdef USE_DEBUG_CODE
                debugIt.Set( 250 );
#endif
/** &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */
            }
            
            //-------------------------------------------------------------------------------------------------------------
            // FINALLY, SET THE OUTPUT PIXEL
            for( unsigned int ch=0; ch<this->GetInput()->GetVectorLength(); ++ch )
                outPixel[ch] = static_cast<ScalarType>( dFiltered[ch] );
            it.Set(   static_cast<OutputPixelType>( outPixel )   );
/** &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */
// DEBUG:
#ifdef USE_DEBUG_CODE
            ++debugIt;
#endif
/** &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */
        }
        
/** &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */
// DEBUG:
#ifdef USE_DEBUG_CODE
        typedef itk::ImageFileWriter<DebugImageType> DebugImageWriterType;
        typedef DebugImageWriterType::Pointer        DebugImageWriterPointer;
        DebugImageWriterPointer debugWriter = DebugImageWriterType::New();
        debugWriter->SetInput( debugImage );
        debugWriter->SetFileName( DEBUG_FILE );
        debugWriter->Update();
#endif
/** &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */
        
        //==================================================================================================================================
        // CLEAR ALL MANUALLY ALLOCATED MEMORY:
        delete[] diff;
        delete[] dSecondAveragedMoment;
        delete[] dSquaredMagnitude;
        delete[] dFiltered;
        delete[] dFourthAveragedMoment;
        delete[] bSqMag;
        delete[] bSqAvg;
        delete[] bRes;
        delete[] dSqMag;
        delete[] dSqAvg;
        delete[] dRes;
    }
    
    /** Smart approximate inversion of C_{M^2M^2} (high SNR case)*/
    template <class TInputImage, class TOutput>
    void LMMSEVectorImageFilter<TInputImage, TOutput>
    ::CMMInversion( const double* measures, const double* squaredAverages, double normal, double* whitened,
                   unsigned int order,
                   unsigned int K ) const
    {
        // Where:
        //     measures: the squared measurements, which is, the original data (one per channel)
        //     squaredAverages: The vector containing the second order moment for each DWI channel
        //     normal: the variance of the second order moment normalised by the square of the second order moment
        //     whitened: the processed signal, which is, C_MM^(-1)*(M^2-E{M^2})
        //     order: the number of iterations, i.e., the order of Taylor series expansion
        // Auxiliar value to precompute constants:
        if( K == 1 )
        {
            double var  = m_Sigma * m_Sigma;
            double aux  = squaredAverages[0];
            aux         = normal * aux * aux + 4.0f * var * aux + 4.0f * var * var;
            whitened[0] = measures[0] / aux;
            return;
        }
        normal     = itk::NumericTraits<double>::One / normal; // For convenience
        double aux = 4.0f * m_Sigma * m_Sigma * normal;
        // The terms in the inverse matrix:
        double  Ad = aux;
        double* Ai = new double[K];
        for( unsigned int k = 0; k < K; ++k )
        {
            Ad   += squaredAverages[k];
            Ai[k] = itk::NumericTraits<double>::One / ( aux * squaredAverages[k] );
        }
        Ad     = -itk::NumericTraits<double>::One / ( aux * Ad );
        // Now, recursively process the output; initiallise w_0 = x
        for( unsigned int k = 0; k < K; ++k )
        {
            whitened[k] = measures[k];
        }
        double cum; // Auxiliar value
        aux *= (m_Sigma * m_Sigma);
        // Iterate: w_{n+1} = x - D^{-1}w_n
        for( unsigned int o = 0; o < order; ++o )        // If order=0, this loop does nothing!
        { // Compute A_d*w
            cum = itk::NumericTraits<double>::Zero;  // Initiallise acumulator
            for( unsigned int k = 0; k < K; ++k )
            {
                cum += whitened[k];
            }
            cum *= Ad;
            // Compute A_i*w
            for( unsigned int k = 0; k < K; ++k )
            {
                whitened[k] = measures[k] - aux * (   Ai[k] * whitened[k] + cum   );
            }
        }
        // Now we have the truncated series of ( Id + D^(-1) )^(-1). It remains to
        // multiplicate by D^(-1):
        // Compute A_d*w
        cum = itk::NumericTraits<double>::Zero; // Initiallise acumulator
        for( unsigned int k = 0; k < K; ++k )
        {
            cum += whitened[k];
        }
        cum *= Ad;
        // Compute A_i*w + A_d*w and normalise
        for( unsigned int k = 0; k < K; ++k )
        {
            whitened[k] = (   Ai[k] * whitened[k] + cum   ) * normal;
        }
        // Delete allocated memory:
        delete[] Ai;
        return;
    }
    
    /** Matrix inversion; the general case */
    template <class TInputImage, class TOutput>
    bool LMMSEVectorImageFilter<TInputImage, TOutput>
    ::ComputeInverseMatrix( const double* measures, const double* squaredAverages, double normal, double* whitened,
                           unsigned int K ) const
    {
        if( K == 1 )
        {
            double var  = m_Sigma * m_Sigma;
            double aux  = squaredAverages[0];
            aux         = normal * aux * aux + 4.0f * var * aux + 4.0f * var * var;
            whitened[0] = measures[0] / aux;
            return true;
        }
        // Compute the matrix to invert
        double* * matrix = new double *[K];
        for( unsigned int j = 0; j < K; ++j )
        {
            matrix[j]    = new double[K];
        }
        for( unsigned int j = 0; j < K; ++j )
        {
            matrix[j][j] = normal * squaredAverages[j] * squaredAverages[j] + 4 * m_Sigma * m_Sigma
            * (squaredAverages[j] + m_Sigma * m_Sigma);
            for( unsigned int k = j + 1; k < K; ++k )
            {
                matrix[j][k] = normal * squaredAverages[j] * squaredAverages[k];
                matrix[k][j] = matrix[j][k];
            }
        }
        // Compute the independent term:
        double* iterm = new double[K];
        for( unsigned int j = 0; j < K; ++j )
        {
            iterm[j] = measures[j];
        }
        // For each column col = 1 to m_Channels-1, we need to make zeros in rows from
        // col+1 to m_Channels (note that in C++ array indices are 0-based):
        for( unsigned int col = 0; col < K - 1; ++col )  // For each column
        { // We need a non-null element in the position (col,col), in order to
            // accomplish gaussian elimination:
            if( fabs(matrix[col][col]) <= 1e-10 )
            {
                // Bad luck! The element is zero. We need to add a complete row to
                // the row in position c, so that the new element in position (c,c)
                // is not null. Find the first row for which the element (row,col)
                // is non-zero:
                unsigned int row = col + 1;
                while( fabs(matrix[row][col]) <= 1e-10 && row < K )
                {
                    ++row;
                }
                
                // If we are not able to find a row satisfying this condition, then
                // the matrix is singular, and this should not be the case; for
                // this reason, we do not perform bound checking, for efficiency. We
                // assume that row is a valid position, and then correct the input
                // and output:
                if( row == K )  // Singular matrix!!!
                {
                    for( unsigned int j = 0; j < K; ++j )
                    {
                        delete[] matrix[j];
                    }
                    delete[] matrix;
                    delete[] iterm;
                    return false;
                }
                for( unsigned int cc = col; cc < K; ++cc )
                {
                    matrix[col][cc]  += matrix[row][cc];
                }
                iterm[col] += iterm[row];
            }
            // At this point, we have a valid (col,col), element. We scale the whole
            // corresponding col-th row so that the pivoting element is simply 1:
            double scale = itk::NumericTraits<double>::One / matrix[col][col];
            for( unsigned int cc = col; cc < K; ++cc )
            {
                matrix[col][cc]  *= scale;
            }
            iterm[col] *= scale;
            // Now, we may perform gaussian elimination for each row:
            for( unsigned int row = col + 1; row < K; ++row )  // For each row
            {
                scale = matrix[row][col]; // This is the scale, since input[col][col] = 1.
                // Once again, for each column, we add the corresponding scaled
                // version of the pivoting element; however, in the input matrix,
                // values at the left of this column are assumed to be already zero:
                for( unsigned int cc = col; cc < K; ++cc ) // Only the columns from col
                {
                    matrix[row][cc] -= scale * matrix[col][cc];
                }
                iterm[row] -= scale * iterm[col];
                // We have completed this row
            }
            // We have completed this column
        }
        // Now we have an upper-triangular matrix, where all diagonal elements are
        // just 1, except for the last one; Now, we may compute the output in a recursive
        // fashion:
        if( fabs(matrix[K - 1][K - 1]) <= 1e-10 )
        {
            for( unsigned int j = 0; j < K; ++j )
            {
                delete[] matrix[j];
            }
            delete[] matrix;
            delete[] iterm;
            return false;
        }
        whitened[K - 1] = iterm[K - 1] / matrix[K - 1][K - 1]; // The last one
        for( int k = K - 2; k >= 0; --k )                      // For each component
        {
            whitened[k] = iterm[k]; // Initiallise
            for( unsigned int j = k + 1; j < K; ++j )
            {
                whitened[k] -= whitened[j] * matrix[k][j];
            }
        }
        // Delete allocated memory:
        for( unsigned int j = 0; j < K; ++j )
        {
            delete[] matrix[j];
        }
        delete[] matrix;
        delete[] iterm;
        // Matrix has been inverted!!
        return true;
    }
    
    template< class TInputImage, class TOutputImage >
    float LMMSEVectorImageFilter<TInputImage, TOutputImage >
    ::ComputeTraceMO0( const InputImageSizeType& rcomp )
    {
        unsigned int size = 1;
        for( unsigned int k=0; k<TInputImage::ImageDimension; ++k )
            size *= (2*rcomp[k]+1);
        typedef itk::ConstNeighborhoodIterator<InputImageType> IteratorType;
        IteratorType bit = IteratorType( rcomp, this->GetInput(), this->GetInput()->GetBufferedRegion()  );
        typename IteratorType::OffsetType idx;
        bit.GoToBegin();
        float norm  = itk::NumericTraits<float>::Zero;
        float trace = itk::NumericTraits<float>::Zero;
        for( unsigned int k=0; k<size/2; ++k ){
            idx    = bit.GetOffset(k);
            float aux = itk::NumericTraits<float>::Zero;
            for( unsigned int j=0; j<TInputImage::ImageDimension; ++j )
                aux += ((float)idx[j])*((float)idx[j]);
            norm  += ::exp(-aux/2);
            trace += ::exp(-aux);
        }
        norm  = 2.0f*norm  + ::exp(-0.5f);
        trace = 2.0f*trace + ::exp(-1.0f);
        return(trace/norm/norm);
    }
    
    
    template< class TInputImage, class TOutputImage >
    float LMMSEVectorImageFilter<TInputImage, TOutputImage >
    ::ComputeTraceMO1( const InputImageSizeType& rcomp )
    {
        unsigned int size = 1;
        for( unsigned int k=0; k<TInputImage::ImageDimension; ++k )
            size *= (2*rcomp[k]+1);
        typedef itk::ConstNeighborhoodIterator<InputImageType> IteratorType;
        IteratorType bit = IteratorType( rcomp, this->GetInput(), this->GetInput()->GetBufferedRegion()  );
        typename IteratorType::OffsetType idx;
        bit.GoToBegin();
        float norm  = itk::NumericTraits<float>::Zero;
        float trace = itk::NumericTraits<float>::Zero;
        for( unsigned int k=0; k<size/2; ++k ){
            idx    = bit.GetOffset(k);
            float aux = itk::NumericTraits<float>::Zero;
            for( unsigned int j=0; j<TInputImage::ImageDimension; ++j )
                aux += ((float)idx[j])*((float)idx[j]);
            norm  += ::exp(-aux/2);
            trace += ::exp(-aux);
        }
        norm  = 2.0f*norm  + ::exp(-0.5f);
        trace = 2.0f*trace + ::exp(-1.0f);
        trace = trace/norm/norm;
        if( TInputImage::ImageDimension==2 )
            trace = 30.0f*trace*trace; 
        else if( TInputImage::ImageDimension==3 )
            trace = 126.6f*trace*trace; 
        else
            trace = 0.1f;
        return(trace);
    }
    
    /** Standard "PrintSelf" method */
    template <class TInputImage, class TOutput>
    void LMMSEVectorImageFilter<TInputImage, TOutput>
    ::PrintSelf( std::ostream& os, Indent indent ) const
    {
        Superclass::PrintSelf( os, indent );
        os << indent << "Radius: "                             << m_Radius                             << std::endl;
        os << indent << "UseAbsoluteValue: "                   << m_UseAbsoluteValue                   << std::endl;
        os << indent << "KeepValue: "                          << m_KeepValue                          << std::endl;
        os << indent << "MinimumNumberOfUsedVoxelsFiltering: " << m_MinimumNumberOfUsedVoxelsFiltering << std::endl;
    }
    
} // end namespace itk

#endif

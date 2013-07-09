/*=========================================================================
 
 Program:   Diffusion Applications
 Language:  C++
 Module:    $HeadURL: http://svn.slicer.org/Slicer4/trunk/Modules/CLI/DiffusionApplications/dwiNoiseFilter/dwiNoiseFilter.cxx $
 Date:      $Date: 2008-11-25 18:46:58 +0100 (Tue, 25 Nov 2008) $
 Version:   $Revision: 7972 $
 
 Copyright (c) Brigham and Women's Hospital (BWH) All Rights Reserved.
 
 See License.txt or http://www.slicer.org/copyright/copyright.txt for details.
 
 ==========================================================================*/
#include <itkMetaDataObject.h>

#ifdef _WIN32
// to pick up M_SQRT2 and other nice things...
#define _USE_MATH_DEFINES
#endif

// DWIJointRicianLMMSEFilter includes
#include "itkLMMSEVectorImageFilter.h"

// Noise estimation includes
#include "NoiseEstimation/itkComputeRestrictedHistogram.h"
#include "NoiseEstimation/itkOtsuStatistics.h"
#include "NoiseEstimation/itkOtsuThreshold.h"
#include "NoiseEstimation/itkThresholdToMaskImageFilter.h"
#include "NoiseEstimation/itkVectorLocalStdImageFilter.h"

#include "SlicerJointRicianAnisotropicLMMSEFilterCLP.h"

// CLI includes
#ifdef THISISASLICERBUILD

#include "itkPluginUtilities.h"

#else

#ifdef SLICERV4
#include "itkPluginUtilities.h"
#else
#include "SlicerExecutionModel/itkPluginUtilities.h"
#endif

#endif

// ITK includes
#include <itkImageFileWriter.h>
#include <itkNrrdImageIO.h>
#include <itkCastImageFilter.h>

#define DIMENSION 3
#define dwiPI 3.141592653589793

template <class PixelType>
int DoIt( int argc, char * argv[], PixelType )
{
    PARSE_ARGS;
    
    // do the typedefs
    
    typedef itk::VectorImage<PixelType, DIMENSION>       DiffusionImageType;
    
    typedef itk::Image<PixelType, DIMENSION>             ScalarImageType;
    
    typedef float                                        PixelTypeFloat;
    typedef itk::VectorImage<PixelTypeFloat, DIMENSION>  FloatDiffusionImageType;
    typedef itk::Image<PixelTypeFloat, DIMENSION>        ScalarFloatImageType;
    
    typedef itk::CovariantVector<double, DIMENSION>      CovariantVectorType;
    
    std::vector<CovariantVectorType> diffusionDirections;
    
    typedef itk::ImageFileReader<DiffusionImageType> FileReaderType;
    typename FileReaderType::Pointer reader = FileReaderType::New();
    reader->SetFileName( inputVolume.c_str() );
    reader->Update();
    
    itk::MetaDataDictionary            imgMetaDictionary = reader->GetMetaDataDictionary();
    std::vector<std::string>           imgMetaKeys = imgMetaDictionary.GetKeys();
    std::string                        metaString;
    
    std::cout << "Number of keys = " << imgMetaKeys.size() << std::endl;
    
    typedef itk::MetaDataDictionary DictionaryType;
    const DictionaryType & dictionary = reader->GetMetaDataDictionary();
    
    typedef itk::MetaDataObject<std::string> MetaDataStringType;
    
    DictionaryType::ConstIterator itr = dictionary.Begin();
    DictionaryType::ConstIterator end = dictionary.End();
    
    double       dBValue = 1000;
    int          iFoundBValue = 0;
    unsigned int channels = 0;
    while( itr != end )
    {
        itk::MetaDataObjectBase::Pointer entry = itr->second;
        MetaDataStringType::Pointer      entryvalue =
        dynamic_cast<MetaDataStringType *>( entry.GetPointer() );
        
        if( entryvalue )
        {
            ::size_t pos = itr->first.find("DWMRI_gradient");
            
            if( pos != std::string::npos )
            {
                std::string tagkey = itr->first;
                std::string tagvalue = entryvalue->GetMetaDataObjectValue();
                
                double dx[DIMENSION];
                std::sscanf(tagvalue.c_str(), "%lf %lf %lf\n", &dx[0], &dx[1], &dx[2]);
                diffusionDirections.push_back( (CovariantVectorType)(dx) );
                ++channels;
            }
            else
            {
                
                // try to find the b-value
                
                ::size_t posB = itr->first.find("DWMRI_b-value");
                
                if( posB != std::string::npos )
                {
                    std::string tagvalue = entryvalue->GetMetaDataObjectValue();
                    std::sscanf(tagvalue.c_str(), "%lf\n", &dBValue );
                    iFoundBValue = 1;
                }
                else
                {
                    // std::cout << itr->first << " " << entryvalue->GetMetaDataObjectValue() << std::endl;
                }
            }
        }
        
        ++itr;
    }
    
    // find the first zero baseline image and use it for the noise estimation
    
    ::size_t iNrOfDWIs = diffusionDirections.size();
    ::size_t iFirstBaseline = std::string::npos;
    for( ::size_t iI = 0; iI < iNrOfDWIs; iI++ )
    {
        
        if( diffusionDirections[iI].GetNorm() == 0 )
        {
            iFirstBaseline = iI;
            std::cout << "First baseline found at index = " << iFirstBaseline << std::endl;
            break;
        }
        
    }
    
    if( iFirstBaseline == std::string::npos )
    {
        
        std::cout << "Did not find an explicit baseline image." << std::endl;
        std::cout << "Treating the first volume as the baseline volume." << std::endl;
        iFirstBaseline = 0;
        
    }
    
    typename ScalarImageType::SizeType indexRadiusE;
    typename ScalarImageType::SizeType indexRadiusF;
    typename ScalarImageType::SizeType indexRadiusFeatures;
    
    indexRadiusF[0] = iRadiusFiltering[0]; // radius along x
    indexRadiusF[1] = iRadiusFiltering[1]; // radius along y
    indexRadiusF[2] = iRadiusFiltering[2]; // radius along z
    
    indexRadiusE[0] = iRadiusEstimation[0]; // radius along x
    indexRadiusE[1] = iRadiusEstimation[1]; // radius along y
    indexRadiusE[2] = iRadiusEstimation[2]; // radius along z
    
    indexRadiusFeatures[0] = iRadiusFeatures[0];
    indexRadiusFeatures[1] = iRadiusFeatures[1];
    indexRadiusFeatures[2] = iRadiusFeatures[2];
    
    //================================================================================
    /** In case a mask is needed, read and apply it */
    typedef itk::Image<unsigned short,DIMENSION>   LabelImageType;
    typedef itk::ImageFileReader<LabelImageType>   LabelReaderType;
    LabelReaderType::Pointer maskReader;
    bool useMask = false;
    if( inputMask.length() > 0 ){ // We will use a mask to save computations in the background
        maskReader = LabelReaderType::New();
        maskReader->SetFileName( inputMask.c_str() );
        try{
            maskReader->Update();
            // Check if this has the appropriate size
            if( maskReader->GetOutput()->GetLargestPossibleRegion() != reader->GetOutput()->GetLargestPossibleRegion() ){
                std::cerr << "The size of the mask image does not seem to match" << std::endl;
                std::cerr << "the size of the input. Running without mask..."    << std::endl;
            }
            else
                useMask = true;
        }
        catch( itk::ExceptionObject& e ){
            std::cerr << "The filter was supposed to use a mask, but it was" << std::endl;
            std::cerr << "impossible to read it. Running without mask..."    << std::endl;
        }
    }
    
    //================================================================================
    
    typedef itk::LMMSEVectorImageFilter<DiffusionImageType, FloatDiffusionImageType> RicianLMMSEImageFilterType;
    typename RicianLMMSEImageFilterType::Pointer ricianFilter = RicianLMMSEImageFilterType::New();
    ricianFilter->SetInput( reader->GetOutput() );
    ricianFilter->SetRadius( indexRadiusF );
    // Do we need to include a mask in the computations?
    if( useMask )
        ricianFilter->SetMask( maskReader->GetOutput() );
    typename RicianLMMSEImageFilterType::GradientType grad;
    unsigned int     nDWI = 0;
    unsigned int     nBaselines = 0;
    std::vector<int> pDWI;
    std::vector<int> pBaselines;
    for( unsigned int iI = 0; iI < channels; ++iI )
    {
        float norm = diffusionDirections[iI].GetNorm();
        if( norm > 1e-3 )
        {
            grad[0] = diffusionDirections[iI][0] / norm;
            grad[1] = diffusionDirections[iI][1] / norm;
            grad[2] = diffusionDirections[iI][2] / norm;
            ricianFilter->AddGradientDirection( grad );
            ++nDWI;
            pDWI.push_back( iI );
        }
        else
        {
            ++nBaselines;
            pBaselines.push_back( iI );
        }
    }
    ricianFilter->SetNDWI( nDWI );
    ricianFilter->SetNBaselines( nBaselines );
    unsigned int* indicator = new unsigned int[nDWI];
    for( unsigned int iI = 0; iI < nDWI; ++iI )
    {
        indicator[iI] = pDWI[iI];
    }
    ricianFilter->SetDWI( indicator );
    delete[] indicator;
    indicator = new unsigned int[nBaselines];
    for( unsigned int iI = 0; iI < nBaselines; ++iI )
    {
        indicator[iI] = pBaselines[iI];
    }
    ricianFilter->SetBaselines( indicator );
    delete[] indicator;
    if( iNumNeighbors == 0 )
    {
        iNumNeighbors = nDWI;
    }
    ricianFilter->SetNeighbours( iNumNeighbors );
    
    // ======================================================================================================
    // Noise estimation is performed only if the flag "overrideNoise" IS NOT selected.
    // Otherwise, the noise level manually introduced by the user is applied
    typedef itk::Image<float, DiffusionImageType::ImageDimension>                       NoiseImageType;
    typedef itk::OtsuStatistics<DiffusionImageType, NoiseImageType>                     StatsType;
    typedef typename StatsType::Pointer                                                 StatsPointer;
    typedef itk::OtsuThreshold<NoiseImageType, NoiseImageType>                          ThresholdType;
    typedef typename ThresholdType::Pointer                                             ThresholdPointer;
    typedef itk::ThresholdToMaskImageFilter<NoiseImageType, LabelImageType>             ThresholdToMaskType;
    typedef typename ThresholdToMaskType::Pointer                                       ThresholdToMaskPointer;
    typedef itk::VectorLocalStdImageFilter<DiffusionImageType,FloatDiffusionImageType>  LocalStdType;
    typedef typename LocalStdType::Pointer                                              LocalStdPointer;
    typedef itk::ComputeRestrictedHistogram<FloatDiffusionImageType, NoiseImageType>    HistogramType;
    typedef typename HistogramType::Pointer                                             HistogramPointer;
    double sigma = itk::NumericTraits<double>::One;
    if( !overrideNoise ){
        StatsPointer           stats     = StatsType::New();
        ThresholdPointer       threshold = ThresholdType::New();
        ThresholdToMaskPointer th2mask   = ThresholdToMaskType::New();
        LocalStdPointer        localStd  = LocalStdType::New();
        HistogramPointer       histogram = HistogramType::New();
        //-------------------------------------------------
        typename StatsType::IndicatorType ind( 1 );
        ind[0] = iFirstBaseline;
        stats->SetInput( reader->GetOutput() );
        stats->SetIndicator( ind );
        stats->SetRadius( indexRadiusE );
        stats->SetChannels( channels );
        stats->SetUseNeighborhoodBaselines();
        stats->Update();
        std::cout << "The extreme values are " << stats->GetMin() << " and " << stats->GetMax() << std::endl;
        //-------------------------------------------------
        threshold->SetMin( stats->GetMin() );
        threshold->SetMax( stats->GetMax() );
        threshold->SetW( 2.0f );
        threshold->SetBins( 2048 );
        threshold->SetInput( stats->GetOutput() );
        threshold->Update();
        double th = threshold->GetThreshold();
        std::cout << "The Otsu threshold in " << th << std::endl;
        //-------------------------------------------------
        th2mask->SetInput( threshold->GetOutput() );
        th2mask->SetThreshold( th );
        th2mask->Update();
        //-------------------------------------------------
        localStd->SetInput( reader->GetOutput() );
        localStd->SetRadius( indexRadiusE );
        localStd->SetMask( th2mask->GetOutput() );
        localStd->Update();
        double meanLocalVar = localStd->GetLocalStdMean();
        std::cout << "Overall mean of local variance: " << meanLocalVar << std::endl;
        //-------------------------------------------------
        histogram->SetInput( localStd->GetOutput() );
        histogram->SetMin(  meanLocalVar / 100.0f );
        histogram->SetMax(  meanLocalVar * 3.0f  );
        histogram->SetBins( 1024 );
        histogram->Update();
        //-------------------------------------------------
        sigma  = histogram->GetMode();
        std::cout << "Estimating noise std between " << meanLocalVar / 100.0f << " and "  << meanLocalVar * 3.0f << ": " << sigma << std::endl;
    }
    else{
        sigma = noiseLevel;
    }
    // ======================================================================================================
    ricianFilter->SetSigma( sigma );
    ricianFilter->SetH( iH );
    ricianFilter->SetSetZeroBck( setZeroBck );
    ricianFilter->SetOnlyUNLM( onlyUNLM );
    ricianFilter->SetFilterOutliers( filterOutliers );
    ricianFilter->SetRadiusFeatures( indexRadiusFeatures );
    ricianFilter->Update();
    
    std::cout << "number of components per pixel" << ricianFilter->GetOutput()->GetNumberOfComponentsPerPixel()
    << std::endl;
    
    // now cast it back to diffusionimagetype
    
    typedef itk::CastImageFilter<FloatDiffusionImageType, DiffusionImageType> CastImageFilterType;
    
    typename CastImageFilterType::Pointer castImageFilter = CastImageFilterType::New();
    castImageFilter->SetInput( ricianFilter->GetOutput() );
    castImageFilter->Update();
    
    itk::MetaDataDictionary metaDataDictionary;
    metaDataDictionary = reader->GetMetaDataDictionary();
    castImageFilter->GetOutput()->SetMetaDataDictionary(metaDataDictionary);
    
    
    // let's write it out
    typedef itk::ImageFileWriter<DiffusionImageType> WriterType;
    typename WriterType::Pointer nrrdWriter = WriterType::New();
    nrrdWriter->UseInputMetaDataDictionaryOn();
    nrrdWriter->SetInput( castImageFilter->GetOutput() );
    nrrdWriter->SetFileName( outputVolume.c_str() );
    if (compressOutput)
    {
        nrrdWriter->UseCompressionOn();
    }
    else
    {
        nrrdWriter->UseCompressionOff();
    }
    try
    {
        nrrdWriter->Update();
    }
    catch( itk::ExceptionObject e )
    {
        std::cerr << argv[0] << ": exception caught !" << std::endl;
        std::cerr << e << std::endl;
        return EXIT_FAILURE;
    }
    
    std::cout << "success = " << EXIT_SUCCESS << std::endl;
    
    return EXIT_SUCCESS;
    
}

int main( int argc, char * argv[] )
{
    
    PARSE_ARGS;
    
    itk::ImageIOBase::IOPixelType     pixelType;
    itk::ImageIOBase::IOComponentType componentType;
    
    // try
    // {
    itk::GetImageType(inputVolume, pixelType, componentType);
    
    // This filter handles all types
    
    switch( componentType )
    {
#ifndef WIN32
        case itk::ImageIOBase::UCHAR:
            return DoIt( argc, argv, static_cast<unsigned char>(0) );
            break;
        case itk::ImageIOBase::CHAR:
            return DoIt( argc, argv, static_cast<char>(0) );
            break;
#endif
        case itk::ImageIOBase::USHORT:
            return DoIt( argc, argv, static_cast<unsigned short>(0) );
            break;
        case itk::ImageIOBase::SHORT:
            return DoIt( argc, argv, static_cast<short>(0) );
            break;
        case itk::ImageIOBase::UINT:
            return DoIt( argc, argv, static_cast<unsigned int>(0) );
            break;
        case itk::ImageIOBase::INT:
            return DoIt( argc, argv, static_cast<int>(0) );
            break;
#ifndef WIN32
        case itk::ImageIOBase::ULONG:
            return DoIt( argc, argv, static_cast<unsigned long>(0) );
            break;
        case itk::ImageIOBase::LONG:
            return DoIt( argc, argv, static_cast<long>(0) );
            break;
#endif
        case itk::ImageIOBase::FLOAT:
            return DoIt( argc, argv, static_cast<float>(0) );
            break;
        case itk::ImageIOBase::DOUBLE:
            std::cout << "DOUBLE type not currently supported." << std::endl;
            break;
        case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
        default:
            std::cout << "unknown component type" << std::endl;
            break;
    }
    
    // }
    
    /*catch( itk::ExceptionObject &excep)
     {
     std::cerr << argv[0] << ": exception caught !" << std::endl;
     std::cerr << excep << std::endl;
     return EXIT_FAILURE;
     }*/
    
    return EXIT_SUCCESS;
}

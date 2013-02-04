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

// ITK includes
#include <itkImageFileReader.h>
#include <itkVectorImage.h>
#include <itkImageRegionConstIterator.h>

#define DIMENSION 3

int main( int argc, char * argv[] )
{
    typedef float                                    PixelType;
    typedef itk::VectorImage<PixelType,DIMENSION>    DiffusionImageType;
    typedef itk::ImageFileReader<DiffusionImageType> FileReaderType;
    
    if( argc<3 ){
        std::cerr << "Wrong number of inputs" << std::endl;
        std::cerr << "Usage: " << argv[0] << " input output [tol=1e-6]" <<std::endl;
        return EXIT_FAILURE;
    }
    
    double tol = 1.0e-6;
    if( argc>3 )
        tol = (double)( ::atof(argv[3]) );
    
    FileReaderType::Pointer reader1 = FileReaderType::New();
    FileReaderType::Pointer reader2 = FileReaderType::New();
    
    try{
        reader1->SetFileName( argv[1] );
        reader1->Update();
        
        reader2->SetFileName( argv[2] );
        reader2->Update();
    }
    catch(...){
        std::cerr << "Unable to read the inputs" << std::endl;
        return EXIT_FAILURE;
    }
    
    if( reader1->GetOutput()->GetLargestPossibleRegion().GetSize()
           != reader2->GetOutput()->GetLargestPossibleRegion().GetSize() ){
        std::cerr << "Input 1 and input 2 have different sizes" << std::endl;
        return EXIT_FAILURE;
    }
    
    if( reader1->GetOutput()->GetVectorLength()
           != reader2->GetOutput()->GetVectorLength() ){
        std::cerr << "Input 1 and input 2 have different number of gradients" << std::endl;
        return EXIT_FAILURE; 
    }
    
    itk::ImageRegionConstIterator<DiffusionImageType> it1( reader1->GetOutput(), reader1->GetOutput()->GetLargestPossibleRegion() );
    itk::ImageRegionConstIterator<DiffusionImageType> it2( reader2->GetOutput(), reader2->GetOutput()->GetLargestPossibleRegion() );
    
    double cumdiff = itk::NumericTraits<double>::Zero;
    for( it1.GoToBegin(),it2.GoToBegin(); !it1.IsAtEnd(); ++it1,++it2 ){
        for( unsigned int k=0; k<reader1->GetOutput()->GetVectorLength(); ++k )
            cumdiff += vcl_abs( it1.Get()[k] - it2.Get()[k] );
    }
    
    if( cumdiff > tol ){
        std::cerr << "The error (" << cumdiff << ") is above the tolerance (" << tol << ")" << std::endl;
        return EXIT_FAILURE;
    }
    
    return EXIT_SUCCESS;
    
}

cmake_minimum_required(VERSION 3.1)
project(RGBDTracking)

set(RGBDTRACKING_MAJOR_VERSION 0)
set(RGBDTRACKING_MINOR_VERSION 1)
set(RGBDTRACKING_VERSION ${RGBDTRACKING_MAJOR_VERSION}.${RGBDTRACKING_MINOR_VERSION})

find_package(OpenCV)
option(WITH_CUDA FALSE)
set(HAVECUDA FALSE)
IF(WITH_CUDA)
find_package(CUDA)
IF(CUDA_FOUND)
  MESSAGE(STATUS "CUDA found")
  set(HAVECUDA TRUE)
ELSE(CUDA_FOUND)
  MESSAGE(STATUS "CUDA not found")
ENDIF(CUDA_FOUND)
ELSE(WITH_CUDA)
MESSAGE(STATUS "Without CUDA") 
ENDIF(WITH_CUDA)

find_package(image QUIET)

FIND_PACKAGE(VISP REQUIRED)
IF(VISP_FOUND)
  INCLUDE(${VISP_USE_FILE})
ENDIF(VISP_FOUND)
LINK_DIRECTORIES(${VISP_LINK_DIRECTORIES})

message(STATUS ${VISP_INCLUDE_DIRS})

find_package(PCL 1.9 REQUIRED COMPONENTS)
link_directories(${PCL_LIBRARY_DIRS})
include_directories(BEFORE ${PCL_INCLUDE_DIRS})
add_definitions(${PCL_DEFINITIONS})

find_package(Boost COMPONENTS system filesystem REQUIRED)

set(ALL_LIBRARIES ${OpenCV_LIBS} freeimage ${FREEGLUT_LIBRARY} ${GLEW_LIBRARY} ${OPENGL_LIBRARIES} ${PCL_LIBRARIES} ${VISP_LIBRARIES})

set(HEADER_FILES
	RGBDDataProcessing.h
	RGBDDataProcessing.inl
	DataIO.h
	DataIO.inl
	MeshProcessing.h
	MeshProcessing.inl
	ClosestPoint.h
	ClosestPoint.inl
	ClosestPointForceField.h
#        FeatureMatchingForceField.h
	RenderTextureAR.h
	RenderTextureAR.inl
	RegistrationRigid.h
	ImageConverter.h
	ImageConverter.inl
	#DataGeneration.h
   	segmentation.h
#	KalmanFilter.h
#	ccd.h
#	p_helper.h
        RenderingManager.h
#	ShapeSpringsForceField.h
	VirtualCamera.h
	VirtualCamera.inl
#        DataSource.hpp
#        serialization.h
#        ZMQCommunication.h
)

set(SOURCE_FILES
	initRGBDTracking.cpp
	RegistrationRigid.cpp
	ImageConverter.cpp
	#DataGeneration.cpp
    	segmentation.cpp
#	KalmanFilter.cpp
#	ccd.cpp
	RGBDDataProcessing.cpp
	DataIO.cpp
	MeshProcessing.cpp
	RenderTextureAR.cpp
	ClosestPoint.cpp
	ClosestPointForceField.cpp
#        FeatureMatchingForceField.cpp
	RenderingManager.cpp
#	ShapeSpringsForceField.cpp
	VirtualCamera.cpp
#        ZMQCommunication.cpp
)

set(README_FILES rgbdtracking.txt)

if (HAVECUDA)

list(APPEND HEADER_FILES cudaSegmentation.h)
list(APPEND SOURCE_FILES cudaSegmentation.cpp)

set(CUDA_SOURCES
GrabcutGMM.cu 
GrabcutHistogram.cu 
GrabcutUtil.cu
#emicp.cu
)

set(${PROJECT_NAME}_INCLUDE_EXPORTS
  ${OpenCV_INCLUDE_DIRS}
  ${VISP_INCLUDE_DIR}
${VISP_INCLUDE_DIRS}
/usr/local/cuda-7.0/include/
/usr/local/cuda-7.0/samples/common/
/usr/local/cuda-7.0/samples/common/inc/
/usr/local/cuda-7.0/samples/7_CUDALibraries/common/UtilNPP/
)

# build options
set(GENCODE -gencode=arch=compute_30,code=sm_30 -gencode=arch=compute_35,code=compute_35)
set(GENCODE ${GENCODE} -gencode=arch=compute_20,code=sm_20)
set(GENCODE ${GENCODE} -gencode=arch=compute_11,code=sm_11)

# debug options
set(CUDA_NVCC_FLAGS_DEBUG "-G;-g")
set(CUDA_NVCC_FLAGS_RELEASE "-O2")
set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} \"-D_FORCE_INLINES\" ")

#Give current directory to nvcc compiler
cuda_include_directories(${CMAKE_CURRENT_SOURCE_DIR})
cuda_include_directories(${SOFA_FRAMEWORK_DIR})
cuda_include_directories(${SOFA_EXTLIBS_DIR})
cuda_include_directories(${SOFA_INC_DIR})
else(HAVECUDA)
set(${PROJECT_NAME}_INCLUDE_EXPORTS
  ${OpenCV_INCLUDE_DIRS}
  ${VISP_INCLUDE_DIR}
${VISP_INCLUDE_DIRS}
)
#Give current directory to nvcc compiler
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${SOFA_FRAMEWORK_DIR})
include_directories(${SOFA_EXTLIBS_DIR})
include_directories(${SOFA_INC_DIR})
endif(HAVECUDA)

set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

set_directory_properties( PROPERTIES COMPILE_DEFINITIONS "" )


## Specify additional locations of header files
## Your package locations should be listed before other locations
include_directories(
INCLUDE_DIRS ${${PROJECT_NAME}_INCLUDE_EXPORTS}
)

## RGBDTracking/config.h
configure_file(config.h.in "${CMAKE_BINARY_DIR}/include/RGBDTracking/config.h")
install(FILES "${CMAKE_BINARY_DIR}/include/RGBDTracking/config.h" DESTINATION "include/RGBDTracking")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

# Config files and install rules for pythons scripts
sofa_set_python_directory(RGBDTracking "python")


#Get nvcc compiler to build the library
if(HAVECUDA)
cuda_add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} ${CUDA_SOURCES} ${README_FILES} OPTIONS ${VERBOSE_PTXAS})
else(HAVECUDA)
add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} ${README_FILES})
endif(HAVECUDA)

if(SOFA-CUDA_CUBLAS)
    CUDA_ADD_CUBLAS_TO_TARGET(${PROJECT_NAME})
    find_package(CUDASparse REQUIRED)
    #message(STATUS CUDA_SPARSE_LIBRARY=${CUDA_SPARSE_LIBRARY})
    target_link_libraries(${PROJECT_NAME} ${CUDA_SPARSE_LIBRARY} ${GENCODE})
endif()

target_link_libraries(${PROJECT_NAME} ${ALL_LIBRARIES} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} image SofaGuiQt SofaMeshCollision SofaMiscCollision SofaBaseCollision SofaGuiCommon SofaBaseVisual SofaExporter SofaLoader SofaMiscForceField SofaGeneralEngine)

if (HAVECUDA)
target_link_libraries(${PROJECT_NAME} ${GENCODE} ${CUDA_LIBRARIES} ${CUDA_nppi_LIBRARY})
endif(HAVECUDA)

if(NOT ${SOFA_NO_OPENGL})
    target_link_libraries(${PROJECT_NAME} SofaOpenglVisual)
endif()

target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>")
target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>")
target_include_directories(${PROJECT_NAME} PUBLIC "$<INSTALL_INTERFACE:include>")
set_target_properties(${PROJECT_NAME} PROPERTIES VERSION ${RGBDTRACKING_VERSION})
set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "-DSOFA_BUILD_RGBDTRACKING")
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${HEADER_FILES}")

set(RGBDTRACKING_COMPILER_FLAGS "-DSOFA_BUILD_RGBDTRACKING")

set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "${RGBDTRACKING_COMPILER_FLAGS}")


## Install rules for the library; CMake package configurations files
sofa_create_package(RGBDTracking ${RGBDTRACKING_VERSION} RGBDTracking RGBDTracking)


install(DIRECTORY examples/ DESTINATION share/sofa/plugins/${PROJECT_NAME})


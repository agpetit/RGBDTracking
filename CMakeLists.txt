cmake_minimum_required(VERSION 3.1)
project(RGBDTracking)

set(RGBDTRACKING_MAJOR_VERSION 0)
set(RGBDTRACKING_MINOR_VERSION 1)
set(RGBDTRACKING_VERSION ${RGBDTRACKING_MAJOR_VERSION}.${RGBDTRACKING_MINOR_VERSION})

find_package(OpenCV)

find_package(image QUIET)

FIND_PACKAGE(VISP)
IF(VISP_FOUND)
  INCLUDE(${VISP_USE_FILE})
ENDIF(VISP_FOUND)
LINK_DIRECTORIES(${VISP_LINK_DIRECTORIES})

message(STATUS ${VISP_INCLUDE_DIRS})

find_package(PCL 1.8 REQUIRED COMPONENTS)
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
	RenderTextureAR.h
	RenderTextureAR.inl
	RegistrationForceFieldCam.h
	RegistrationRigid.h
	ImageConverter.h
	ImageConverter.inl
	#DataGeneration.h
   	segmentation.h
	KalmanFilter.h
	ccd.h
	p_helper.h
        RenderingManager.h
	ShapeSpringsForceField.h

)

set(SOURCE_FILES
	initRGBDTracking.cpp
        RegistrationForceFieldCam.cpp
	RegistrationRigid.cpp
	ImageConverter.cpp
	#DataGeneration.cpp
    	segmentation.cpp
	KalmanFilter.cpp
	ccd.cpp
	RGBDDataProcessing.cpp
	DataIO.cpp
	MeshProcessing.cpp
	RenderTextureAR.cpp
	ClosestPoint.cpp
	RenderingManager.cpp
	ShapeSpringsForceField.cpp
)

set(README_FILES rgbdtracking.txt)

set(${PROJECT_NAME}_INCLUDE_EXPORTS
  ${OpenCV_INCLUDE_DIRS}
  ${VISP_INCLUDE_DIR}
${VISP_INCLUDE_DIRS}
${POWERBALLCTRL_INCLUDE_DIR}
)

if (MSVC)
  set(CMAKE_CXX_FLAGS_DEBUG "/DEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE "/O2")
else()
  set(CMAKE_CXX_FLAGS_DEBUG "-g")
  set(CMAKE_CXX_FLAGS_RELEASE "-O2")
endif()


set_directory_properties( PROPERTIES COMPILE_DEFINITIONS "" )

#Give current directory to nvcc compiler
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${SOFA_FRAMEWORK_DIR})
include_directories(${SOFA_EXTLIBS_DIR})
include_directories(${SOFA_INC_DIR})

## Specify additional locations of header files
## Your package locations should be listed before other locations
include_directories(
INCLUDE_DIRS ${${PROJECT_NAME}_INCLUDE_EXPORTS}
)

## RGBDTracking/config.h
configure_file(config.h.in "${CMAKE_BINARY_DIR}/include/RGBDTracking/config.h")
install(FILES "${CMAKE_BINARY_DIR}/include/RGBDTracking/config.h" DESTINATION "include/RGBDTracking")

# Config files and install rules for pythons scripts
sofa_set_python_directory(RGBDTracking "python")

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} ${README_FILES})

target_link_libraries(${PROJECT_NAME} ${ALL_LIBRARIES} ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY} image SofaGuiQt SofaMeshCollision SofaMiscCollision SofaBaseCollision SofaGuiCommon SofaBaseVisual SofaExporter SofaLoader SofaMiscForceField SofaGeneralEngine)

if(NOT ${SOFA_NO_OPENGL})
    target_link_libraries(${PROJECT_NAME} SofaOpenglVisual)
endif()

target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>")
target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>")
set_target_properties(${PROJECT_NAME} PROPERTIES VERSION ${RGBDTRACKING_VERSION})
set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "-DSOFA_BUILD_RGBDTRACKING")
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${HEADER_FILES}")

## Install rules for the library; CMake package configurations files
sofa_create_package(RGBDTracking ${RGBDTRACKING_VERSION} RGBDTracking RGBDTracking)
install(DIRECTORY sofa
       DESTINATION include/${PROJECT_NAME}
       FILES_MATCHING PATTERN "*.h" PATTERN "*hpp" PATTERN "*.inl"
)



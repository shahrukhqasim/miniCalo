# $Id: CMakeLists.txt 104871 2017-06-23 14:17:13Z gcosmo $

#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
MESSAGE(${CMAKE_VERSION})
project(B4a)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
#option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
#if(WITH_GEANT4_UIVIS)
#  MESSAGE(WITH UI TOOLS)
#option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
find_package(Geant4 REQUIRED ui_all vis_all)
#else()
#  find_package(Geant4 REQUIRED)
#endif()
#
#
##find_package(root REQUIRED)
#get_cmake_property(_variableNames VARIABLES)
#list (SORT _variableNames)
#foreach (_variableName ${_variableNames})
#  message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()


find_package(Boost 1.56 REQUIRED COMPONENTS json)

include_directories(${Boost_INCLUDE_DIRS})


#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(/Users/shahrukhqasim/Workspace/NextCal/pythia8306/include)
include_directories(/Users/shahrukhqasim/Workspace/NextCal/fastjet-3.4.0/fastjet-install/include)
include_directories(/usr/local/include/root/)
MESSAGE(${Geant4_INCLUDE_DIRS})
include_directories(${Geant4_INCLUDE_DIRS})


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
#file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
#file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)




#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(exampleB4a exampleB4a.cc
        src/B4DetectorConstruction.cc
        src/B4JetGeneratorAction.cc
        src/B4PrimaryGeneratorAction.cc
        src/B4aActionInitialization.cc
        src/B4aSteppingAction.cc
        src/B4PartGeneratorBase.cc
        src/B4RunAction.cc
        src/B4aEventAction.cc
        src/sensorContainer.cc)


add_executable(visualize src/visualize.cc
        src/B4DetectorConstruction.cc
        src/B4JetGeneratorAction.cc
        src/B4PrimaryGeneratorAction.cc
        src/B4aActionInitialization.cc
        src/B4aSteppingAction.cc
        src/B4PartGeneratorBase.cc
        src/B4RunAction.cc
        src/B4aEventAction.cc
        src/sensorContainer.cc)



add_executable(json_test src/json_reading_test.cc)

MESSAGE(${Geant4_LIBRARIES})
target_link_libraries(exampleB4a Boost::json ${Geant4_LIBRARIES} /Users/shahrukhqasim/Workspace/NextCal/pythia8306/lib/libpythia8.dylib /Users/shahrukhqasim/Workspace/NextCal/fastjet-3.4.0/fastjet-install/lib/libfastjet.dylib /usr/local/Cellar/root/6.24.04/lib/root/libCore.so;/usr/local/Cellar/root/6.24.04/lib/root/libImt.so;/usr/local/Cellar/root/6.24.04/lib/root/libRIO.so;/usr/local/Cellar/root/6.24.04/lib/root/libNet.so;/usr/local/Cellar/root/6.24.04/lib/root/libHist.so;/usr/local/Cellar/root/6.24.04/lib/root/libGraf.so;/usr/local/Cellar/root/6.24.04/lib/root/libGraf3d.so;/usr/local/Cellar/root/6.24.04/lib/root/libGpad.so;/usr/local/Cellar/root/6.24.04/lib/root/libROOTDataFrame.so;/usr/local/Cellar/root/6.24.04/lib/root/libTree.so;/usr/local/Cellar/root/6.24.04/lib/root/libTreePlayer.so;/usr/local/Cellar/root/6.24.04/lib/root/libRint.so;/usr/local/Cellar/root/6.24.04/lib/root/libPostscript.so;/usr/local/Cellar/root/6.24.04/lib/root/libMatrix.so;/usr/local/Cellar/root/6.24.04/lib/root/libPhysics.so;/usr/local/Cellar/root/6.24.04/lib/root/libMathCore.so;/usr/local/Cellar/root/6.24.04/lib/root/libThread.so;/usr/local/Cellar/root/6.24.04/lib/root/libMultiProc.so;/usr/local/Cellar/root/6.24.04/lib/root/libROOTVecOps.so)
target_link_libraries(visualize Boost::json ${Geant4_LIBRARIES} /Users/shahrukhqasim/Workspace/NextCal/pythia8306/lib/libpythia8.dylib /Users/shahrukhqasim/Workspace/NextCal/fastjet-3.4.0/fastjet-install/lib/libfastjet.dylib /usr/local/Cellar/root/6.24.04/lib/root/libCore.so;/usr/local/Cellar/root/6.24.04/lib/root/libImt.so;/usr/local/Cellar/root/6.24.04/lib/root/libRIO.so;/usr/local/Cellar/root/6.24.04/lib/root/libNet.so;/usr/local/Cellar/root/6.24.04/lib/root/libHist.so;/usr/local/Cellar/root/6.24.04/lib/root/libGraf.so;/usr/local/Cellar/root/6.24.04/lib/root/libGraf3d.so;/usr/local/Cellar/root/6.24.04/lib/root/libGpad.so;/usr/local/Cellar/root/6.24.04/lib/root/libROOTDataFrame.so;/usr/local/Cellar/root/6.24.04/lib/root/libTree.so;/usr/local/Cellar/root/6.24.04/lib/root/libTreePlayer.so;/usr/local/Cellar/root/6.24.04/lib/root/libRint.so;/usr/local/Cellar/root/6.24.04/lib/root/libPostscript.so;/usr/local/Cellar/root/6.24.04/lib/root/libMatrix.so;/usr/local/Cellar/root/6.24.04/lib/root/libPhysics.so;/usr/local/Cellar/root/6.24.04/lib/root/libMathCore.so;/usr/local/Cellar/root/6.24.04/lib/root/libThread.so;/usr/local/Cellar/root/6.24.04/lib/root/libMultiProc.so;/usr/local/Cellar/root/6.24.04/lib/root/libROOTVecOps.so)

target_link_libraries(visualize Boost::json)
#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B4a. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(EXAMPLEB4A_SCRIPTS
  exampleB4a.out
  exampleB4.in
  gui.mac
  init_vis.mac
  plotHisto.C
  run1.mac
  run2.mac
  vis.mac
  batchrun.mac
  condorSubmit
  batchscript.sh
  env.sh
  )

foreach(_script ${EXAMPLEB4A_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
#install(TARGETS exampleB4a DESTINATION bin)


#add_executable(helloworld src/helloworld.cc)

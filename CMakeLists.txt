# $Id: CMakeLists.txt 104871 2017-06-23 14:17:13Z gcosmo $

#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(B4a)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
#option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
#if(WITH_GEANT4_UIVIS)
find_package(Geant4 REQUIRED ui_all vis_all)
#else()
#  find_package(Geant4 REQUIRED)
#endif()


include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(/Users/shahrukhqasim/Workspace/NextCal/pythia8306/include)
include_directories(/Users/shahrukhqasim/Workspace/NextCal/fastjet-3.4.0/fastjet-install/include)
include_directories(/usr/local/include/root/)
include_directories(g4h/include)
MESSAGE(${Geant4_INCLUDE_DIRS})
include_directories(${Geant4_INCLUDE_DIRS})

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(exampleB4a exampleB4a.cc ${sources} ${headers})

#add_library()


set(plibs /Users/shahrukhqasim/Workspace/NextCal/pythia8306/lib/libpythia8.dylib;/Users/shahrukhqasim/Workspace/NextCal/fastjet-3.4.0/fastjet-install/lib/libfastjet.dylib;/usr/local/Cellar/root/6.24.04/lib/root/libCore.so;/usr/local/Cellar/root/6.24.04/lib/root/libImt.so;/usr/local/Cellar/root/6.24.04/lib/root/libRIO.so;/usr/local/Cellar/root/6.24.04/lib/root/libNet.so;/usr/local/Cellar/root/6.24.04/lib/root/libHist.so;/usr/local/Cellar/root/6.24.04/lib/root/libGraf.so;/usr/local/Cellar/root/6.24.04/lib/root/libGraf3d.so;/usr/local/Cellar/root/6.24.04/lib/root/libGpad.so;/usr/local/Cellar/root/6.24.04/lib/root/libROOTDataFrame.so;/usr/local/Cellar/root/6.24.04/lib/root/libTree.so;/usr/local/Cellar/root/6.24.04/lib/root/libTreePlayer.so;/usr/local/Cellar/root/6.24.04/lib/root/libRint.so;/usr/local/Cellar/root/6.24.04/lib/root/libPostscript.so;/usr/local/Cellar/root/6.24.04/lib/root/libMatrix.so;/usr/local/Cellar/root/6.24.04/lib/root/libPhysics.so;/usr/local/Cellar/root/6.24.04/lib/root/libMathCore.so;/usr/local/Cellar/root/6.24.04/lib/root/libThread.so;/usr/local/Cellar/root/6.24.04/lib/root/libMultiProc.so;/usr/local/Cellar/root/6.24.04/lib/root/libROOTVecOps.so)
#set(customgeantlibs /Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4FR.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4GMocren.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4OpenGL.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4OpenInventor.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4RayTracer.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4Tree.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4VRML.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4analysis.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4clhep.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4digits_hits.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4error_propagation.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4event.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4geometry.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4gl2ps.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4global.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4graphics_reps.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4intercoms.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4interfaces.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4materials.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4modeling.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4parmodels.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4particles.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4persistency.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4physicslists.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4processes.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4ptl.0.0.2.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4ptl.0.0.2.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4ptl.0.0.2.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4readout.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4run.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4tasking.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4track.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4tracking.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4visHepRep.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4visQt3D.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4visXXX.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4vis_management.dylib;/Users/shahrukhqasim/Workspace/NextCal/geant4.10.07.p02/make/BuildProducts/lib/libG4zlib.dylib)
add_library(g4libs g4h/src/CustomPDGParser.cc
        g4h/src/ExRhadDetectorMessenger.cc
        g4h/src/ExRhadPrimaryGeneratorMessenger.cc
        g4h/src/G4ProcessHelper.cc
        g4h/src/MyStackingAction.cc
        g4h/src/CustomParticle.cc
        g4h/src/ExRhadEventAction.cc
        g4h/src/ExRhadRunAction.cc
        g4h/src/HadronicProcessHelper.cc
        g4h/src/MyTrackerHit.cc
        g4h/src/CustomParticleFactory.cc
        g4h/src/ExRhadEventActionMessenger.cc
        g4h/src/ExRhadSteppingAction.cc
        g4h/src/MyAnalysis.cc
        g4h/src/MyTrackingAction.cc
        g4h/src/Decay3Body.cc
        g4h/src/ExRhadPhysicsList.cc
        g4h/src/ExRhadSteppingVerbose.cc
        g4h/src/MyG4ReactionDynamics.cc
        g4h/src/ToyModelHadronicProcess.cc
        g4h/src/ExRhadDetectorConstruction.cc
        g4h/src/ExRhadPrimaryGeneratorAction.cc
        g4h/src/FullModelHadronicProcess.cc
        g4h/src/MySensitiveTracker.cc)


MESSAGE(HELLO)
MESSAGE(${Geant4_LIBRARIES})
MESSAGE(HELLO2)

target_link_libraries(g4libs ${Geant4_LIBRARIES} ${plibs})
target_link_libraries(exampleB4a ${Geant4_LIBRARIES} ${plibs} g4libs)

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
  osxenv.sh
  )

foreach(_script ${EXAMPLEB4A_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

set(PARTICLE_FILES
  particles.txt
  decays.txt
  PhysicsConfiguration.txt
  ProcessList.txt
  )

set(G4RHADRONS g4h)

foreach(_file ${PARTICLE_FILES})
  configure_file(
    ${G4RHADRONS}/${_file}
    ${PROJECT_BINARY_DIR}/${_file}
    COPYONLY
    )
endforeach()

file(COPY ${G4RHADRONS}/data DESTINATION ${PROJECT_BINARY_DIR}/)


#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
#install(TARGETS exampleB4a DESTINATION bin)

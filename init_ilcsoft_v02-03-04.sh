export ILCSOFT=/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04

# -------------------------------------------------------------------- ---

# ---  Use the same compiler and python as used for the installation   ---

# -------------------------------------------------------------------- ---
. /cvmfs/sft.cern.ch/lcg/releases/gcc/13.1.0-b3d18/x86_64-el9/setup.sh

export PATH=/cvmfs/sft.cern.ch/lcg/releases/Python/3.11.9-2924c/x86_64-el9-gcc13-opt/bin:${PATH}
export LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/releases/Python/3.11.9-2924c/x86_64-el9-gcc13-opt/lib:${LD_LIBRARY_PATH}

export PYTHONPATH=/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python/examples:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/podio/v01-02/install/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/podio/v01-02/install/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python/examples:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/podio/v01-02/install/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python/examples:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/podio/v01-02/install/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python/examples:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/podio/v01-02/install/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python/examples:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/podio/v01-02/install/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python/examples:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/podio/v01-02/install/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python/examples:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/podio/v01-02/install/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python/examples:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python/examples:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/podio/v01-02/install/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python/examples:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/podio/v01-02/install/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python/examples:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/podio/v01-02/install/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python/examples:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/podio/v01-02/install/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python/examples:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/podio/v01-02/install/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python/examples:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python/examples:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python/examples:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python/examples:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/podio/v01-02/install/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python/examples:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/podio/v01-02/install/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python/examples:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/podio/v01-02/install/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python/examples:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/podio/v01-02/install/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python/examples:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/podio/v01-02/install/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python/examples:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python/examples:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/podio/v01-02/install/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python/examples:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/podio/v01-02/install/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python/examples:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/podio/v01-02/install/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python/examples:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/podio/v01-02/install/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python/examples:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/podio/v01-02/install/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python/examples:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/podio/v01-02/install/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python/examples:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/podio/v01-02/install/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python/examples:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/podio/v01-02/install/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python/examples:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/podio/v01-02/install/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/podio/v01-02/install/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python/examples:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python/examples:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/python/examples:/cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/lib/python3.11/site-packages

export CXX=g++
export CC=gcc
export GIT_EXEC_PATH=/cvmfs/sft.cern.ch/lcg/releases/git/2.29.2-daa17/x86_64-el9-gcc13-opt/libexec/git-core

#--------------------------------------------------------------------------------
#     LCCD
#--------------------------------------------------------------------------------
export LCCD="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lccd/v01-05-03"
# --- additional LCCD commands ------- 


#--------------------------------------------------------------------------------
#     CondDBMySQL
#--------------------------------------------------------------------------------
export CondDBMySQL="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/CondDBMySQL/CondDBMySQL_ILC-0-9-7"
export COND_DB_DEBUGLOG="/dev/stdout"
# --- additional CondDBMySQL commands ------- 
export LD_LIBRARY_PATH="$CondDBMySQL/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     LCIO
#--------------------------------------------------------------------------------
export LCIO="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05"
# --- additional LCIO commands ------- 
export PATH="$LCIO/tools:$LCIO/bin:$PATH"
export LD_LIBRARY_PATH="$LCIO/lib:$LCIO/lib64:$LD_LIBRARY_PATH"
export PYTHONPATH="$LCIO/python:$LCIO/python/examples:$PYTHONPATH"
export ROOT_INCLUDE_PATH="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcio/v02-22-05/include:$ROOT_INCLUDE_PATH"
# --- additional ROOT commands ------- 
test -r /afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/root/6.32.12/bin/thisroot.sh && . /afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/root/6.32.12/bin/thisroot.sh


#--------------------------------------------------------------------------------
#     CMake
#--------------------------------------------------------------------------------
# --- additional CMake commands ------- 
export PATH="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/CMake/3.31.6/bin:$PATH"


#--------------------------------------------------------------------------------
#     ILCUTIL
#--------------------------------------------------------------------------------
export ilcutil="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/ilcutil/v01-07-03"
# --- additional ILCUTIL commands ------- 
export LD_LIBRARY_PATH="$ilcutil/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     Marlin
#--------------------------------------------------------------------------------
export MARLIN="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/Marlin/v01-19-05"
# --- additional Marlin commands ------- 
export PATH="$MARLIN/bin:$PATH"
export MARLIN_DLL="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/LCFIVertex/v00-09/install/lib/libLCFIVertexProcessors.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/LCFIPlus/v00-11/install/lib/libLCFIPlus.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/MarlinDD4hep/v00-06-02/install/lib/libMarlinDD4hep.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/DDMarlinPandora/v00-13/install/lib/libDDMarlinPandora.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/MarlinKinfit/v00-06-01/install/lib/libMarlinKinfit.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/MarlinReco/v01-36-02/install/lib/libMarlinReco.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/PandoraAnalysis/v02-00-01/install/lib/libPandoraAnalysis.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/CEDViewer/v01-20/install/lib/libCEDViewer.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/Overlay/v00-23/install/lib/libOverlay.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/MarlinFastJet/v00-05-03/install/lib/libMarlinFastJet.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/LCTuple/v01-14/install/lib/libLCTuple.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/MarlinTrkProcessors/v02-12-07/install/lib/libMarlinTrkProcessors.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/MarlinKinfitProcessors/v00-05/install/lib/libMarlinKinfitProcessors.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/ILDPerformance/v01-12/install/lib/libILDPerformance.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/Clupatra/v01-03-01/install/lib/libClupatra.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/FCalClusterer/v01-00-06/install/lib/libFCalClusterer.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/ForwardTracking/v01-14-02/install/lib/libForwardTracking.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/ConformalTracking/v01-12/install/lib/libConformalTracking.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/LICH/v00-01/install/lib/libLICH.so:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/Garlic/v03-01/install/lib/libGarlic.so:$MARLIN_DLL"


#--------------------------------------------------------------------------------
#     CLHEP
#--------------------------------------------------------------------------------
export CLHEP="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/CLHEP/2.4.7.1"
export CLHEP_BASE_DIR="$CLHEP"
export CLHEP_INCLUDE_DIR="$CLHEP/include"
# --- additional CLHEP commands ------- 
export PATH="$CLHEP_BASE_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$CLHEP_BASE_DIR/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     RAIDA
#--------------------------------------------------------------------------------
export RAIDA_HOME="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/RAIDA/v01-11"
# --- additional RAIDA commands ------- 
export PATH="$RAIDA_HOME/bin:$PATH"


#--------------------------------------------------------------------------------
#     GEAR
#--------------------------------------------------------------------------------
export GEAR="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/gear/v01-09-05"
# --- additional GEAR commands ------- 
export PATH="$GEAR/tools:$GEAR/bin:$PATH"
export LD_LIBRARY_PATH="$GEAR/lib:$LD_LIBRARY_PATH"
# --- additional MarlinDD4hep commands ------- 


#--------------------------------------------------------------------------------
#     DD4hep
#--------------------------------------------------------------------------------
export DD4HEP="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/DD4hep/v01-31"
export DD4HEP_ENVINIT="${DD4HEP}/bin/thisdd4hep.sh"
export DD4hepExamples="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/DD4hep/v01-31/examples"
# --- additional DD4hep commands ------- 
test -r ${DD4HEP_ENVINIT} && . ${DD4HEP_ENVINIT}


#--------------------------------------------------------------------------------
#     Geant4
#--------------------------------------------------------------------------------
export G4INSTALL="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/geant4/11.3.1"
export G4SYSTEM="unknown-g++"
export G4ENV_INIT="$G4INSTALL/bin/geant4.sh"
# --- additional Geant4 commands ------- 
test -r ${G4ENV_INIT} && { cd $(dirname ${G4ENV_INIT}) ; . ./$(basename ${G4ENV_INIT}) ; cd $OLDPWD ; }


#--------------------------------------------------------------------------------
#     XercesC
#--------------------------------------------------------------------------------
export XercesC_HOME="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/xercesc/v3.2.3"
# --- additional XercesC commands ------- 
export PATH="$XercesC_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$XercesC_HOME/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     Boost
#--------------------------------------------------------------------------------
export BOOST_ROOT="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/boost/1.87.0"
# --- additional Boost commands ------- 


#--------------------------------------------------------------------------------
#     edm4hep
#--------------------------------------------------------------------------------
export EDM4HEP_DIR="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/edm4hep/v00-99-01/install"
# --- additional edm4hep commands ------- 
export LD_LIBRARY_PATH="$EDM4HEP_DIR/lib:$EDM4HEP_DIR/lib64:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     podio
#--------------------------------------------------------------------------------
export podio_DIR="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/podio/v01-02/install"
# --- additional podio commands ------- 
export PATH="$podio_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$podio_DIR/lib:$podio_DIR/lib64:$LD_LIBRARY_PATH"
export PYTHONPATH="$podio_DIR/python:$PYTHONPATH"
export ROOT_INCLUDE_PATH="$podio_DIR/include:$ROOT_INCLUDE_PATH"


#--------------------------------------------------------------------------------
#     SIO
#--------------------------------------------------------------------------------
export SIO_DIR="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/sio/v00-02"
# --- additional SIO commands ------- 
export PATH="$SIO_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$SIO_DIR/lib:$SIO_DIR/lib64:$LD_LIBRARY_PATH"
# --- additional DDKalTest commands ------- 


#--------------------------------------------------------------------------------
#     KalTest
#--------------------------------------------------------------------------------
export KALTEST="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/KalTest/v02-05-02"
# --- additional KalTest commands ------- 
export LD_LIBRARY_PATH="$KALTEST/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     aidaTT
#--------------------------------------------------------------------------------
export AIDATT="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/aidaTT/v00-10"
# --- additional aidaTT commands ------- 
export PATH="$AIDATT/bin:$PATH"
export LD_LIBRARY_PATH="$AIDATT/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     GSL
#--------------------------------------------------------------------------------
export GSL_HOME="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/gsl/2.7"
# --- additional GSL commands ------- 
export PATH="$GSL_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$GSL_HOME/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     GBL
#--------------------------------------------------------------------------------
export GBL="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/GBL/V02-02-01"
# --- additional GBL commands ------- 
export LD_LIBRARY_PATH="$GBL/lib:$LD_LIBRARY_PATH"
# --- additional Eigen commands ------- 


#--------------------------------------------------------------------------------
#     MarlinUtil
#--------------------------------------------------------------------------------
# --- additional MarlinUtil commands ------- 
export LD_LIBRARY_PATH="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/MarlinUtil/v01-18-02/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     CED
#--------------------------------------------------------------------------------
# --- additional CED commands ------- 
export PATH="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/CED/v01-10/bin:$PATH"
export LD_LIBRARY_PATH="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/CED/v01-10/lib:$LD_LIBRARY_PATH"
# --- additional DDMarlinPandora commands ------- 


#--------------------------------------------------------------------------------
#     PandoraPFANew
#--------------------------------------------------------------------------------
export PANDORAPFANEW="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/PandoraPFANew/v04-15-00"
# --- additional PandoraPFANew commands ------- 
export LD_LIBRARY_PATH="$PANDORAPFANEW/lib:$LD_LIBRARY_PATH"
# --- additional MarlinTrk commands ------- 
# --- additional KalDet commands ------- 
# --- additional MarlinReco commands ------- 
# --- additional MarlinKinfit commands ------- 


#--------------------------------------------------------------------------------
#     PandoraAnalysis
#--------------------------------------------------------------------------------
export PANDORA_ANALYSIS_DIR="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/PandoraAnalysis/v02-00-01"
# --- additional PandoraAnalysis commands ------- 
export PATH="$PANDORA_ANALYSIS_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$PANDORA_ANALYSIS_DIR/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     LCFIVertex
#--------------------------------------------------------------------------------
# --- additional LCFIVertex commands ------- 
export LD_LIBRARY_PATH="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/LCFIVertex/v00-09install/lib:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     CEDViewer
#--------------------------------------------------------------------------------
# --- additional CEDViewer commands ------- 
export PATH="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/CEDViewer/v01-20/bin:$PATH"
export LD_LIBRARY_PATH="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/CEDViewer/v01-20/lib:$LD_LIBRARY_PATH"
# --- additional Overlay commands ------- 
# --- additional MarlinFastJet commands ------- 


#--------------------------------------------------------------------------------
#     FastJet
#--------------------------------------------------------------------------------
export FastJet_HOME="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/FastJet/3.4.3"
# --- additional FastJet commands ------- 
export PATH="$FastJet_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$FastJet_HOME/lib:$LD_LIBRARY_PATH"
# --- additional LCTuple commands ------- 
# --- additional KiTrack commands ------- 
# --- additional KiTrackMarlin commands ------- 
# --- additional MarlinTrkProcessors commands ------- 
# --- additional MarlinKinfitProcessors commands ------- 
# --- additional ILDPerformance commands ------- 
# --- additional Clupatra commands ------- 
# --- additional Physsim commands ------- 
# --- additional FCalClusterer commands ------- 


#--------------------------------------------------------------------------------
#     LCFIPlus
#--------------------------------------------------------------------------------
# --- additional LCFIPlus commands ------- 
export LD_LIBRARY_PATH="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/LCFIPlus/v00-11/install/lib:$LD_LIBRARY_PATH"
export ROOT_INCLUDE_PATH="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/LCFIPlus/v00-11/include:$ROOT_INCLUDE_PATH"
# --- additional ForwardTracking commands ------- 
# --- additional ConformalTracking commands ------- 
# --- additional LICH commands ------- 
# --- additional Garlic commands ------- 


#--------------------------------------------------------------------------------
#     lcgeo
#--------------------------------------------------------------------------------
export lcgeo_DIR="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcgeo/v00-21"
export lcgeo_ENVINIT="${lcgeo_DIR}/bin/thislcgeo.sh"
export k4geo_DIR="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/lcgeo/v00-21"
# --- additional lcgeo commands ------- 
test -r ${lcgeo_ENVINIT} && . ${lcgeo_ENVINIT}
test -r ${k4geo_DIR}/bin/thisk4geo.sh && . ${k4geo_DIR}/bin/thisk4geo.sh


#--------------------------------------------------------------------------------
#     k4edm4hep2lcioconv
#--------------------------------------------------------------------------------
# --- additional k4edm4hep2lcioconv commands ------- 
export PATH="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/k4edm4hep2lcioconv/v00-11/install/bin:$PATH"
export LD_LIBRARY_PATH="/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/k4edm4hep2lcioconv/v00-11/install/lib:/afs/desy.de/project/ilcsoft/sw/x86_64_gcc131_el9/v02-03-04/k4edm4hep2lcioconv/v00-11/install/lib64:$LD_LIBRARY_PATH"


#--------------------------------------------------------------------------------
#     MySQL
#--------------------------------------------------------------------------------
export MYSQL_HOME="/cvmfs/sft.cern.ch/lcg/releases/mysql/10.5.20-7d082/x86_64-el9-gcc13-opt"
export MYSQL="$MYSQL_HOME"
export MYSQL_PATH="$MYSQL_HOME"
export MYSQL_LIBDIR="$MYSQL_HOME/lib64/mysql"
# --- additional MySQL commands ------- 
export PATH="$MYSQL_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$MYSQL_HOME/lib64:$MYSQL_HOME/lib:$MYSQL_HOME/lib64/mysql:$MYSQL_HOME/lib/mysql:$LD_LIBRARY_PATH"

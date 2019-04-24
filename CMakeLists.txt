cmake_minimum_required (VERSION 3.1)
project (MeshClust2)

include_directories(src/exception src/nonltr src/utility src/cluster src/prediction src/clutil src/fastcar)
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})


add_library(Fastcar
  ${CMAKE_SOURCE_DIR}/src/fastcar/FC_Runner.cpp
)

add_library(ClusterUtil
  ${CMAKE_SOURCE_DIR}/src/clutil/DivergencePoint.cpp
  ${CMAKE_SOURCE_DIR}/src/clutil/Histogram.cpp
  ${CMAKE_SOURCE_DIR}/src/clutil/Loader.cpp
  ${CMAKE_SOURCE_DIR}/src/clutil/SingleFileLoader.cpp
  ${CMAKE_SOURCE_DIR}/src/clutil/Progress.cpp
  ${CMAKE_SOURCE_DIR}/src/clutil/Datatype.cpp
  ${CMAKE_SOURCE_DIR}/src/clutil/Clock.cpp
)

add_library(Predict
  ${CMAKE_SOURCE_DIR}/src/predict/Feature.cpp
  ${CMAKE_SOURCE_DIR}/src/predict/GLM.cpp
  ${CMAKE_SOURCE_DIR}/src/predict/HandleSeq.cpp
  ${CMAKE_SOURCE_DIR}/src/predict/Matrix.cpp
  ${CMAKE_SOURCE_DIR}/src/predict/MultiMute.cpp
  ${CMAKE_SOURCE_DIR}/src/predict/Predictor.cpp
  ${CMAKE_SOURCE_DIR}/src/predict/SingMute.cpp
  ${CMAKE_SOURCE_DIR}/src/predict/FeatureSelector.cpp
  ${CMAKE_SOURCE_DIR}/src/predict/GreedySelector.cpp
  ${CMAKE_SOURCE_DIR}/src/predict/BestFirstSelector.cpp
)

add_library(Cluster
  ${CMAKE_SOURCE_DIR}/src/cluster/ClusterFactory.cpp
  ${CMAKE_SOURCE_DIR}/src/cluster/CRunner.cpp
  ${CMAKE_SOURCE_DIR}/src/cluster/Trainer.cpp
  ${CMAKE_SOURCE_DIR}/src/cluster/bvec.cpp
  ${CMAKE_SOURCE_DIR}/src/cluster/bvec_iterator.cpp

)

add_library(Exception
  ${CMAKE_SOURCE_DIR}/src/exception/FileDoesNotExistException.cpp
  ${CMAKE_SOURCE_DIR}/src/exception/InvalidInputException.cpp
  ${CMAKE_SOURCE_DIR}/src/exception/InvalidOperationException.cpp
  ${CMAKE_SOURCE_DIR}/src/exception/InvalidOrderOfOperationsException.cpp
  ${CMAKE_SOURCE_DIR}/src/exception/InvalidScoreException.cpp
  ${CMAKE_SOURCE_DIR}/src/exception/InvalidStateException.cpp
)

add_library(Nonltr
  ${CMAKE_SOURCE_DIR}/src/nonltr/ChromDetectorMaxima.cpp
  ${CMAKE_SOURCE_DIR}/src/nonltr/ChromListMaker.cpp
  ${CMAKE_SOURCE_DIR}/src/nonltr/Chromosome.cpp
  ${CMAKE_SOURCE_DIR}/src/nonltr/ChromosomeOneDigit.cpp
  ${CMAKE_SOURCE_DIR}/src/nonltr/ChromosomeOneDigitDna.cpp
  ${CMAKE_SOURCE_DIR}/src/nonltr/ChromosomeOneDigitProtein.cpp
  ${CMAKE_SOURCE_DIR}/src/nonltr/ChromosomeRandom.cpp
  ${CMAKE_SOURCE_DIR}/src/nonltr/DetectorMaxima.cpp
  ${CMAKE_SOURCE_DIR}/src/nonltr/HMM.cpp
  ${CMAKE_SOURCE_DIR}/src/nonltr/LocationList.cpp
  ${CMAKE_SOURCE_DIR}/src/nonltr/LocationListCollection.cpp
  ${CMAKE_SOURCE_DIR}/src/nonltr/Scanner.cpp
  ${CMAKE_SOURCE_DIR}/src/nonltr/Scorer.cpp
  ${CMAKE_SOURCE_DIR}/src/nonltr/TableBuilder.cpp
  ${CMAKE_SOURCE_DIR}/src/nonltr/Trainer.cpp
)

add_library(Utility
  ${CMAKE_SOURCE_DIR}/src/utility/EmptyLocation.cpp
  ${CMAKE_SOURCE_DIR}/src/utility/GlobAlignE.cpp
  ${CMAKE_SOURCE_DIR}/src/utility/Location.cpp
  ${CMAKE_SOURCE_DIR}/src/utility/Util.cpp
)

target_include_directories(Exception PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(Nonltr PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(Utility PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(Cluster PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(Fastcar PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(ClusterUtil PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(Predict PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set (HEADER_FILES
  ${CMAKE_SOURCE_DIR}/src/nonltr/KmerHashTable.h
  ${CMAKE_SOURCE_DIR}/src/nonltr/EnrichmentMarkovView.h
  ${CMAKE_SOURCE_DIR}/src/nonltr/TableBuilder.h
)

set (CMAKE_CXX_COMPILER g++)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_FLAGS "-fopenmp -g -O3 -march=native")

target_compile_definitions(Cluster PRIVATE VERSION="2.3.0")
target_compile_definitions(Fastcar PRIVATE VERSION="0.7.1")

add_executable(Red ${CMAKE_SOURCE_DIR}/src/nonltr/RepeatsDetector.cpp )
add_executable(meshclust2 ${CMAKE_SOURCE_DIR}/src/cluster/meshclust2.cpp)
add_executable(fastcar ${CMAKE_SOURCE_DIR}/src/fastcar/fastcar.cpp)

target_link_libraries(Red Exception Nonltr Utility ${HEADER_FILES})
target_link_libraries(Utility Exception ${HEADER_FILES})
target_link_libraries(Nonltr Utility Exception ${HEADER_FILES})
target_link_libraries(ClusterUtil Nonltr ${HEADER_FILES})
target_link_libraries(Predict ClusterUtil Nonltr ${HEADER_FILES})
target_link_libraries(meshclust2 Cluster Nonltr ClusterUtil Predict ${HEADER_FILES})
target_link_libraries(fastcar Nonltr ClusterUtil Fastcar Predict ${HEADER_FILES})

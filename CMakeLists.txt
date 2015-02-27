# Copyright: 2015 Fondazione Istituto Italiano di Tecnologia
# Author: Silvio Traversaro
# CopyPolicy: Released under the terms of the GNU GPL v2.0.
#

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11)

SET(PROJECTNAME icubArmRegressor)

PROJECT(${PROJECTNAME})

FIND_PACKAGE(YARP REQUIRED)
FIND_PACKAGE(iDynTree REQUIRED)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${YARP_MODULE_PATH})

# SET(folder_source main.cpp multitaskRecursiveLinearEstimator.cpp multitaskSVDLinearEstimator.cpp)
SET(folder_source main_icubdyn.cpp multitaskRecursiveLinearEstimator.cpp multitaskSVDLinearEstimator.cpp)

SOURCE_GROUP("Source Files" FILES ${folder_source})

message(STATUS "iDynTree_INCLUDE_DIRS is ${iDynTree_INCLUDE_DIRS}")
message(STATUS "iDynTree_LIBRARIES is ${iDynTree_LIBRARIES}")

INCLUDE_DIRECTORIES(${iDynTree_INCLUDE_DIRS}
                    ${YARP_INCLUDE_DIRS}
                    ${PROJECT_SOURCE_DIR}/include
                    )


# ADD_EXECUTABLE(${PROJECTNAME} ${folder_source})
ADD_EXECUTABLE(${PROJECTNAME}_icubdyn ${folder_source})

TARGET_LINK_LIBRARIES(${PROJECTNAME}_icubdyn ${YARP_LIBRARIES}
                                     ${iDynTree_LIBRARIES}
                                     )                                     



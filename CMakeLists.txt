cmake_minimum_required(VERSION 3.14)
project(magnetic_wire_control)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-Wno-deprecated-declarations -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../")


find_package(Eigen3 REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

include_directories(
        ${EIGEN_INCLUDE_DIRS}
        ${OPENGL_INCLUDE_DIRS}
        ${GLUT_INCLUDE_DIRS}
        )


add_executable(optimizeDER
        # environment files
        src_optimization/main.cpp
        src_optimization/world.cpp

        # mechanics
        src_optimization/elasticRod.cpp
        src_optimization/elasticStretchingForce.cpp
        src_optimization/elasticBendingForce.cpp
        src_optimization/elasticTwistingForce.cpp
        src_optimization/externalMagneticForce.cpp

        # inertia
        src_optimization/inertialForce.cpp
        
        # external forces
        src_optimization/dampingForce.cpp
        src_optimization/externalGravityForce.cpp
        src_optimization/externalContactForce.cpp

        # solver
        src_optimization/timeStepper.cpp

        # utils
        src_optimization/setInput.cpp
        )


add_executable(simDER
        # environment files
        src_simulation/main.cpp
        src_simulation/world.cpp

        # mechanics
        src_simulation/elasticRod.cpp
        src_simulation/elasticStretchingForce.cpp
        src_simulation/elasticBendingForce.cpp
        src_simulation/elasticTwistingForce.cpp
        src_simulation/externalMagneticForce.cpp

        # inertia
        src_simulation/inertialForce.cpp
        
        # external forces
        src_simulation/dampingForce.cpp
        src_simulation/externalGravityForce.cpp
        src_simulation/externalContactForce.cpp

        # solver
        src_simulation/timeStepper.cpp

        # utils
        src_simulation/setInput.cpp
        )



target_include_directories(optimizeDER PUBLIC
        ${CMAKE_SOURCE_DIR}/src_optimization
)


target_include_directories(simDER PUBLIC
        ${CMAKE_SOURCE_DIR}/src_simulation
)

target_link_libraries(optimizeDER
        ${OPENGL_LIBRARIES}
        ${GLUT_LIBRARIES}
        Eigen3::Eigen
        lapack
        gfortran
        pthread
        rt
        m)

target_link_libraries(simDER
        ${OPENGL_LIBRARIES}
        ${GLUT_LIBRARIES}
        Eigen3::Eigen
        lapack
        gfortran
        pthread
        rt
        m)
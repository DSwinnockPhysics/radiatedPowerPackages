include_directories(/home/phrnxbhome/Code/C++/Libraries/BFieldModels/include)
include_directories(/home/phrnxbhome/Code/C++/Libraries/Eigen3.4.0/eigen-3.4.0)
add_library(ElectronDynamics BaseField.cxx QTNMFields.cxx BorisSolver.cxx TrajectoryGen.cxx MyCustomBorisSolver.cxx QTNMApproxFields.cxx)
target_link_libraries(ElectronDynamics PUBLIC /home/phrnxbhome/Code/C++/Libraries/BFieldModels/MyBFieldModels/build/libBFieldModels.a BasicFunctions ${ROOT_LIBRARIES} ${Boost_MATH_LIBRARY})

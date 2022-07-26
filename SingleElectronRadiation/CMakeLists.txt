include_directories(/home/phrnxbhome/Code/C++/Libraries/Eigen3.4.0/eigen-3.4.0)
add_library(SingleElectronRadiation SingleElectronRadiation.cxx)
target_link_libraries(SingleElectronRadiation PUBLIC ${ROOT_LIBRARIES} ${ROOTFFTW_LIB})

add_library(libhl SHARED FastMath.cpp FastMath.h FastLibpano13.h FastLibpano13.cpp)
target_link_libraries(libhl ${PANO_LIBRARIES})
install(TARGETS libhl LIBRARY DESTINATION ${LIBDIR}/hugin NAMELINK_SKIP)

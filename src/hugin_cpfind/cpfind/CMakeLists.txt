add_executable(cpfind PanoDetector.cpp PanoDetectorLogic.cpp TestCode.cpp Utils.cpp main.cpp ImageImport.h
                         KDTree.h KDTreeImpl.h PanoDetector.h PanoDetectorDefs.h TestCode.h Tracer.h Utils.h
)
set(cpfind_deps localfeatures ${ZTHREAD_LIBRARIES} ${image_libs} ${common_libs} celeste)

IF(NOT WIN32)
    set(cpfind_deps ${cpfind_deps} pthread)
    IF(FLANN_FOUND)
        set(cpfind_deps ${cpfind_deps} ${FLANN_LIBRARIES})
    ENDIF(FLANN_FOUND)
    IF(ENABLE_HUGIN_LITE)
        set(cpfind_deps ${cpfind_deps} libhl)
    ENDIF(ENABLE_HUGIN_LITE)
ENDIF(NOT WIN32)
target_link_libraries(cpfind ${cpfind_deps})

install(TARGETS cpfind DESTINATION ${BINDIR})


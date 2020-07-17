# In:
#  OPENVDB_ROOT
#
# Out:
#  OPENVDB_FOUND
#  OPENVDB_INCLUDE_DIRS
#  OPENVDB_LIBRARY_DIRS
#  OPENVDB_DEFINITIONS

if(NOT OPENVDB_ROOT)
	if (IS_DIRECTORY "/usr/local/opt/openvdb/") #  mac/homebrew (version independent))
		set(OPENVDB_ROOT "/usr/local/opt/openvdb/")
	elseif(IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/openVDB")
		set(OPENVDB_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/openVDB")
	endif()
else()
	set(OPENVDB_ROOT "" CACHE PATH "OpenVDB base path")
endif()

if(OPENVDB_ROOT)
	set(OPENVDB_INCLUDE_DIR "${OPENVDB_ROOT}/include" CACHE PATH "Include files for openVDB")
	set(OPENVDB_LIBRARY_DIR "${OPENVDB_ROOT}/lib" CACHE PATH "OpenVDB libraries")

	set(OPENVDB_DEFINITIONS "")
	list(APPEND OPENVDB_DEFINITIONS -DOPENVDB_DLL)

	if(APPLE)
		set(OPENVDB_LIBRARY "openvdb")
	else()
		set(OPENVDB_LIBRARY "${OPENVDB_LIBRARY_DIR}/openVDB.lib" CACHE FILEPATH "OpenVDB library filepath")
		list(APPEND OPENVDB_DEFINITIONS -DOPENVDB_3_ABI_COMPATIBLE)
	endif()

	set(OPENVDB_FOUND TRUE)
endif()
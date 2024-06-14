if( HAVE_CPU)
  add_subdirectory( cpu )
endif()
if( HAVE_GPU )
  add_subdirectory( gpu )
endif()

# install headers
file( GLOB etrans_interface include/biper/* include/etrans/*)
install(
  FILES        ${etrans_interface}
  DESTINATION  include/ectrans
)

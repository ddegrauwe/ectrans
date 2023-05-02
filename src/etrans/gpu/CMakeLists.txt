# Create library etrans_${prec}.so



if(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")

  # Compile setup_trans with pinned memory to improve data movement performance. 
  ectrans_add_compile_options(
      SOURCES external/setup_trans.F90
      #FLAGS   "-gpu=pinned,deepcopy,fastmath,nordc")
      FLAGS   "-gpu=pinned,fastmath")
  # TODO: check if it is sufficient to only set "-gpu=pinned" which appends rather than overwrites

endif()

# list of source files
ecbuild_list_add_pattern( LIST etrans_src
                          GLOB
                                external/*.F90
								internal/*.F90
								aux/*.F90
                          QUIET
                        )#

#set_source_files_properties( internal/ftinv_mod.F90 PROPERTIES COMPILE_OPTIONS "-O2" )
#set_source_files_properties( internal/ftdir_mod.F90 PROPERTIES COMPILE_OPTIONS "-O2" )


foreach( prec sp dp )
  if( HAVE_${prec} )
    foreach( gpumethod acc)
      if( HAVE_${gpumethod})

        # We build an object library first. And then use these objects to create a shared and static library.

        ectrans_add_library(
          TARGET           etrans_gpu_object_${prec}_${gpumethod}
          TYPE             OBJECT
          SOURCES          ${etrans_src}
          PUBLIC_INCLUDES 
				$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include/etrans>
				$<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
				$<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}/trans_gpu_${prec}_${gpumethod}>
				$<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}/etrans_gpu_${prec}_${gpumethod}>
				$<INSTALL_INTERFACE:include>
				$<INSTALL_INTERFACE:include/ectrans>

          PUBLIC_LIBS
						  fiat
						  parkind_${prec}
						  
          PRIVATE_LIBS    hip::hipfft
                          roc::hipblas
			                    roc::rocblas
			                    roc::rocfft
                          MPI::MPI_Fortran
                          mpifort
                          omptarget
                          #MPI::MPI_C
             			  MPI::MPI_CXX
                          #nvhpcwrapnvtx
                          #${LAPACK_LIBRARIES}
        )

        target_link_options   ( etrans_gpu_object_${prec}_${gpumethod} INTERFACE "$<$<LINK_LANG_AND_ID:Fortran,NVHPC>:${${gpumethod}_link_flags}>" )
        target_compile_options( etrans_gpu_object_${prec}_${gpumethod} PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${${gpumethod}_flags}>)

        ectrans_target_fortran_module_directory(
          TARGET            etrans_gpu_object_${prec}_${gpumethod}
          #MODULE_DIRECTORY  ${PROJECT_BINARY_DIR}/module/etrans_gpu_${prec}_${gpumethod}
		  MODULE_DIRECTORY  ${CMAKE_BINARY_DIR}/module/etrans_gpu_${prec}_${gpumethod}
          INSTALL_DIRECTORY module/etrans_gpu_${prec}_${gpumethod}
        )

        if( HAVE_OMP )
          target_link_libraries( etrans_gpu_object_${prec}_${gpumethod} PRIVATE OpenMP::OpenMP_Fortran )
        endif()

        if( prec STREQUAL sp )
          target_compile_definitions( etrans_gpu_object_${prec}_${gpumethod} PRIVATE TRANS_SINGLE PARKINDTRANS_SINGLE )
        endif()

        if( gpumethod STREQUAL acc )
          target_compile_definitions( etrans_gpu_object_${prec}_${gpumethod} PRIVATE ACCGPU )
          target_link_libraries( etrans_gpu_object_${prec}_${gpumethod} PRIVATE OpenACC::OpenACC_Fortran )
        endif()

        if( gpumethod STREQUAL omp )
          target_compile_definitions( etrans_gpu_object_${prec}_${gpumethod} PRIVATE OMPGPU )
        endif()

        if( HAVE_GPU_AWARE_MPI )
          target_compile_definitions( etrans_gpu_object_${prec}_${gpumethod} PRIVATE USE_CUDA_AWARE_MPI_FT )
        endif()

        ectrans_add_library(
          TARGET           etrans_gpu_shared_${prec}_${gpumethod}
          OUTPUT_NAME      etrans_gpu_${prec}_${gpumethod}
          TYPE             SHARED
          LINKER_LANGUAGE  Fortran
          PUBLIC_LIBS      etrans_gpu_object_${prec}_${gpumethod}
        )

        ectrans_add_library(
          TARGET           etrans_gpu_${prec}_${gpumethod}
          TYPE             STATIC
          LINKER_LANGUAGE  Fortran
          PUBLIC_LIBS      etrans_gpu_object_${prec}_${gpumethod}
        )
      endif()
    endforeach()
  endif()
endforeach()

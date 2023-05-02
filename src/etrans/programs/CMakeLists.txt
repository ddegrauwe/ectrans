foreach( prec sp dp )

if ( HAVE_CPU )

	  if( HAVE_${prec} )

		ecbuild_add_executable(TARGET  aatestprog_${prec}
							   SOURCES aatestprog.F90
							   LIBS
								 etrans_${prec}
								 trans_${prec}
								 fiat
								 parkind_${prec}
							   LINKER_LANGUAGE Fortran
							  )	
		target_link_libraries( aatestprog_${prec} ${OpenMP_Fortran_FLAGS} OpenMP::OpenMP_Fortran )

	  endif()
		
endif() # HAVE_CPU

if ( HAVE_GPU )
	if ( HAVE_${prec} )
      foreach( gpumethod acc )
        if( HAVE_${gpumethod} )
          if( gpumethod STREQUAL acc )
            ecbuild_add_executable(TARGET  aatestprog-gpu-${prec}-${gpumethod}
                                  SOURCES aatestprog.F90
				                          LINKER_LANGUAGE Fortran
                                  LIBS
                                    etrans_gpu_${prec}_${gpumethod}
                                    trans_gpu_${prec}_${gpumethod}
                                    fiat
                                    parkind_${prec}
                                    OpenACC::OpenACC_Fortran
                                  )
            target_compile_definitions( aatestprog-gpu-${prec}-${gpumethod} PRIVATE ACCGPU )
          endif()

          if( gpumethod STREQUAL omp )

            ecbuild_add_executable(TARGET  aatestprog-gpu-${prec}-${gpumethod}
                                  SOURCES aatestprog.F90
				                          LINKER_LANGUAGE Fortran
                                  LIBS
                                    etrans_gpu_${prec}_${gpumethod}
                                    trans_gpu_${prec}_${gpumethod}
                                    fiat
                                    parkind_${prec}
                                  )
            target_compile_definitions( aatestprog-gpu-${prec}-${gpumethod} PRIVATE OMPGPU )
          endif()
                    
          target_link_options   ( aatestprog-gpu-${prec}-${gpumethod} INTERFACE "$<$<LINK_LANG_AND_ID:Fortran,NVHPC>:${${gpumethod}_link_flags}>" )
          target_compile_options( aatestprog-gpu-${prec}-${gpumethod} PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${${gpumethod}_flags}>)

        endif()
      endforeach()
	endif() # HAVE_${prec}
endif() # HAVE_GPU

endforeach() # prec 

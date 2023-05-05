# Limited-Area Model spectral transforms on LUMI

## Usage (aka the most interesting part)

### (re)Compilation

    source ~dadegrau2/ENV_lumi
    rm -rf ${BUILDDIR}/ectrans ${INSTALLDIR}/ectrans
    mkdir -p ${BUILDDIR}/ectrans
    cd ${BUILDDIR}/ectrans
    ecbuild --prefix=${INSTALLDIR}/ectrans -Dfiat_ROOT=${INSTALLDIR}/fiat -DBUILD_SHARED_LIBS=OFF -DENABLE_FFTW=OFF -DENABLE_GPU=ON -DENABLE_OMPGPU=OFF -DENABLE_ACCGPU=ON -DENABLE_TESTS=OFF -DENABLE_GPU_AWARE_MPI=ON -DENABLE_CPU=ON -DENABLE_ETRANS=ON ${SOURCEDIR}/ectrans
    make -j32
	make -j32   # note: dependencies aren't worked out entirely correctly by cmake, therefore a second make is necessary.
    make install

### Test run on single GPU

Allocate GPU resource with

    salloc --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --gpus-per-node=1 --account=project_462000140 --partition=standard-g --time=04:00:00 --mem=0

Recompile/run with 

    cd ${BASEDIR}/test/
    make -j16 -C ${BUILDDIR}/ectrans/ install

	args="--truncation 79 --nproma 32 --vordiv --scders --uvders --nfld 1 --nlev 10 --norms --check 10"
    srun ${INSTALLDIR}/ectrans/bin/ectrans-benchmark-sp  ${args}             # run global benchmark on CPU
    srun ${INSTALLDIR}/ectrans/bin/ectrans-benchmark-gpu-sp-acc ${args}      # run global benchmark on GPU
	
	args="--nlon 128 --nlat 128 --nproma 32 --vordiv --scders --uvders --nfld 5 --nlev 10 --norms --check 10"
	srun ${INSTALLDIR}/ectrans/bin/ectrans-lam-benchmark-sp ${args}          # run LAM benchmark on CPU
	srun ${INSTALLDIR}/ectrans/bin/ectrans-lam-benchmark-gpu-sp-acc ${args}  # run LAM benchmark on GPU

Note: recompiling like this may not be sufficient when modifying e.g. hip.cc files; do a full recompilation (above) in this case.

Note: a bug still resides in the global gpu runs when setting nfld>1.

### Test run on multiple GPUs

Allocate GPU resources with

    salloc --nodes=1 --ntasks-per-node=8 --gpus-per-node=8 --account=project_462000140 --partition=standard-g --time=04:00:00 --mem=0

Make sure to set

    export MPICH_GPU_SUPPORT_ENABLED=1
    CPU_BIND="map_cpu:48,56,16,24,1,8,32,40"

Then launch global cpu/gpu runs with
	
    args="--truncation 79 --nproma 32 --vordiv --scders --uvders --nfld 1 --nlev 10 --norms --check 10"
    srun --cpu-bind=${CPU_BIND} ./select_gpu ${INSTALLDIR}/ectrans/bin/ectrans-benchmark-sp  ${args}
    srun --cpu-bind=${CPU_BIND} ./select_gpu ${INSTALLDIR}/ectrans/bin/ectrans-benchmark-gpu-sp-acc ${args}

while the LAM cases are launched with

    args="--nlon 256 --nlat 256 --nproma 32 --niter 1 --vordiv --scders --uvders --nfld 5 --nlev 80 --norms --check 10 --dump-values"
    srun --cpu-bind=${CPU_BIND} ./select_gpu ${INSTALLDIR}/ectrans/bin/ectrans-lam-benchmark-sp ${args}
    srun --cpu-bind=${CPU_BIND} ./select_gpu ${INSTALLDIR}/ectrans/bin/ectrans-lam-benchmark-gpu-sp-acc ${args}

The `select_gpu` wrapper can be found on the lumi documentation pages.

Small note: the `nfld>1` case gives a crash in the global run, both on cpu and on gpu.

## Prerequisites: ecbuild and fiat

### ecbuild installation

    cd ${SOURCEDIR}
    git clone https://github.com/ecmwf/ecbuild.git
    cd ecbuild
    git checkout master
    sed -i -e "s/-Gfast//" cmake/compiler_flags/Cray_Fortran.cmake # remove obsolete switch -Gfast
    mkdir -p ${BUILDDIR}/ecbuild
    cd ${BUILDDIR}/ecbuild
    ${SOURCEDIR}/ecbuild/bin/ecbuild --prefix=${INSTALLDIR}/ecbuild ${SOURCEDIR}/ecbuild
    make
    make install

### fiat installation
With Cray compiler on lumi, one gets into trouble with OpenMP: for some reason, during linking the openmp library isn't found... This is solved by adding `${OpenMP_C_FLAGS}` in `programs/CMakeLists.txt`: `target_link_libraries( fiat-printbinding ${OpenMP_C_FLAGS} OpenMP::OpenMP_C )`

    cd ${SOURCEDIR}
    git clone https://github.com/ecmwf-ifs/fiat
    cd fiat
    git checkout main
	# ADD OpenMP::OpenMP_C here!
    rm -rf ${BUILDDIR}/fiat
    mkdir -p ${BUILDDIR}/fiat
    cd ${BUILDDIR}/fiat
    ecbuild -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=${INSTALLDIR}/fiat -DENABLE_MPI=ON -DBUILD_SHARED_LIBS=OFF -DENABLE_TESTS=OFF ${SOURCEDIR}/fiat
    make -j16
    make install

## Code organization and data layout

Just for reference: a somewhat simplified overview of routines and data layout

### Original CPU code


Inverse transforms:
```
    einv_trans                 # input: PSPVOR(NLEV,NSPEC2), PSPDIV(NLEV,NSPEC2), PSPSC3A(NFLD*NLEV,NSPEC2)
      einv_trans_ctl
        eltinv_ctl
          DO JM=1,NX_l       # loop over x-wavenumbers
            eltinv             # north-south transforms
              eprfi1b          # transpose to (NY,NFLD)
              eleinv
                plan_fft
                execute_fft    # small-batched transform along y with stride=1, distance=ny+2, lot=nfld
              easre1b          # transpose to FOUBUF_IN(NFLD,NY)
          ENDDO
          trmtol               # inter-GPU communications to FOUBUF(NFLD,NY)
          
        eftinv_ctl
          DO JGL=1,NY_l   # loop over latitudes
            fourier_in         # copy FOUBUF to ZGTF(NFLD,NX,JGL)
            eftinv
              plan_fft
              execute_fft      # small-batched transform along x with stride=nfld,distance=1, lot=NFLD
          ENDDO
          trltog               # inter-GPU comms and transpose to PGP3A(NPROMA,NLEV*NFLD,NBLK),PGPUV(NPROMA,NLEV*2,NBLK)
```

Direct transforms:
```
    edir_trans_ctl            # input PGP3A(NPROMA,NFLD*NLEV,NBLK)
      eftdir_ctl
        trgtol                # inter-GPU comms and transposition to ZGTF(NFLD,NX)
        DO JGL=1,NY_l    # loop over latitudes
          ftdir
            plan_fft
            execute_fft       # batched transform along x with stride=NFLD, distance=1, lot=NFLD
          fourier_out         # copy to FOUBUF_IN(NFLD,NX)
        ENDDO
      eltdir_ctl
        trltom                # inter-GPU comms to FOUBUF(NFLD,NY)
        DO JM=1,NX_l        # loop over x-wavenumbers
          eltdir
            eprfi2b           # transpose to PFFT(NY,NFLD)
            eledir
              plan_fft
              execute_fft     # batched transform along y with stride=1, distance=ny+2, lot=NFLD
        eupdsp
          eupdspb             # transpose to PSPSC3A(NFLD,NSPEC2)
```

Two aspects of this code and data organization are quite striking: (i) 2x3 transpositions are performed. On CPU this doesn't seem to be very expensive; (ii) the Fourier transforms are performed inside loops over JGL and JM. This means that the batch size of the transforms is quite limited.

### GPU code

To unleash the computing power of GPUs, as much as possible parallellism should be exposed. In the case of the spectral transforms, this means the batch size should be as large as possible. Instead of using the loops over JGL and JM and small batch sizes as for the GPU code, the following organization is taken, which removes the loops and increases the batch sizes.

Inverse transforms:
```
    einv_trans               # input: PSPVOR(NLEV,NSPEC2), PSPDIV(NLEV,NSPEC2), PSPSC3A(NFLD*NLEV,NSPEC2)
      einv_trans_ctl
        eltinv_ctl
          eltinv             # north-south transforms
            eprfi1b          # transpose to (NY,NFLD)
            eleinv
              plan_fft
              execute_fft    # batched transform along y with stride=1, distance=ny+2
            easre1b          # transpose to FOUBUF_IN(NFLD,NY)
          trmtol             # inter-GPU communications to FOUBUF(NFLD,NY)
          
        eftinv_ctl
          efourier_in        # transpose FOUBUF to ZGTF(NX,NFLD)
          eftinv
            plan_fft
            execute_fft      # batched transform along x with stride=1,distance=nx+2
          trltog             # inter-GPU comms to PGP3A(NPROMA,NLEV*NFLD,NBLK),PGPUV(NPROMA,NLEV*2,NBLK)
```

Direct transforms:
```
    edir_trans_ctl           # input PGP3A(NPROMA,NFLD*NLEV,NBLK)
      eftdir_ctl
        trgtol               # inter-GPU comms to ZGTF(NX,NY_l*NFLD)
        eftdir
          plan_fft
          execute_fft        # batched transform along x with stride=1, distance=nx+2, lot=NY_l*NFLD
        efourier_out         # transpose to FOUBUF_IN(NFLD,NY_l,NX)
      eltdir_ctl
        trltom               # inter-GPU comms to FOUBUF(NFLD*NX_l,NY)
        eltdir
          eprfi2b            # transpose to PFFT(NY,NX_l*NFLD)
          eledir
            plan_fft
            execute_fft      # batched transform along y with stride=1, distance=ny+2, lot=NX_l*NFLD
        eupdsp
          eupdspb            # transpose to PSPSC3A(NFLD,NSPEC2)
```

Although the number of transpositions remains the same (2x3), they are performed at a different place. This is necessary because a batched Fourier transform must be performed either on the first dimension, or on the last dimension. As a consequence, the input to `trltog` and the output from `trgtol` is transposed w.r.t. the CPU version and a nasty switch `LDTRANSPOSED` is necessary in those routines.


## Optimizations

Just documenting some ideas...

### Tiling transpositions

The most expensive routines (especially on GPU) are those where data are transposed, i.e. where the leading dimension changes. The reason is that a transposition is very much nonlocal in terms of data access: coalesced data access is only possible for the input or for the output, but not for both.

One possibility to solve this is to use a tiled approach, either coded explicitly, as explained [here](https://developer.nvidia.com/blog/efficient-matrix-transpose-cuda-cc/), or by using the OpenACC `tile` directive. How well this directive is treated by the compiler (as compared to using shared memory) is to be tested. Also the optimal tile size in the `tile` directive is machine-dependent and could be optimized.

### Avoiding transpositions

Considering the direct transforms, at least one transposition is necessary, because the leading dimension of the input is `NPROMA` ($x$), while the leading dimension of the output is `NFLD`. But the other two transpositions (in `eprfi2b` and `eupdspb`) may be avoided. As a matter of fact, in the nvidia-optimized branch by Lukas, they are avoided.

The problem is that these transpositions are necessary to maximize the batch size of the Fourier transforms. A different possibility is to keep the loop over `JM` as in the CPU code, but to use different GPU streams to treat the different x-wavenumbers. According to Lukas, the performance of a batched FFT is better than that of multiplexed FFTs in different streams, but if this avoids two transpositions, it may well be worth the effort.

This can be tested in a small toy application, where the performance of taking FFTs over different dimensions is considered:

* batched along leading dimension:
```
Z(NY,NX_l,NFLD)                   # input/output data
plan_fft           
execute_fft(Z)                    # stride=1, distance=NY, lot=NX_l*NFLD
```
* batched along trailing dimension:
```
Z(NX_l,NFLD,NY)                   # input/output data
plan_fft           
execute_fft(Z)                    # stride=NX_l*NFLD, distance=1, lot=NX_l*NFLD
```
* streamed along trailing dimension
```
Z(NX_l,NY,NFLD)                   # input/output data
DO JFLD=1,NFLD
  hipfft_setstream(JFLD)          # ... or something like that
  plan_fft           
  execute_fft(Z(1,1,JFLD))        # stride=NX_l, distance=1, lot=NX_l
ENDDO
```
* streamed along leading dimension
```
Z(NX_l,NY,NFLD)                   # input/output data
DO JX=1,NX_l
  hipfft_setstream(JX)          # ... or something like that
  plan_fft           
  execute_fft(Z(JX,1,1))        # stride=NX_l, distance=NX_l*NY, lot=NFLD
ENDDO
```

Realistic dimension values are something like `NX_l=128`, `NFLD=240`, `NY=1024`.

Besides avoiding transpositions, if the multiplexed transforms perform well, it would also allow to remove the nasty `LDTRANSPOSED` switch in `trltog` and `trgtol`.

### Revise triple loop in trgtol/trltog

In `trltog` and `trgtol`, data layout is changed from `(NPROMA,NFLD,NBLK)` to `(NX,NFLD)`. It may be better to loop over `NX` and calculate the indices `JLON` and `JBLK`, instead of the current approach which uses a double loop over `JLON` and `JBLK`, calculating `JX` on the go. This is also how it's done in the nvidia-opt branch.

### Overlapping host-device transfers and communications

Considering the direct transforms, in the current code, copying of output data from the device to the host only starts after all calculations are finished: `edir_trans`, there is

```
!$ACC data copyin (PGP  ) if (present (PGP  ))
!$ACC data copyin (PGPUV) if (present (PGPUV))
!$ACC data copyin (PGP3A) if (present (PGP3A))
!$ACC data copyin (PGP3B) if (present (PGP3B))
!$ACC data copyin (PGP2 ) if (present (PGP2 ))
!$ACC data copyout (PSPVOR   ) if (present (PSPVOR   ))
!$ACC data copyout (PSPDIV   ) if (present (PSPDIV   ))
!$ACC data copyout (PSPSCALAR) if (present (PSPSCALAR))
!$ACC data copyout (PSPSC3A  ) if (present (PSPSC3A  ))
!$ACC data copyout (PSPSC3B  ) if (present (PSPSC3B  ))
!$ACC data copyout (PSPSC2   ) if (present (PSPSC2   ))
CALL EDIR_TRANS_CTL(IF_UV_G,IF_SCALARS_G,IF_GP,IF_FS,IF_UV,IF_SCALARS,&
 & PSPVOR,PSPDIV,PSPSCALAR,KVSETUV,KVSETSC,PGP,&
 & PSPSC3A,PSPSC3B,PSPSC2,KVSETSC3A,KVSETSC3B,KVSETSC2,PGPUV,PGP3A,PGP3B,PGP2,&
 & PMEANU,PMEANV,AUX_PROC)
!$ACC end data
!$ACC end data
!$ACC end data
!$ACC end data
!$ACC end data
!$ACC end data
!$ACC end data
!$ACC end data
!$ACC end data
!$ACC end data
!$ACC end data
```

However, looking inside the code deeper, it becomes apparent that calculations on `PSPVOR` have already finished well before the end of `edir_trans_ctl`. So the transfer of `PSPVOR `from the device to the host could already start earlier, while processing of the other fields `PSPDIV`, `PSPSC3A`, `PSPSC2` still continues. In terms of code, this would require changing the `copyout` clauses in `edir_trans` to `create` clauses, and to put an `update host async(1)` clause right after the processing of `PSPVOR` in `eupdsp` (and same for the other fields).

Completely symmetric actions could be taken for the input data of the inverse transforms in `eltinv`.

## History

### Different existing repositories contain elements for this:
* IAL/ectrans_withlam:withlam
integration of CPU lam sources in ectrans
* anmrde/ectrans:gpu_omp
OpenACC/(OpenMP)/hip branch of global with Cray compiler
* ddegrauwe/etrans:gpu_daand_lam
OpenACC/cufft branch for Nvidia compiler
* ddegrauwe/ectrans:gpu_daand_lam
modifs to global Nvidia-targeted code for lam:
    * LDGW switch to work on transposed arrays e.g. PGLAT in TRGTOL
    * generalizing FFTC plans to have stride
    * switch LALLOPERM2 (only used in etrans, not in ectrans)
* ddegrauwe/ectrans:gpu_lumi
OpenACC/rocfft branch for Cray compiler

### Merge plan:
1. Start from anmrde/ectrans:gpu_omp
2. integrate LAM (CPU) by *merging* in IAL/ectrans_withlam:withlam
3. move to GPU by *merging* ddegrauwe/ectrans:gpu_daand_lam
4. introduce optimizations by Lukas M.

### 1. Andreas' branch

When enabling GPU-aware MPI communications, Crya/Lumi complains about quite a lot of scalars not being present in OpenACC regions. Fixes for this were committed here.


### 2. Integrate LAM sources (CPU)

I created a new branch gpu-omp-daand-lam for this. `etrans` sources were taken from git@github.com:ACCORD-NWP/ectrans_withlam.git

Some incompatibilities due to version differences were solved as follows:
* NSTACK_MEMORY_TR isn't present in TPM_GEN (due to ectrans not being the latest IAL version). This was removed from eftdir_ctl, eftinv_ctl, eftdir_ctlad, eftinv_ctlad
* egath_spec was rewritten; not compatible with gath_spec_control that's in ectrans; I put back the version of egath_spec from cy43t2.
* same for edist_grid

### 3. LAM-GPU changes
Put the cpu-specific sources in `cpu`, and put OpenACC/cuFFT sources from Thomas B. in `gpu`.

Changes done for hipfft:
* remove all references to cuda, including tpm_fftc, fftc (cuda fft data type)
* removed LALLOPERM2 (assuming .FALSE.)
* removed LDGW, transposing ZGTF everywhere.
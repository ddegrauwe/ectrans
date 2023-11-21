#include "hip/hip_runtime.h"
#include "hipfft.h"
#ifdef TRANS_SINGLE
typedef hipfftReal HIP_DATA_TYPE_REAL;
#else
typedef hipfftDoubleReal HIP_DATA_TYPE_REAL;
#endif

#define hipBlock_DIM 4
typedef HIP_DATA_TYPE_REAL Real;

__global__ void hip_easre1b_hiph_kernel(Real* pfft, Real* foubuf_in, int* d_npntgtb1, int k_npntgtb1, int r_ndgl, int d_num, int kjgl, int kjgm, int kfield_2)
{
    unsigned int jgl  = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
    unsigned int jm   = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
    unsigned int jfld = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;
    const int offset_istan = 0;

    if (jgl < r_ndgl && jm < d_num && jfld < kfield_2)
    {
        int offset_in       = jgl + kjgl * jm + kjgl * kjgm * jfld;
        int offset_npntgtb1 = jm + jgl * k_npntgtb1; 
        int offset_out      = jfld + (offset_istan+d_npntgtb1[offset_npntgtb1])*kfield_2;

        foubuf_in[offset_out] = pfft[offset_in];
    }   

}

__global__ void hip_easre1b_hiph_lds_kernel(Real* pfft, Real* foubuf_in, int* d_npntgtb1, int k_npntgtb1, int r_ndgl, int d_num, int kjgl, int kjgm, int kfield_2)
{
     __shared__ float hipBlock[hipBlock_DIM][hipBlock_DIM+1][hipBlock_DIM+1];

    unsigned int jgl  = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
    unsigned int jm   = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;
    unsigned int jfld = hipBlockDim_z * hipBlockIdx_z + hipThreadIdx_z;
    const int offset_istan = 0;

    if (jgl < r_ndgl && jm < d_num && jfld < kfield_2)
    {
        int offset_in       = jgl + kjgl * jm + kjgl * kjgm * jfld;

        hipBlock[hipThreadIdx_x][hipThreadIdx_y][hipThreadIdx_z] = pfft[offset_in];

    }   
 
    __syncthreads();
  
    if (jgl < r_ndgl && jm < d_num && jfld < kfield_2)
    {
        int offset_npntgtb1 = jm + jgl * k_npntgtb1; 
        int offset_out      = jfld + (offset_istan+d_npntgtb1[offset_npntgtb1])*kfield_2;

        foubuf_in[offset_out] = hipBlock[hipThreadIdx_x][hipThreadIdx_y][hipThreadIdx_z] ;
    }   

}

extern "C"
void easre1b_hiph_(Real* pfft, Real* foubuf_in, int* d_npntgtb1, int k_npntgtb1, int r_ndgl, int d_nump, int kjgl, int kjgm, int kfield, bool lds)
{
    int kfield_2 = 2 * kfield;
    dim3 grid(r_ndgl / hipBlock_DIM, d_nump / hipBlock_DIM, kfield_2 / hipBlock_DIM);
    dim3 threads(hipBlock_DIM, hipBlock_DIM, hipBlock_DIM);

    if (lds)
    {
        hipLaunchKernelGGL((hip_easre1b_hiph_lds_kernel), grid,threads, 0, 0, pfft, foubuf_in, d_npntgtb1, k_npntgtb1, r_ndgl, d_nump, kjgl, kjgm, kfield_2);
    }
    else
    {
        hipLaunchKernelGGL((hip_easre1b_hiph_kernel), grid,threads, 0, 0, pfft, foubuf_in, d_npntgtb1, k_npntgtb1, r_ndgl, d_nump, kjgl, kjgm, kfield_2);
    }


    if (hipDeviceSynchronize() != hipSuccess)
    {
        fprintf(stderr, "Hip error: Failed to synchronize\n");
    }
    hipError_t err0 = hipPeekAtLastError ();
    if (err0 != hipSuccess)
    {
        fprintf(stderr, "Hip error\n");
    }
}

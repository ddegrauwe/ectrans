#define hipfftSafeCall(err) __hipfftSafeCall(err, __FILE__, __LINE__)
#include "hip/hip_runtime.h"
#include "hipfft.h"
#include "stdio.h"
#include "execute_plan_ffth.hip.h"

#ifdef TRANS_SINGLE
typedef hipfftComplex HIP_DATA_TYPE_COMPLEX;
typedef hipfftReal HIP_DATA_TYPE_REAL;
#else
typedef hipfftDoubleComplex HIP_DATA_TYPE_COMPLEX;
typedef hipfftDoubleReal HIP_DATA_TYPE_REAL;
#endif


static const char *_hipGetErrorEnum(hipfftResult error)
{
    switch (error)
    {
        case HIPFFT_SUCCESS:
            return "HIPFFT_SUCCESS";

        case HIPFFT_INVALID_PLAN:
            return "HIPFFT_INVALID_PLAN";

        case HIPFFT_ALLOC_FAILED:
            return "HIPFFT_ALLOC_FAILED";

        case HIPFFT_INVALID_TYPE:
            return "HIPFFT_INVALID_TYPE";

        case HIPFFT_INVALID_VALUE:
            return "HIPFFT_INVALID_VALUE";

        case HIPFFT_INTERNAL_ERROR:
            return "HIPFFT_INTERNAL_ERROR";

        case HIPFFT_EXEC_FAILED:
            return "HIPFFT_EXEC_FAILED";

        case HIPFFT_SETUP_FAILED:
            return "HIPFFT_SETUP_FAILED";

        case HIPFFT_INVALID_SIZE:
            return "HIPFFT_INVALID_SIZE";

        case HIPFFT_UNALIGNED_DATA:
            return "HIPFFT_UNALIGNED_DATA";

        case HIPFFT_INCOMPLETE_PARAMETER_LIST:
            return "HIPFFT_INCOMPLETE_PARAMETER_LIST";

        case HIPFFT_INVALID_DEVICE:
            return "HIPFFT_INVALID_DEVICE";

        case HIPFFT_PARSE_ERROR:
            return "HIPFFT_PARSE_ERROR";

        case HIPFFT_NO_WORKSPACE:
            return "HIPFFT_NO_WORKSPACE";

        case HIPFFT_NOT_IMPLEMENTED:
            return "HIPFFT_NOT_IMPLEMENTED";

        case HIPFFT_NOT_SUPPORTED:
            return "HIPFFT_NOT_SUPPORTED";
    }

    return "<unknown>";
}

inline void __hipfftSafeCall(hipfftResult err, const char *file, const int line)
{
    if( HIPFFT_SUCCESS != err) {
        fprintf(stderr, "HIPFFT error at 1\n");
        fprintf(stderr, "HIPFFT error in file '%s'\n",__FILE__);
        fprintf(stderr, "HIPFFT error at 2\n");
        /*fprintf(stderr, "HIPFFT error line '%s'\n",__LINE__);*/
        fprintf(stderr, "HIPFFT error at 3\n");
        /*fprintf(stderr, "HIPFFT error in file '%s', line %d\n %s\nerror %d: %s\nterminating!\n",__FILE__, __LINE__,err, \
          _hipGetErrorEnum(err)); \*/
        fprintf(stderr, "HIPFFT error %d: %s\nterminating!\n",err,_hipGetErrorEnum(err)); \
            hipDeviceReset(); return; \
    }
}

__global__ void debug(int varId, int N, HIP_DATA_TYPE_COMPLEX *x) {
    //printf("Hello from GPU\n");
    for (int i = 0; i < N; i++)
    {
        HIP_DATA_TYPE_COMPLEX a = x[i];
        double b = (double)a.x;
        double c = (double)a.y;
        if (varId == 0) printf("GPU: input[%d]=(%2.4f,%2.4f)\n",i+1,b,c);
        if (varId == 1) printf("GPU: output[%d]=(%2.4f,%2.4f)\n",i+1,b,c);
    }}

__global__ void debugFloat(int varId, int N, HIP_DATA_TYPE_REAL *x) {
    //printf("Hello from GPU\n");
    for (int i = 0; i < N; i++)
    {
        double a = (double)x[i];
        if (varId == 0) printf("GPU: input[%d]=%2.4f\n",i+1,a);
        if (varId == 1) printf("GPU: output[%d]=%2.4f\n",i+1,a);
    }}

extern "C" {

    void 
        execute_plan_ffth_c_(int ISIGNp, int N, DATA_TYPE *data_in_host, DATA_TYPE *data_out_host, long *iplan)
        //void hipfunction(int ISIGNp, int N, DATA_TYPE *data_in_host, DATA_TYPE *data_out_host, long *iplan)
        {
            HIP_DATA_TYPE_COMPLEX *data_in = reinterpret_cast<HIP_DATA_TYPE_COMPLEX*>(data_in_host);
            HIP_DATA_TYPE_COMPLEX *data_out = reinterpret_cast<HIP_DATA_TYPE_COMPLEX*>(data_out_host);
            hipfftHandle* PLANp = reinterpret_cast<hipfftHandle*>(iplan);
            //fprintf(stderr, "execute_plan_ffth_c_: plan-address = %p\n",PLANp);
            //abort();
            hipfftHandle plan = *PLANp;
            int ISIGN = ISIGNp;

            // Check variables on the GPU:
            /*int device_count = 0;
              hipGetDeviceCount(&device_count);
              for (int i = 0; i < device_count; ++i) {
              hipSetDevice(i);
              hipLaunchKernelGGL(debug, dim3(1), dim3(1), 0, 0, 0, N, data_in);
              hipDeviceSynchronize();
              }*/

            /*if (hipDeviceSynchronize() != hipSuccess){
              fprintf(stderr, "Hip error: Failed to synchronize\n");
              return;	
              }*/

            if( ISIGN== -1 ){
#ifdef TRANS_SINGLE
                hipfftSafeCall(hipfftExecR2C(plan, (HIP_DATA_TYPE_REAL*)data_in, data_out));
#else
                hipfftSafeCall(hipfftExecD2Z(plan, (HIP_DATA_TYPE_REAL*)data_in, data_out));
#endif	
            }
            else if( ISIGN== 1){
#ifdef TRANS_SINGLE
                hipfftSafeCall(hipfftExecC2R(plan, data_in, (HIP_DATA_TYPE_REAL*)data_out));
#else
                hipfftSafeCall(hipfftExecZ2D(plan, data_in, (HIP_DATA_TYPE_REAL*)data_out));
#endif	
            }
            else {
                abort();
            }

            hipDeviceSynchronize();

            /*for (int i = 0; i < device_count; ++i) {
              hipSetDevice(i);
              hipLaunchKernelGGL(debugFloat, dim3(1), dim3(1), 0, 0, 1, N, (HIP_DATA_TYPE_REAL*)data_out);
              hipDeviceSynchronize();
              }*/

            //if (hipDeviceSynchronize() != hipSuccess){
            //	fprintf(stderr, "Hip error: Failed to synchronize\n");
            //	return;	
            //}


        }

}

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

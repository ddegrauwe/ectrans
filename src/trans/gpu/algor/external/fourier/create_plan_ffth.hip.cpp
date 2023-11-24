#define hipfftSafeCall(err) __hipfftSafeCall(err, __FILE__, __LINE__)
#include "stdio.h"
#include <hip/hip_runtime.h>
#include "hipfft.h"
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
			fprintf(stderr, "HIPFFT error in file '%s', line %d; error %d: %s\nterminating!\n", \
				__FILE__, __LINE__,err, \
				_hipGetErrorEnum(err)); \
			fprintf(stderr, "HIPFFT error %d: %s\nterminating!\n",err,_hipGetErrorEnum(err)); \
			hipDeviceReset(); abort(); \
		}
    }


static int allocatedWorkspace=0;
static void* planWorkspace;
static int planWorkspaceSize=100*1024*1024; //100MB
 
extern "C"
void
create_plan_ffth_(hipfftHandle * * plan_ptr_ptr, int *ISIGNp, int *Np, int *LOTp, int * NONSTRIDEDp)
{
int ISIGN = *ISIGNp;
int N = *Np;
int LOT = *LOTp;
int NONSTRIDED = *NONSTRIDEDp;

if (hipDeviceSynchronize() != hipSuccess){
	fprintf(stderr, "Hip error: Failed to synchronize\n");
	return;	
}

hipfftHandle * plan_ptr = new hipfftHandle;
plan_ptr_ptr[0]=plan_ptr;

// create plan
hipfftSafeCall(hipfftCreate(plan_ptr));

// disable auto allocation so we can re-use a single workspace (created above)
hipfftSafeCall(hipfftSetAutoAllocation(*plan_ptr, false));

//create a single re-usable workspace
if(!allocatedWorkspace){
  allocatedWorkspace=1;
  //allocate plan workspace
  hipMalloc(&planWorkspace,planWorkspaceSize);
}


// plan parameters
int embed[1];
int stride;
int cdist, rdist;

#ifdef TRANS_SINGLE
hipfftType hipfft_1 = HIPFFT_R2C;
hipfftType hipfft_2 = HIPFFT_C2R;
#else
hipfftType hipfft_1 = HIPFFT_D2Z;
hipfftType hipfft_2 = HIPFFT_Z2D;
#endif

embed[0] = 0;
if (NONSTRIDED==0) {
	// for global and LAM zonal
	stride   = LOT;
	cdist     = 1;
	rdist     = 1;
} else {
	// for LAM meridional
	stride=1;
	cdist=N/2+1;
	rdist=N+2;
}


fprintf(stderr,"CreatePlan hipfft for \n");
fprintf(stderr,"  %s %p \n","plan address=",plan_ptr);
fprintf(stderr,"  %s %d \n","LOT=",LOT);
fprintf(stderr,"  %s %d \n","stride=",stride);
fprintf(stderr,"  %s %d \n","rdist=",rdist);
fprintf(stderr,"  %s %d \n","cdist=",cdist);
fprintf(stderr,"  %s %p \n","embed address=",embed);
fprintf(stderr,"  %s %d \n","ISIGN=",ISIGN);
fprintf(stderr,"  %s %d \n","N=",N);
fprintf(stderr,"  %s %p \n","N address=",&N);

size_t workSize=123456;

if( ISIGN== -1 ){
  hipfftSafeCall(hipfftMakePlanMany(*plan_ptr, 1, &N,
                 embed, stride, rdist, 
                 embed, stride, cdist, 
                 hipfft_1, LOT, &workSize));
}
else if( ISIGN== 1){
  hipfftSafeCall(hipfftMakePlanMany(*plan_ptr, 1, &N,
                 embed, stride, cdist, 
                 embed, stride, rdist, 
                 hipfft_2, LOT, &workSize));
}
else {
  abort();
}

// use our reusaable work area for the plan
hipfftSafeCall(hipfftSetWorkArea(*plan_ptr,planWorkspace));

// print worksize returned from hipfftMakePlan
fprintf(stderr,"  %s %d \n","workSize from hipfftMakePlan=",workSize);

// get worksize
hipfftSafeCall(hipfftGetSize(*plan_ptr, &workSize));
fprintf(stderr,"  %s %d \n","workSize from hipfftGetSize=",workSize);

// abort if we don't have enough space for the work area in the re-usable workspace
if(workSize > planWorkspaceSize){
	fprintf(stderr,"create_plan_ffth: plan workspace size not large enough - aborting\n");
	abort();
}

if (hipDeviceSynchronize() != hipSuccess){
	fprintf(stderr, "Hip error: Failed to synchronize\n");
	return;	
}

return;


}


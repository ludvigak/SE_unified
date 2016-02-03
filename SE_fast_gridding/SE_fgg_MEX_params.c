#include "SE_fgg.h"
#include "mex.h" 

// Get field s from matlab struct p. abort if field is missing
void* get_arg(const mxArray* p, const char* s)
{
    mxArray* x=mxGetField(p,0,s);
    if(x) return mxGetData(x);
    else
    {
	mexErrMsgTxt("Missing mandatory parameter in struct");
	return (void*) NULL;
    }
}

// get relevant parameters from mxArray and populate FGG params struct
void SE_FGG_MEX_params(SE_FGG_params* params, const mxArray* OPT, int N)
{

#ifdef THREE_PERIODIC

    const double* m    = (double*) get_arg(OPT,"M");
    const double* p    = (double*) get_arg(OPT,"P");
    const double* c    = (double*) get_arg(OPT,"c");
    const double* box  = (double*) get_arg(OPT,"box");

    params->N = N;
    params->P = (int) *p;
    params->P_half=half( (int) *p );
    params->c = *c;
    params->d = pow(params->c/PI,1.5);
    params->h = box[0]/m[0];
    params->a = -FGG_INF;

    params->dims[0] = (int) m[0];
    params->dims[1] = (int) m[1];
    params->dims[2] = (int) m[2];

    params->npdims[0] = params->dims[0]+params->P;
    params->npdims[1] = params->dims[1]+params->P;
    params->npdims[2] = params->dims[2]+params->P;

#endif

#ifdef TWO_PERIODIC

    const double* m    = (double*) get_arg(OPT,"M");
    const double* mz   = (double*) get_arg(OPT,"Mz");
    const double* p    = (double*) get_arg(OPT,"P");
    const double* c    = (double*) get_arg(OPT,"c");
    const double* box  = (double*) get_arg(OPT,"box");
    const double* a    = (double*) get_arg(OPT,"a"); /* z-dir offset. RENAME */

    params->N = N;
    params->P = (int) *p;
    params->P_half=half( (int) *p );
    params->c = *c;
    params->d = pow(params->c/PI,1.5);
    params->h = box[0]/m[0];
    params->a = a[0];

    params->dims[0] = (int)  m[0];
    params->dims[1] = (int)  m[0];
    params->dims[2] = (int) mz[0];

    params->npdims[0] = params->dims[0]+params->P;
    params->npdims[1] = params->dims[1]+params->P;
    params->npdims[2] = params->dims[2];

#endif
}

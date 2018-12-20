//// HEADER GUARD ///////////////////////////
// If automatically generated, keep above
// comment as first line in file.  
#ifndef __MRG_MYSA_H__
#define __MRG_MYSA_H__
//// HEADER GUARD ///////////////////////////
// DO NOT EDIT THIS SOURCE CODE FILE
// ANY CHANGES TO THIS FILE WILL BE OVERWRITTEN!!!!
//   If you need to make a change, do what you have to do and send Rob
//   an email with what you had to change and why.  He'll fix the translator
//   so your fix will be recorded the next time the translator runs.

#define MRG_MYSA_REQDAT Vm_DATA_FLAG
#define MRG_MYSA_MODDAT Iion_DATA_FLAG

typedef struct MRG_MYSA_params {

} MRG_MYSA_Params;

typedef struct MRG_MYSA_state {
    char dummy; //PGI doesn't allow a structure of 0 bytes.

} MRG_MYSA_state;

#ifdef __CUDACC__
extern "C" {
#endif

void initialize_params_MRG_MYSA(ION_IF *);
void construct_tables_MRG_MYSA(ION_IF *);
void destroy_MRG_MYSA(ION_IF *);
void initialize_sv_MRG_MYSA(ION_IF *, GlobalData_t**);
GLOBAL void compute_MRG_MYSA(int, int, ION_IF *, GlobalData_t**);

#ifdef __CUDACC__
};
#endif

//// HEADER GUARD ///////////////////////////
#endif
//// HEADER GUARD ///////////////////////////

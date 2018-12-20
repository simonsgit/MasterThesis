//// HEADER GUARD ///////////////////////////
// If automatically generated, keep above
// comment as first line in file.  
#ifndef __MRG_NODE_H__
#define __MRG_NODE_H__
//// HEADER GUARD ///////////////////////////
// DO NOT EDIT THIS SOURCE CODE FILE
// ANY CHANGES TO THIS FILE WILL BE OVERWRITTEN!!!!
//   If you need to make a change, do what you have to do and send Rob
//   an email with what you had to change and why.  He'll fix the translator
//   so your fix will be recorded the next time the translator runs.

#define MRG_Node_REQDAT Vm_DATA_FLAG
#define MRG_Node_MODDAT Iion_DATA_FLAG

typedef struct MRG_Node_params {

} MRG_Node_Params;

typedef struct MRG_Node_state {
    Gatetype h;
    Gatetype m;
    Gatetype p;
    Gatetype s;

} MRG_Node_state;

#ifdef __CUDACC__
extern "C" {
#endif

void initialize_params_MRG_Node(ION_IF *);
void construct_tables_MRG_Node(ION_IF *);
void destroy_MRG_Node(ION_IF *);
void initialize_sv_MRG_Node(ION_IF *, GlobalData_t**);
GLOBAL void compute_MRG_Node(int, int, ION_IF *, GlobalData_t**);

#ifdef __CUDACC__
};
#endif

//// HEADER GUARD ///////////////////////////
#endif
//// HEADER GUARD ///////////////////////////

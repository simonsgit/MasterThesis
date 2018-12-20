
// DO NOT EDIT THIS SOURCE CODE FILE
// ANY CHANGES TO THIS FILE WILL BE OVERWRITTEN!!!!
//   If you need to make a change, do what you have to do and send Rob
//   an email with what you had to change and why.  He'll fix the translator
//   so your fix will be recorded the next time the translator runs.

#include "ION_IF.h"
#include "NumComp.h"
#include "MRG_MYSA.h"


#ifdef __CUDACC__
#define pow powf
#define log logf
#endif

void trace_MRG_MYSA(ION_IF* IF, int node, FILE* file, GlobalData_t** impdata);

void destroy_MRG_MYSA( ION_IF *IF )
{
  destroy_luts( IF );
  SV_free( &IF->sv_tab );
  // rarely need to do anything else
}

// Define all constants
#define Cm (GlobalData_t)(2.)
#define E_R (GlobalData_t)(-80.)
#define V_init (GlobalData_t)(-80.)
#define g_L (GlobalData_t)(1.)
#define E_L (GlobalData_t)((E_R-(10.)))



void initialize_params_MRG_MYSA( ION_IF *IF )
{
  cell_geom *region = &IF->cgeom;
  MRG_MYSA_Params *p = (MRG_MYSA_Params *)IF->params;

  // Compute the regional constants
  {
  }
  // Compute the regional initialization
  {
  }

}


// Define the parameters for the lookup tables
enum Tables {
  V_TAB,

  N_TABS
};

// Define the indices into the lookup tables.
enum V_TableIndex {
  Iion_idx,
  NROWS_V
};



void construct_tables_MRG_MYSA( ION_IF *IF )
{
  GlobalData_t dt = IF->dt*1e0;
  cell_geom *region = &IF->cgeom;
  MRG_MYSA_Params *p = (MRG_MYSA_Params *)IF->params;

  IF->numLUT = N_TABS;
  IF->tables = (LUT *)IMP_malloc( N_TABS, sizeof(LUT) );

  // Define the constants that depend on the parameters.
  
  // Create the V lookup table
  LUT* V_tab = &IF->tables[V_TAB];
  LUT_alloc(V_tab, NROWS_V, -200, 200, 0.05, "MRG_MYSA V");    
  for (int __i=V_tab->mn_ind; __i<=V_tab->mx_ind; __i++) {
    double V = V_tab->res*__i;
    LUT_data_t* V_row = V_tab->tab[__i];
    double i_L = (g_L*(V-(E_L)));
    V_row[Iion_idx] = (i_L/Cm);
  }
  check_LUT(V_tab);
  

}



void    initialize_sv_MRG_MYSA( ION_IF *IF, GlobalData_t **impdata )
{
  GlobalData_t dt = IF->dt*1e0;
  cell_geom *region = &IF->cgeom;
  MRG_MYSA_Params *p = (MRG_MYSA_Params *)IF->params;

  SV_alloc( &IF->sv_tab, IF->numNode, sizeof(MRG_MYSA_state) );
  MRG_MYSA_state *sv_base = (MRG_MYSA_state *)IF->sv_tab.y;
  GlobalData_t t = 0;
  // Define the constants that depend on the parameters.
  //Prepare all the public arrays.
  GlobalData_t *Iion_ext = impdata[Iion];
  GlobalData_t *V_ext = impdata[Vm];
  //Prepare all the private functions.

  //set the initial values
  for(int __i=0; __i<IF->sv_tab.numSeg; __i++ ){
    MRG_MYSA_state *sv = sv_base+__i;

    // Initialize nodal variables that have been declared with param
    //Initialize the external vars to their current values
    GlobalData_t Iion = Iion_ext[__i];
    GlobalData_t V = V_ext[__i];
    //Change the units of external variables as appropriate.
    
    
    // Initialize the rest of the nodal variables
    V = V_init;
    double i_L = (g_L*(V-(E_L)));
    Iion = (i_L/Cm);
    //Change the units of external variables as appropriate.
    
    
    //Save all external vars
    Iion_ext[__i] = Iion;
    V_ext[__i] = V;

  }

}

/** compute the  current
 *
 * param start   index of first node
 * param end     index of last node
 * param IF      IMP
 * param plgdata external data needed by IMP
 */
GLOBAL void compute_MRG_MYSA(int start, int end, ION_IF *IF, GlobalData_t **impdata )
{
  GlobalData_t dt = IF->dt*1e0;
  cell_geom *region = &IF->cgeom;
  MRG_MYSA_Params *p  = (MRG_MYSA_Params *)IF->params;
  MRG_MYSA_state *sv_base = (MRG_MYSA_state *)IF->sv_tab.y;

  GlobalData_t t = IF->tstp.cnt*dt;

  // Define the constants that depend on the parameters.
  //Prepare all the public arrays.
  GlobalData_t *Iion_ext = impdata[Iion];
  GlobalData_t *V_ext = impdata[Vm];
  //Prepare all the private functions.

#ifdef __CUDACC__
  int __i = blockDim.x * blockIdx.x + threadIdx.x;
  if ( __i>=end ) { return; }
  {
#else
#pragma omp parallel for schedule(static)
  for (int __i=start; __i<end; __i++) {
#endif
    MRG_MYSA_state *sv = sv_base+__i;

    //Initialize the external vars to their current values
    GlobalData_t Iion = Iion_ext[__i];
    GlobalData_t V = V_ext[__i];
    //Change the units of external variables as appropriate.
    
    
    //Compute lookup tables for things that have already been defined.
    LUT_data_t V_row[NROWS_V];
    LUT_interpRow(&IF->tables[V_TAB], V, __i, V_row);
    
    
    //Compute storevars and external modvars
    Iion = V_row[Iion_idx];
    
    
    //Complete Forward Euler Update
    
    
    //Complete Rush Larsen Update
    
    
    //Complete RK2 Update
    {
      GlobalData_t t = t + dt/2;
    }
    
    
    //Complete RK4 Update
    {
      GlobalData_t t = t + dt/2;
    }
    {
      GlobalData_t t = t + dt/2;
    }
    {
      GlobalData_t t = t + dt;
    }
    
    
    //Complete Sundnes Update
    
    
    //Complete Markov Backward Euler method
    
    
    //Complete Rosenbrock Update
    
    
    //Complete Cvode Update
    
    
    //Finish the update
    if (0) {}
    else {Iion = Iion;}
    //Change the units of external variables as appropriate.
    
    
    //Save all external vars
    Iion_ext[__i] = Iion;
    V_ext[__i] = V;

  }    


}

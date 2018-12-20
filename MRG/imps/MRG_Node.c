
// DO NOT EDIT THIS SOURCE CODE FILE
// ANY CHANGES TO THIS FILE WILL BE OVERWRITTEN!!!!
//   If you need to make a change, do what you have to do and send Rob
//   an email with what you had to change and why.  He'll fix the translator
//   so your fix will be recorded the next time the translator runs.

#include "ION_IF.h"
#include "NumComp.h"
#include "MRG_Node.h"


#ifdef __CUDACC__
#define pow powf
#define log logf
#endif

void trace_MRG_Node(ION_IF* IF, int node, FILE* file, GlobalData_t** impdata);

void destroy_MRG_Node( ION_IF *IF )
{
  destroy_luts( IF );
  SV_free( &IF->sv_tab );
  // rarely need to do anything else
}

// Define all constants
#define Cm (GlobalData_t)(2.)
#define E_R (GlobalData_t)(-80.)
#define V_init (GlobalData_t)(-80.)
#define g_K (GlobalData_t)(80.)
#define g_L (GlobalData_t)(20.)
#define g_Naf (GlobalData_t)(3000.)
#define g_Nap (GlobalData_t)(10.)
#define E_K (GlobalData_t)((E_R-(10.)))
#define E_L (GlobalData_t)((E_R-(10.)))
#define E_Na (GlobalData_t)((E_R+130.))



void initialize_params_MRG_Node( ION_IF *IF )
{
  cell_geom *region = &IF->cgeom;
  MRG_Node_Params *p = (MRG_Node_Params *)IF->params;

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
  h_rush_larsen_A_idx,
  h_rush_larsen_B_idx,
  i_L_idx,
  m_rush_larsen_A_idx,
  m_rush_larsen_B_idx,
  p_rush_larsen_A_idx,
  p_rush_larsen_B_idx,
  s_rush_larsen_A_idx,
  s_rush_larsen_B_idx,
  NROWS_V
};



void construct_tables_MRG_Node( ION_IF *IF )
{
  GlobalData_t dt = IF->dt*1e0;
  cell_geom *region = &IF->cgeom;
  MRG_Node_Params *p = (MRG_Node_Params *)IF->params;

  IF->numLUT = N_TABS;
  IF->tables = (LUT *)IMP_malloc( N_TABS, sizeof(LUT) );

  // Define the constants that depend on the parameters.
  
  // Create the V lookup table
  LUT* V_tab = &IF->tables[V_TAB];
  LUT_alloc(V_tab, NROWS_V, -200, 200, 0.05, "MRG_Node V");    
  for (int __i=V_tab->mn_ind; __i<=V_tab->mx_ind; __i++) {
    double V = V_tab->res*__i;
    LUT_data_t* V_row = V_tab->tab[__i];
    double a_h = (((fabs(((V+114.)/11.)))<0.000001) ? (0.38*11.) : ((0.38*(-(V+114.)))/(1.-((exp(((V+114.)/11.)))))));
    double a_m = (((fabs(((-(V+21.4))/10.3)))<0.000001) ? (7.11*10.3) : ((7.11*(V+21.4))/(1.-((exp(((-(V+21.4))/10.3)))))));
    double a_p = (((fabs(((V+27.)/10.2)))<0.000001) ? (0.0382*10.2) : ((0.0382*(V+27.))/(1.-((exp(((-(V+27.))/10.2)))))));
    double a_s = (0.3348/(1.+(exp(((V+53.)/-5.)))));
    double b_h = (14.1/(1.+(exp(((V+31.8)/-13.4)))));
    double b_m = (((fabs(((V+25.7)/9.16)))<0.000001) ? (0.329*9.16) : ((0.329*(-(V+25.7)))/(1.-((exp(((V+25.7)/9.16)))))));
    double b_p = (((fabs(((V+34.)/10.)))<0.000001) ? (0.000955*10.) : ((0.000955*(-(V+34.)))/(1.-((exp(((V+34.)/10.)))))));
    double b_s = (0.03348/(1.+(exp(((V+90.)/-1.)))));
    V_row[i_L_idx] = (g_L*(V-(E_L)));
    V_row[h_rush_larsen_A_idx] = (((-a_h)/(a_h+b_h))*(expm1(((-dt)*(a_h+b_h)))));
    V_row[h_rush_larsen_B_idx] = (exp(((-dt)*(a_h+b_h))));
    V_row[m_rush_larsen_A_idx] = (((-a_m)/(a_m+b_m))*(expm1(((-dt)*(a_m+b_m)))));
    V_row[m_rush_larsen_B_idx] = (exp(((-dt)*(a_m+b_m))));
    V_row[p_rush_larsen_A_idx] = (((-a_p)/(a_p+b_p))*(expm1(((-dt)*(a_p+b_p)))));
    V_row[p_rush_larsen_B_idx] = (exp(((-dt)*(a_p+b_p))));
    V_row[s_rush_larsen_A_idx] = (((-a_s)/(a_s+b_s))*(expm1(((-dt)*(a_s+b_s)))));
    V_row[s_rush_larsen_B_idx] = (exp(((-dt)*(a_s+b_s))));
  }
  check_LUT(V_tab);
  

}



void    initialize_sv_MRG_Node( ION_IF *IF, GlobalData_t **impdata )
{
  GlobalData_t dt = IF->dt*1e0;
  cell_geom *region = &IF->cgeom;
  MRG_Node_Params *p = (MRG_Node_Params *)IF->params;

  SV_alloc( &IF->sv_tab, IF->numNode, sizeof(MRG_Node_state) );
  MRG_Node_state *sv_base = (MRG_Node_state *)IF->sv_tab.y;
  GlobalData_t t = 0;
  // Define the constants that depend on the parameters.
  //Prepare all the public arrays.
  GlobalData_t *Iion_ext = impdata[Iion];
  GlobalData_t *V_ext = impdata[Vm];
  //Prepare all the private functions.

  //set the initial values
  for(int __i=0; __i<IF->sv_tab.numSeg; __i++ ){
    MRG_Node_state *sv = sv_base+__i;

    // Initialize nodal variables that have been declared with param
    //Initialize the external vars to their current values
    GlobalData_t Iion = Iion_ext[__i];
    GlobalData_t V = V_ext[__i];
    //Change the units of external variables as appropriate.
    
    
    // Initialize the rest of the nodal variables
    V = V_init;
    double a_h = (((fabs(((V+114.)/11.)))<0.000001) ? (0.38*11.) : ((0.38*(-(V+114.)))/(1.-((exp(((V+114.)/11.)))))));
    double a_m = (((fabs(((-(V+21.4))/10.3)))<0.000001) ? (7.11*10.3) : ((7.11*(V+21.4))/(1.-((exp(((-(V+21.4))/10.3)))))));
    double a_p = (((fabs(((V+27.)/10.2)))<0.000001) ? (0.0382*10.2) : ((0.0382*(V+27.))/(1.-((exp(((-(V+27.))/10.2)))))));
    double a_s = (0.3348/(1.+(exp(((V+53.)/-5.)))));
    double b_h = (14.1/(1.+(exp(((V+31.8)/-13.4)))));
    double b_m = (((fabs(((V+25.7)/9.16)))<0.000001) ? (0.329*9.16) : ((0.329*(-(V+25.7)))/(1.-((exp(((V+25.7)/9.16)))))));
    double b_p = (((fabs(((V+34.)/10.)))<0.000001) ? (0.000955*10.) : ((0.000955*(-(V+34.)))/(1.-((exp(((V+34.)/10.)))))));
    double b_s = (0.03348/(1.+(exp(((V+90.)/-1.)))));
    double i_L = (g_L*(V-(E_L)));
    double h_init = (a_h/(a_h+b_h));
    double m_init = (a_m/(a_m+b_m));
    double p_init = (a_p/(a_p+b_p));
    double s_init = (a_s/(a_s+b_s));
    sv->h = h_init;
    sv->m = m_init;
    sv->p = p_init;
    sv->s = s_init;
    double i_K = ((g_K*sv->s)*(V-(E_K)));
    double i_Naf = (((g_Naf*((sv->m*sv->m)*sv->m))*sv->h)*(V-(E_Na)));
    double i_Nap = ((g_Nap*((sv->p*sv->p)*sv->p))*(V-(E_Na)));
    Iion = (-((-(((i_Naf+i_Nap)+i_K)+i_L))/Cm));
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
GLOBAL void compute_MRG_Node(int start, int end, ION_IF *IF, GlobalData_t **impdata )
{
  GlobalData_t dt = IF->dt*1e0;
  cell_geom *region = &IF->cgeom;
  MRG_Node_Params *p  = (MRG_Node_Params *)IF->params;
  MRG_Node_state *sv_base = (MRG_Node_state *)IF->sv_tab.y;

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
    MRG_Node_state *sv = sv_base+__i;

    //Initialize the external vars to their current values
    GlobalData_t Iion = Iion_ext[__i];
    GlobalData_t V = V_ext[__i];
    //Change the units of external variables as appropriate.
    
    
    //Compute lookup tables for things that have already been defined.
    LUT_data_t V_row[NROWS_V];
    LUT_interpRow(&IF->tables[V_TAB], V, __i, V_row);
    
    
    //Compute storevars and external modvars
    GlobalData_t i_K = ((g_K*sv->s)*(V-(E_K)));
    GlobalData_t i_Naf = (((g_Naf*((sv->m*sv->m)*sv->m))*sv->h)*(V-(E_Na)));
    GlobalData_t i_Nap = ((g_Nap*((sv->p*sv->p)*sv->p))*(V-(E_Na)));
    Iion = (-((-(((i_Naf+i_Nap)+i_K)+V_row[i_L_idx]))/Cm));
    
    
    //Complete Forward Euler Update
    
    
    //Complete Rush Larsen Update
    GlobalData_t h_rush_larsen_A = V_row[h_rush_larsen_A_idx];
    GlobalData_t h_rush_larsen_B = V_row[h_rush_larsen_B_idx];
    GlobalData_t m_rush_larsen_A = V_row[m_rush_larsen_A_idx];
    GlobalData_t m_rush_larsen_B = V_row[m_rush_larsen_B_idx];
    GlobalData_t p_rush_larsen_A = V_row[p_rush_larsen_A_idx];
    GlobalData_t p_rush_larsen_B = V_row[p_rush_larsen_B_idx];
    GlobalData_t s_rush_larsen_A = V_row[s_rush_larsen_A_idx];
    GlobalData_t s_rush_larsen_B = V_row[s_rush_larsen_B_idx];
    GlobalData_t h_new = h_rush_larsen_A+h_rush_larsen_B*sv->h;
    GlobalData_t m_new = m_rush_larsen_A+m_rush_larsen_B*sv->m;
    GlobalData_t p_new = p_rush_larsen_A+p_rush_larsen_B*sv->p;
    GlobalData_t s_new = s_rush_larsen_A+s_rush_larsen_B*sv->s;
    
    
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
    if (0) {}
    else {sv->h = h_new;}
    if (0) {}
    else {sv->m = m_new;}
    if (0) {}
    else {sv->p = p_new;}
    if (0) {}
    else {sv->s = s_new;}
    //Change the units of external variables as appropriate.
    
    
    //Save all external vars
    Iion_ext[__i] = Iion;
    V_ext[__i] = V;

  }    


}

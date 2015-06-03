#include<omp.h>
#include <iostream>
#include <complex>
#include <stdio.h>

// Run it with g++ -std=gnu++11 hrchl.cpp -lm -o hrc

#include<vector>
#include<math.h>
#include "hrchl_tree.h"
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>

using namespace std;



int main(){
    
    
    // inputs from tree structure
    int Dim, NUM_shls_tree, NUM_shls_goy, N_nds_tree, N_nds_goy, N_nds;
    
    double alpha,k0,g,alphabar,q0,Fp,mubg,musm,mubgFac,musmFac;
    int fi; // shell where to apply force.  
    Physics phys;
    Tree tree;
    
    int i,n, j,status,Nav;

    complex<double> I (0.0,1.0);
    
    FILE *rs,*rsZF,*rscmp,*inpPlot;
    
    rs=fopen("rslts.txt","w");
    rsZF=fopen("rsltsZF.txt","w");
    rscmp=fopen("rsltscmp.txt","w");
    inpPlot=fopen("inputsPlot.txt","w"); //creating file for inputs for plot
    
    
    FILE *fInp;
    /*** values from input file ****/


    double kmax;
    

    
    fInp = fopen("INPUT","r");
    Dim = (int) get_data_NAME(fInp,"Dimension");
    fclose(fInp);
    
    fInp = fopen("INPUT","r");
    NUM_shls_tree =(int) get_data_NAME(fInp,"Number of layers in tree"); 
    fclose(fInp);
    
    fInp = fopen("INPUT","r");
    NUM_shls_goy = (int) get_data_NAME(fInp,"Number of layers as shell"); 
    fclose(fInp);
    
    fInp = fopen("INPUT","r");
    k0 = get_data_NAME(fInp,"k min");    
    fclose(fInp);
    
    fInp = fopen("INPUT","r");
    kmax = get_data_NAME(fInp,"k max");    
    fclose(fInp);
    
   
    
    fInp = fopen("INPUT","r");
    musmFac = (double) get_data_NAME(fInp,"musmFac");   
    fclose(fInp);
    
    fInp = fopen("INPUT","r");
    mubgFac = (double) get_data_NAME(fInp,"mubgFac");
    fclose(fInp);
    
    fInp = fopen("INPUT","r");
    fi = get_data_NAME(fInp,"fi shell");    
    fclose(fInp);
    
    fInp = fopen("INPUT","r");
    Fp = get_data_NAME(fInp,"forcing magnitude");
    fclose(fInp);


    N_nds_tree =  (int) ((int) pow( (double)Dim,(NUM_shls_tree))-1)/(Dim-1);
    N_nds_goy = NUM_shls_goy*((int) pow((double)Dim,NUM_shls_tree-1));  
    N_nds = N_nds_tree + N_nds_goy;
    g = exp(log(kmax/k0)/((double) NUM_shls_tree+NUM_shls_goy-1));
    alpha = g*g;

    musm = musmFac*pow(k0,6);     
    mubg = mubgFac*pow(k0*pow(g,NUM_shls_goy+NUM_shls_tree-1),-4)
    ;
    

    
    alphabar = 0; q0 = 0.01;
    Nav = 1;
    
    
    tree.setupTree(Dim,NUM_shls_tree,NUM_shls_goy);
    phys.setupPhysics(alpha,k0,g,alphabar,q0,mubg,musm,Fp,fi,tree); 

    
   
    ParamForODE prms = {&tree,&phys};
    
    
    
    
    int dimension = sizeof(complex<double>)*N_nds;
    
    double eps_abs = 1.0e-12;     // absolute err requested                     
    double eps_rel = 1.e-6;       // relative error requested                     
    
    // define the type of routine for making steps:                               
    
    
    const gsl_odeiv2_step_type *type_ptr = gsl_odeiv2_step_rk8pd;
    
    
    gsl_odeiv2_step *step_ptr = gsl_odeiv2_step_alloc (type_ptr, dimension);
    gsl_odeiv2_control *control_ptr  = gsl_odeiv2_control_standard_new (eps_abs, eps_rel, 1.0,0.0);
    gsl_odeiv2_evolve *evolve_ptr = gsl_odeiv2_evolve_alloc (dimension);
    gsl_odeiv2_system my_system;
    
    double t, t_next,delta_t;   // current and next independent variable                                                                      
    double tmin, tmax;          // range of t and step size for output                     
    
    
    double h = 1.e-12;            // starting step size for ode solver            
    
    int szy;
    
    
    szy = sizeof(complex<double>)*N_nds;
    double y[szy];
    
    my_system.function = func_hmdsi;      // the right-hand-side functions dy[i]/dt                                 
    my_system.jacobian = NULL;    // the Jacobian df[i]/dy[j]                     
    my_system.dimension = dimension;      // number of diffeq's                   
    my_system.params = &prms;       // parameters to pass to rhs and jacobian       
    


    fInp = fopen("INPUT","r");
    delta_t = get_data_NAME(fInp,"dt");
    fclose(fInp);
    fInp = fopen("INPUT","r");
    tmax = get_data_NAME(fInp,"tmax");
    fclose(fInp);

    tmin = 0;
    t=tmin;
    complex<double> *phi = (complex<double> *)&y[0];
    double phiAv[N_nds];
    int v;
    v=0;
    for(i=0;i<N_nds;i++){ 
        phi[i] =1.e-12*exp(rand()*2.*M_PI*I);
        phiAv[i] = 0;            
    }
    
    
    for (t_next = tmin + delta_t; t_next <= tmax; t_next += delta_t)
    {
        v = v+1;
        while (t < t_next)        // evolve from t to t_next       
        {
            status= gsl_odeiv2_evolve_apply (evolve_ptr, control_ptr, step_ptr, &my_system, &t, t_next, &h, y);
            
            if (status != GSL_SUCCESS)
                break;
        }
        if(v>Nav){
            v=0;
        }
        std::scientific;
        printf("%e\n",t/tmax*100.,"\%");
        
        n=0;
        
        std::scientific;
        for(i=0;i<N_nds;i++){
            
            phiAv[i] =  ((double) abs(phi[i])*abs(phi[i]))/((double) Nav);
            
            if((i>phys.k_idx_max[n-1]) || i==0){
                fprintf(rs,"%e \t %e \n", phys.k_n[n], abs(phi[i])*abs(phi[i]));
                n = n+1;
            }
        }
        
        //Print the entire wavelet to file: averaged over Nav;
        if(v==0){
            for(i=phys.k_idx_min[NUM_shls_tree+NUM_shls_goy-1];i<phys.k_idx_max[NUM_shls_tree+NUM_shls_goy-1]+1;i++){
                j=i;
                while(tree.nds[j].has_pr!=0){
                    fprintf(rscmp,"%e \t", phiAv[j]);
                    j=tree.nds[j].i_pr;
                }
                fprintf(rscmp,"%e \n", phiAv[j]);
                
                
            }
            
            for(i=0;i<N_nds;i++)
                phiAv[i] = 0;
        }
        
        
    }
    
    /// input for plotting 
    
    int iniT,NtAv,LstN;
    double allT;
    
    fInp = fopen("INPUT","r");
    iniT = (int) get_data_NAME(fInp,"t1");    
    fclose(fInp);
    
    fInp = fopen("INPUT","r");
    NtAv = (int) get_data_NAME(fInp,"Nav");    
    fclose(fInp);
    
    fInp = fopen("INPUT","r");
    LstN = (int) get_data_NAME(fInp,"LastN");    
    fclose(fInp);
    
    allT = tmax/delta_t;
    
    
    if(LstN !=-1){
        iniT = (allT - LstN ? LstN<allT   : 0);
        NtAv = (LstN ? LstN<allT  : allT);
    }
    
    
    fprintf(inpPlot,"%d \t %d \t %d \t",NUM_shls_goy+NUM_shls_tree,iniT,NtAv);
    
    // all done; free up the gsl_odeiv stuff                                                                                 
    gsl_odeiv2_evolve_free (evolve_ptr);
    gsl_odeiv2_control_free (control_ptr);
    gsl_odeiv2_step_free (step_ptr);
    
    
    
    //    printTree(phys,phys.nds);
    
    
    fclose(inpPlot);
    fclose(rsZF);
    fclose(rscmp);
    fclose(rs);
    tree.delete_tree();
    phys.delete_physics();
    
    return 0;
}

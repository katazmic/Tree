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
    
    double alpha,k0,g,alphabar,q0,Fp,mubg,musm;
    int fi; // shell where to apply force.  
    Physics phys;
    Tree tree;
    
    int i,n, j,status,Nav;

    complex<double> I (0.0,1.0);
    
    FILE *rs,*rsZF,*rscmp;
    
    rs=fopen("rslts.txt","w");
    rsZF=fopen("rsltsZF.txt","w");
    rscmp=fopen("rsltscmp.txt","w");
    
    
    
    
    k0 = 1.; g = 1.45; alpha = g*g;alphabar = 0; q0 = 0.01; Fp = .5; 
    //  musm = 1.e3; mubg =1.e-28;

    fi = 8;
    Nav = 1;
    
    
    tree.setupTree();

    
    Dim = tree.Dim;
    NUM_shls_tree = tree.NUM_shls_tree ; // for the fierarchical part! 1 is teh minimum  makes it GOY shell model
    NUM_shls_goy = tree.NUM_shls_goy; // total 30 for reasonable results. 
    
    //com from input
    N_nds_tree =  (int) ((int) pow( (double)Dim,(NUM_shls_tree))-1)/(Dim-1);
    N_nds_goy = NUM_shls_goy*((int) pow((double)Dim,NUM_shls_tree-1));  
    N_nds = N_nds_tree + N_nds_goy;
    
    musm = 10.*pow(k0,6);  1.e-18;//instead of -13.. 10*pow(VIn->kn[0],6);          
    mubg = 10.*pow(k0*pow(g,NUM_shls_goy+NUM_shls_tree-1),-4); 1.e-24; // instead of -15;//100*pow(VIn->kn[VIn->Nshls-1],-4);                   

   phys.setupPhysics(alpha,k0,g,alphabar,q0,mubg,musm,Fp,fi,tree); 
    //for(i=0;i<N_nds;i++)
    //   cout<<find_n(i)<<"\n";
    
    
    
    ParamForODE prms = {&tree,&phys};
    
    
    //prms.tree->printTree(*prms.phys);
    
   // abort();
    
    int dimension = sizeof(complex<double>)*N_nds;
    
    double eps_abs = 1.0e-12;     // absolute err requested                     
    double eps_rel = 1.e-6;       // relative error requested                     
    
    // define the type of routine for making steps:                               
    
    
    const gsl_odeiv2_step_type *type_ptr = gsl_odeiv2_step_rk8pd;
    
    
    gsl_odeiv2_step *step_ptr = gsl_odeiv2_step_alloc (type_ptr, dimension);
    //gsl_odeiv2_control *control_ptr = gsl_odeiv2_control_y_new (eps_abs, eps_rel);                                                                             
    gsl_odeiv2_control *control_ptr  = gsl_odeiv2_control_standard_new (eps_abs, eps_rel, 1.0,0.0);
    gsl_odeiv2_evolve *evolve_ptr = gsl_odeiv2_evolve_alloc (dimension);
    
    gsl_odeiv2_system my_system;
    
    double t, t_next,delta_t;             // current and next independent variable                                                                      
    double tmin, tmax; // range of t and step size for output                     
    
    
    double h = 1.e-12;            // starting step size for ode solver            
    
    int szy;
    
    
    szy = sizeof(complex<double>)*N_nds;
    double y[szy];
    
    my_system.function = func_hmdsi;      // the right-hand-side functions dy[i]/dt                                 
    my_system.jacobian = NULL;    // the Jacobian df[i]/dy[j]                     
    my_system.dimension = dimension;      // number of diffeq's                   
    my_system.params = &prms;       // parameters to pass to rhs and jacobian       
    
    tmin = 0;
    tmax = 200;
    delta_t = 1.e-5;
    
    t=tmin;
    complex<double> *phi = (complex<double> *)&y[0];
    double phiAv[N_nds];
    int v;
    v=0;
    for(i=0;i<N_nds;i++){ 
        phi[i] =1.e-12*exp(rand()*2.*M_PI*I);//*exp(-5*pow(pow(V1.g, (i+1))-pow(V1.g,3),2));              
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
        printf("%0.2f\n",t);
        std::scientific;
        n=0;
        
        std::scientific;
        for(i=0;i<N_nds;i++){
            
            phiAv[i] =  ((double) abs(phi[i])*abs(phi[i]))/((double) Nav);
            
            if((i>phys.k_idx_max[n-1]) || i==0){
                std::cout << phys.k_n[n]<<"\t"<<abs(phi[i])*abs(phi[i]) <<"\n";
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
                    // cout<<phys.k_n[find_n(j)]<<","<<phiAv[j]<<"\t";
                    j=tree.nds[j].i_pr;
                }
                fprintf(rscmp,"%e \n", phiAv[j]);
                //cout<<phys.k_n[find_n(j)]<<","<<phiAv[j]<<"\n";
                
                
            }
            
            for(i=0;i<N_nds;i++)
                phiAv[i] = 0;
        }
        
        
    }
    
    
    // all done; free up the gsl_odeiv stuff                                                                                 
    gsl_odeiv2_evolve_free (evolve_ptr);
    gsl_odeiv2_control_free (control_ptr);
    gsl_odeiv2_step_free (step_ptr);
    
    
    
    //    printTree(phys,phys.nds);
    
    
    fclose(rsZF);
    fclose(rscmp);
    fclose(rs);
    tree.delete_tree();
    phys.delete_physics();
    
    return 0;
}

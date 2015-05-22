//
//  hrchl_tree.h
//  
//
//  Created by katy ghantous on 5/14/15.
//  Copyright 2015 __MyCompanyName__. All rights reserved.
//

#ifndef _hrchl_tree_h
#define _hrchl_tree_h


#include<omp.h>
#include <iostream>
#include <complex>
#include <stdio.h>

// Run it with g++ -std=gnu++11 hrchl.cpp -lm -o hrc

#include<vector>
#include<math.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>


using namespace std;

class Physics;
class Tree;


struct Node{
    int has_pr, has_ch;
    int i_pr;
    int *i_ch;
};


class Tree{
    
public:
    
    int Dim, NUM_shls_tree, NUM_shls_goy;
    int N_nds_tree,N_nds_goy, N_nds;
    
    Node *nds;
    
    int find_n(int i);
    void setupNodes();    
    void setupTree();
    void printTree(Physics phys); 
    void delete_tree(); //frees all the nodes. 
    
    Tree();
    ~Tree();
    
    
};


class Physics{
    
public:
    
    double alpha, beta,gamma, k0,g; //for spectrum
    double alphabar,q0; //for zonal flow
    int fi;
    double Fp;
    double musm,mubg;
    
    int *k_idx_min,*k_idx_max;
    double *k_n,*a_n,*b_n,*c_n;
    
    void setupPhysics(double alphai, double k0i, double gi, double alphabari, double q0i, double mubgi, double musmi, double Fpi, int fii,Tree tree);
    void delete_physics();
    
    Physics();
    ~Physics();
};

struct ParamForODE{
    
    Tree *tree;
    Physics *phys;
};

int func_hmdsi(double t, const double y[], double dydt[], void *params);


#endif

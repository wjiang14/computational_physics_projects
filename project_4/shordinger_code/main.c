//
//  main.c
//  shordinger_eq
//
//  Created by wei on 15/5/2.
//  Copyright (c) 2015年 wei. All rights reserved.
//

#include <stdio.h>
#include "math.h"

double x, x_left, x_right;
int E, n, i = 1;

double potential(double x){
    E = 1;
    double pot = (E - pow(x, 2));
    return pot;
}

double differential(double (*func)(double x), double x, double h){
    double diff1st = (-func(x + 2*h) + 8*func(x + h) - 8*func(x - h) + func(x - 2*h))/(12*h);
    return diff1st;
}

double differential2nd(double (*func)(double x), double x, double h){
    double diff2nd = (func(x + h) + func(x - h) - 2*func(x))/pow(h, 2);
    return diff2nd;
}


FILE *fp;

double psi[500];

double solve_odd(double h, double x_right){
    fp = fopen( "output.txt", "w" );
    /* boundary condition */
    
    x = x_left;
    psi[0] = 0;
    psi[1] = h + pow(h, 3)/6 * potential(x_left);
    printf( "%e %e\n%e %e\n ", x_left, psi[0], x_left + h, psi[1] );
    fprintf(fp, "%e %e\n%e %e\n ", x_left, psi[0], x_left + h, psi[1] );
    
    /* recursion */
    for (i = 1; i < 500; i++) {
        x = x_left + i*h;
        psi[i+1] =  (2 - 5/6 * pow(h, 2) * potential(x))*psi[i];
        psi[i+1] -= - (1 - pow(h, 2)/12 * potential(x_left))*psi[0];
        psi[i+1] /= 1/(1 - pow(h, 2)/12 * potential(x + h));
        printf( "%e %e\n ", x + h, psi[2]);
        fprintf(fp, "%e %e\n ", x + h, psi[2]);
    }
    fclose( fp );
    return 0;
}

double solve_even(double h, double x_left){
    fp = fopen( "output.txt", "w" );
    /* boundary condition */
    psi[0] = 1;
    psi[1] = psi[0] + pow(h, 2)/2*potential(x_left)*psi[0] + pow(h, 4)/24*(differential2nd(potential, x_left, h)*psi[0] + pow(potential(x_left), 2)*psi[0]);
    printf( "%e %e\n%e %e\n ", x_left, psi[0], x_left + h, psi[1] );
    fprintf(fp, "%e %e\n%e %e\n ", x_left, psi[0], x_left + h, psi[1] );
    /* recursion */
    for (int i = 1; i < 500; i ++) {
        x = x_left + i*h;
        psi[i+1] = (2 + 5/6 * pow(h, 2) * potential(x))*psi[i];
        psi[i+1] -= (1 - pow(h, 2)/12 * potential(x_left))*psi[i-1];
        psi[i+1] /= (1 - pow(h, 2)/12 * potential(x + h));
        printf( "%e %e\n ", x + h, psi[i+1]);
        fprintf(fp, "%e %e\n ", x + h, psi[i+1]);
    }
    fclose( fp );
    return 0;
}
int main(int argc, const char * argv[]) {
    double x_left = 0;
    double h = 0.05;
    
    
    solve_even(h, x_left);
    
    return 0;
}

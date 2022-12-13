#include <iostream>
#include "rezonans.hpp"
#include <math.h>
#include <fstream>
#include <iostream>

Simulation::Simulation() {
    N = 100;
    xmin = 0;
    xmax = 1;
    deltax = (xmax-xmin)/N;
    deltat = 0.0001;
    S = 500000;
    n = 1;
    kappa = 1;
    omega = 3*M_PI*M_PI/2;
    tau = 0;

}

void Simulation::initializeX() {
    for (double x = xmin; x <= xmax; x += deltax )
        xVec.push_back(x);
    xVec.push_back(xmax);


}

void Simulation::initializePsi() {
    for(auto x : xVec) {
        if (x == xmax) 
            psiR.push_back(0);
        else
            psiR.push_back(sqrt(2)*sin(n*M_PI*x));

        psiI.push_back(0);
    }
}

double Simulation::HR(double k, double t) {
    if(k == 0 || k == N)
        return 0;
    else
        return (psiR[k+1]+psiR[k-1]-2.0*psiR[k])/(-2.0*deltax*deltax) + kappa*(xVec[k]-0.5)*psiR[k]*sin(omega*t);
}

double Simulation::HI(double k, double t) {
    if(k == 0 || k == N)
        return 0;
    else
        return (psiI[k+1]+psiI[k-1]-2.0*psiI[k])/(-2.0*deltax*deltax) + kappa*(xVec[k]-0.5)*psiI[k]*sin(omega*t);
}

void Simulation::simulation(char *name) {
    double t = 0;
    std::ofstream file(name);
    for(int step = 0; step < S; step++) {

        for(int k = 0; k < xVec.size(); k++)  
            psiR[k] = psiR[k] + HI(k,t)*deltat/2;

        for(int k = 0; k < xVec.size(); k++)  
            psiI[k] = psiI[k] - HR(k,t + deltat/2.0)*deltat;
                    
        for(int k = 0; k < xVec.size(); k++)  
            psiR[k] = psiR[k] + HI(k,t + deltat)*deltat/2;

        t += deltat;

        if(step%10 == 0){
            norm = 0;
            sr = 0;
            en = 0;
            for(int k = 0; k < xVec.size(); k++) {
                norm += psiR[k]*psiR[k] + psiI[k]*psiI[k];
                sr += xVec[k]*(psiR[k]*psiR[k] + psiI[k]*psiI[k]);
                en += psiR[k]*HR(k,t)+psiI[k]*HI(k,t);
            }
            norm *= deltax;
            sr *= deltax;
            en *= deltax;

            file<<t<<"\t"<<norm<<"\t"<<en<<"\t"<<sr<<"\n";
        }



        

    }

    file.close();

}

int main(int argc, char *argv[]) {
    Simulation *sim = new Simulation();
    sim->initializeX();
    sim->initializePsi();
    sim->simulation(argv[1]);

    return 0;
}
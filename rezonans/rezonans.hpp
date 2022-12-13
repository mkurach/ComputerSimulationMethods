#ifndef _rezonans_hpp_
#define _rezonans_hpp_
#include <vector>

class Simulation {
    public:
        Simulation();
        ~Simulation();
        void initializeX();
        void initializePsi();
        double HR(double k, double t);
        double HI(double k, double t);
        void simulation(char *name);

        int N;
        double xmin,xmax, deltax, deltat, S, n, kappa, omega, tau;
        std::vector<double> xVec, psiR, psiI,tVec;
        double norm, sr, en;

};

#endif
#ifndef _rezonans_hpp_
#define _rezonans_hpp_
#include <vector>

class Simulation {
    public:
        Simulation(double mnoznik);
        ~Simulation();
        void initializeX();
        void initializePsi();
        double HR(double k, double t);
        double HI(double k, double t);
        void simulation(char *name);

        int N;
        double mnoznik;
        double xmin,xmax, deltax, deltat, S, n, kappa, omega0, omega, tau;
        std::vector<double> xVec, psiR, psiI,tVec;
        double norm, sr, en;
        std::vector<double> enVec;

};

//0.9,0.92,0.94,0.95,0.96,0.98,1,1.02,1.04,1.05,1.06,1.08,1.1
#endif
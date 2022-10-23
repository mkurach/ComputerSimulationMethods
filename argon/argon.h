#ifndef _argon_h
#define _argon_h

using namespace std;

struct Parameters{
    double n;
    double m;
    double e;
    double R;
    double f;
    double L;
    double a;
    double T_0;
    double tau;
    int S_o;
    int S_d;
    int S_out;
    int S_xyz;
    
};

struct State{
    double V;
    double P;
    double H;
    double T;
    vector<vector<double>> r;
    vector<vector<double>> p;
    vector<vector<double>> F;
};

void readInput(string fileName, struct Parameters &parameters);
void printState(struct State &state, string fileName, struct Parameters &params);
void initialConditions(struct State &state, struct Parameters &params);

//naglowki funkcji
//struktury

#endif
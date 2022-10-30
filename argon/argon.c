#include <stdio.h>
#include <iostream> 
#include <cmath>
#include <fstream>
#include <time.h>
#include <string.h>
#include <vector>
#include <cmath>
#include <time.h>
#include "argon.h"
using namespace std;

#define kBoltz 8.31e-3

void readInput(string fileName, struct Parameters &parameters) {


    FILE *file = fopen(fileName.data(),"r");
    if (NULL == file) {
        cout<<"Unable to open input file"<<endl;
        exit(1);
    }

    char line[256];
    char smieci[256];
    double tmp;
    int tmpi;

    fgets(line, sizeof(line), file);
    sscanf(line,"%lf %s",&tmp,smieci);
    parameters.n = tmp;
    //cout<<"n :"<<parameters.n<<endl;
    
    fgets(line, sizeof(line), file);
    sscanf(line,"%lf %s",&tmp,smieci);
    parameters.m = tmp;
    //cout<<"m :"<<parameters.m<<endl;

    fgets(line, sizeof(line), file);
    sscanf(line,"%lf %s",&tmp,smieci);
    parameters.e = tmp;
    //cout<<"e :"<<parameters.e<<endl;

    fgets(line, sizeof(line), file);
    sscanf(line,"%lf %s",&tmp,smieci);
    parameters.R = tmp;
    //cout<<"R :"<<parameters.R<<endl;

    fgets(line, sizeof(line), file);
    sscanf(line,"%lf %s",&tmp,smieci);
    parameters.f = tmp;
    //cout<<"f :"<<parameters.f<<endl;

    fgets(line, sizeof(line), file);
    sscanf(line,"%lf %s",&tmp,smieci);
    parameters.L = tmp;
    //cout<<"L :"<<parameters.L<<endl;

    fgets(line, sizeof(line), file);
    sscanf(line,"%lf %s",&tmp,smieci);
    parameters.a = tmp;
    //cout<<"a :"<<parameters.a<<endl;

    fgets(line, sizeof(line), file);
    sscanf(line,"%lf %s",&tmp,smieci);
    parameters.T_0 = tmp;
    //cout<<"T_0 :"<<parameters.T_0<<endl;

    fgets(line, sizeof(line), file);
    sscanf(line,"%lf %s",&tmp,smieci);
    parameters.tau = tmp;
    //cout<<"tau :"<<parameters.tau<<endl;

    fgets(line, sizeof(line), file);
    sscanf(line,"%d %s",&tmpi,smieci);
    parameters.S_o = tmpi;
    //cout<<"S_o :"<<parameters.S_o<<endl;

    fgets(line, sizeof(line), file);
    sscanf(line,"%d %s",&tmpi,smieci);
    parameters.S_d = tmpi;
    //cout<<"S_d :"<<parameters.S_d<<endl;

    fgets(line, sizeof(line), file);
    sscanf(line,"%d %s",&tmpi,smieci);
    parameters.S_out = tmpi;
    //cout<<"S_out :"<<parameters.S_out<<endl;

    fgets(line, sizeof(line), file);
    sscanf(line,"%d %s",&tmpi,smieci);
    parameters.S_xyz = tmpi;
    //cout<<"S_xyz :"<<parameters.S_xyz<<endl;
    
    fclose(file);

    /*ifstream file;
    file.open(fileName.data());

    if(!file) {
        cout<<"Unable to open input file"<<endl;
        exit(1);
    }
    
    file>>parameters.n;
    file>>parameters.m;
    file>>parameters.e;
    file>>parameters.R;
    file>>parameters.f;
    file>>parameters.L;
    file>>parameters.a;
    file>>parameters.T_0;
    file>>parameters.tau;
    file>>parameters.S_o;
    file>>parameters.S_d;
    file>>parameters.S_out;
    file>>parameters.S_xyz;

    file.close();*/
}

void initialConditions(struct State &state, struct Parameters &params){
    
    double b0[3] = {params.a, 0, 0};
    double b1[3] = {params.a/2,  params.a/2*sqrt(3),0};
    double b2[3] = {params.a/2, params.a/6*sqrt(3),params.a*sqrt(2)/sqrt(3)};

    //POLOZENIA

    for(int i = 0; i < params.n; i++){ 
        for(int j = 0; j < params.n; j++) {
            for(int k = 0; k < params.n; k++) { //dla kazdego atomu

                state.r.push_back({(i-(params.n-1)/2)*b0[0]+(j-(params.n-1)/2)*b1[0]+(k-(params.n-1)/2)*b2[0],
                                   (i-(params.n-1)/2)*b0[1]+(j-(params.n-1)/2)*b1[1]+(k-(params.n-1)/2)*b2[1],
                                   (i-(params.n-1)/2)*b0[2]+(j-(params.n-1)/2)*b1[2]+(k-(params.n-1)/2)*b2[2]});
            }
        }
    }

    

    double N = params.n*params.n*params.n;
    int tmp, znak;
    srand(time(NULL));
    state.V = 0;
    state.P = 0;
    state.H = 0;
    state.T = 0;
    double ri, rij, R12, R6;
    double pedCalyx = 0;
    double pedCalyy = 0;
    double pedCalyz = 0;

    for(int i = 0; i < N; i++) { //kazdy atom

        //PEDY  
        tmp = rand();
        if(tmp%2==0)    znak = 1;
        else            znak = -1;

        state.p.push_back({znak*sqrt(2*params.m*(-kBoltz)/2*params.T_0*log((double)(rand()%1000+1)/1000)),
                            znak*sqrt(2*params.m*(-kBoltz)/2*params.T_0*log((double)(rand()%1000+1)/1000)),
                            znak*sqrt(2*params.m*(-kBoltz)/2*params.T_0*log((double)(rand()%1000+1)/1000))});
        
        //pedCaly += sqrt(state.p[i][0]*state.p[i][0] + state.p[i][1]*state.p[i][1] + state.p[i][2]*state.p[i][2]);
        pedCalyx += state.p[i][0];
        pedCalyy += state.p[i][1];
        pedCalyz += state.p[i][2];
        //POTENCJALY I SILY OD SCIANEK
        
        ri = sqrt(state.r[i][0]*state.r[i][0] + state.r[i][1]*state.r[i][1] + state.r[i][2]*state.r[i][2]);
        
        if(ri >= params.L) {
            state.V += params.f/2*(ri-params.L)*(ri-params.L);
            state.F.push_back({params.f*(params.L-ri)/ri*state.r[i][0],
                               params.f*(params.L-ri)/ri*state.r[i][1],
                               params.f*(params.L-ri)/ri*state.r[i][2]});
        }
        else    
            state.F.push_back({0,0,0});
        
        //CISNIENIE
        state.P += sqrt(state.F[i][0]*state.F[i][0] + state.F[i][1]*state.F[i][1] + state.F[i][2]*state.F[i][2])/(4*3.14*params.L*params.L);

        //POTENCJALY I SILY OD INNYCH 
        if(i != 0) {
            for(int j = 0; j < i; j ++){
                rij = sqrt((state.r[i][0]-state.r[j][0])*(state.r[i][0]-state.r[j][0]) + 
                           (state.r[i][1]-state.r[j][1])*(state.r[i][1]-state.r[j][1]) +
                           (state.r[i][2]-state.r[j][2])*(state.r[i][2]-state.r[j][2]));

                R6 = params.R/rij * params.R/rij * params.R/rij * params.R/rij * params.R/rij * params.R/rij;
                R12 = R6 * R6;
                
                state.V += params.e*(R12-2*R6);
                for(int k = 0; k < 3; k++) {
                    state.F[i][k] += 12*params.e*(R12-R6)/rij/rij*(state.r[i][k]-state.r[j][k]);
                    state.F[j][k] += (-12)*params.e*(R12-R6)/rij/rij*(state.r[i][k]-state.r[j][k]);
                }
            }
        }
      
        
    }

    pedCalyx = (double)pedCalyx/N;
    pedCalyy = (double)pedCalyy/N;
    pedCalyz = (double)pedCalyz/N;

    //cout<<pedCalyx<<endl;
    //cout<<pedCalyy<<endl;
    //cout<<pedCalyz<<endl;

    for(int i = 0; i < N; i++) {
        //for(int j = 0; j < 3; j++)
           // state.p[i][j] = state.p[i][j]-pedCaly;
        state.p[i][0] -= pedCalyx;
        state.p[i][1] -= pedCalyy;
        state.p[i][2] -= pedCalyz;
    }



    /*ofstream file;
    file.open("pedy.txt");
    for(auto vec : state.p) {
        file<<vec[0]<<"\t"<<vec[1]<<"\t"<<vec[2]<<"\n";
        //cout<<vec[0]<<"\t"<<vec[1]<<"\t"<<vec[2]<<"\n";
    }
    file.close();*/



}



void simulation(struct State &state, struct Parameters &params, string fileNameMain,string fileNameXYZ) {
    ofstream fileOut, fileOutXYZ;
    fileOut.open(fileNameMain.data());
    fileOutXYZ.open(fileNameXYZ.data());
    
    int N = params.n*params.n*params.n;
    double ri, rij, R12, R6, pi;

    int klatkiCale = params.S_o + params.S_d;
    double Tav = 0;
    double Pav = 0;
    double Hav = 0;

    for(int klatka = 0; klatka < klatkiCale; klatka++) {
        
        state.P = 0;
        state.V = 0;
        state.T = 0;
        state.H = 0;

        
        for(int i = 0; i < N; i++) { //AKTUALIZUJE POLOZENIA I PEDY, ZERUJE SILY
            for(int j = 0; j < 3; j++) { 
                state.p[i][j] += state.F[i][j]/2.0*params.tau;
                state.r[i][j] += state.p[i][j]*params.tau/params.m;
                state.F[i][j] = 0; 
            }
        }
    
        
        for(int i = 0; i < N; i++) {
            //POTENCJALY I SILY OD SCIANEK
            ri = sqrt(state.r[i][0]*state.r[i][0] + state.r[i][1]*state.r[i][1] + state.r[i][2]*state.r[i][2]);
            
            if(ri >= params.L) {
                state.V += params.f/2*(ri-params.L)*(ri-params.L);
                for(int j = 0; j < 3; j ++)
                    state.F[i][j] = params.f*(params.L-ri)/ri*state.r[i][j];
            }
            
            //CISNIENIE
            state.P += sqrt(state.F[i][0]*state.F[i][0] + state.F[i][1]*state.F[i][1] + state.F[i][2]*state.F[i][2])/(4.0*3.14*params.L*params.L);

            
            for(int j = 0; j < i; j ++) { //POTENCJALY I SILY OD INNYCH 

                rij = sqrt((state.r[i][0]-state.r[j][0])*(state.r[i][0]-state.r[j][0]) + 
                        (state.r[i][1]-state.r[j][1])*(state.r[i][1]-state.r[j][1]) +
                        (state.r[i][2]-state.r[j][2])*(state.r[i][2]-state.r[j][2]));

                R6 = params.R/rij * params.R/rij * params.R/rij * params.R/rij * params.R/rij * params.R/rij;
                R12 = R6 * R6;
                state.V += params.e*(R12-2*R6);

                for(int k = 0; k < 3; k++) {
                    state.F[i][k] += 12.0*params.e*(R12-R6)/rij/rij*(state.r[i][k]-state.r[j][k]);
                    state.F[j][k] += (-12.0)*params.e*(R12-R6)/rij/rij*(state.r[i][k]-state.r[j][k]);
                }
            }
            
        }

        for(int i = 0; i < N; i++) {
            //AKTUALIZUJE PEDY 
            for(int j = 0; j < 3; j++)
                state.p[i][j] += state.F[i][j]*params.tau/2;

            //TEMPERATURA I ENERGIA
            pi = sqrt(state.p[i][0]*state.p[i][0] + state.p[i][1]*state.p[i][1] + state.p[i][2]*state.p[i][2]);
            state.T += 2.0/3.0/N/kBoltz/params.m/2.0*pi*pi;
            state.H += pi*pi/2.0/params.m;
        
        }

        state.H += state.V;


        if(klatka%params.S_out == 0) {
            fileOut<<klatka*params.tau<<"\t\t"<<state.H<<"\t"<<state.V<<"\t"<<state.T<<"\t"<<state.P<<"\n";
        }
        if(klatka%params.S_xyz == 0) {
            fileOutXYZ<<N<<"\n\n";
            for(auto vec : state.r)
                fileOutXYZ<<"Ar\t"<<vec[0]<<"\t"<<vec[1]<<"\t"<<vec[2]<<"\n";
            cout<<"Postep: "<<klatka*100/klatkiCale<<"%"<<endl;
    
        }
        if(klatka > params.S_o) {
            Tav += state.T/params.S_d;
            Pav += state.P/params.S_d;
            Hav += state.H/params.S_d;
        }
    } //kazda klatka

    fileOut.close();
    fileOutXYZ.close();




}

int main(int argc, char *argv[]) {
    struct Parameters params;
    struct State state;

    if (argc != 4) {
        cout<<"zla ilosc plikow"<<endl;
        exit(1);
    }
    else {
        readInput(argv[1],params);
        initialConditions(state,params);
        simulation(state, params, argv[2],argv[3]);

    }


    
    
}
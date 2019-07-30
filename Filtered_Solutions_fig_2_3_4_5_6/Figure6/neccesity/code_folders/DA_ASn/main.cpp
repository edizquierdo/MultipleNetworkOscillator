// =============================================================
// Evolution of Integrated Neuromechanical Forward Locomotion
// Eduardo Izquierdo
// Indiana University
// February, 2015
// =============================================================
#include <iostream>
#include <iomanip>  // cout precision
#include <math.h>
#include "TSearch.h"
#include "VectorMatrix.h"
#include "Worm.h"

int skip_steps = 10;
using namespace std;

// Integration parameters
const double Duration = 20.0;       //
const double Transient = 50.0;       //
const double StepSize = 0.005;
const int N_curvs = 23;             // Number of cuvature points

double OSCT = 0.25 * Duration; // Cap for oscillation evaluation
const double agarfreq = 0.44;

// Genotype -> Phenotype Mapping (Ventral cord)
const double	BiasRange				= 15.0;
const double    SCRange                 = 15.0;
const double    CSRange                 = 15.0;
const double    TauMin                 = 0.1;
const double    TauMax                 = 2.5;
const double    ESRange                 = 2.0;
const double    NMJmax                  = 1.2;
const double    IIRange                 = 15.0;

// Fitness
const double    AvgSpeed = 0.00022;             // Average speed of the worm in meters per seconds
const double    BBCfit = AvgSpeed*Duration;

// Size of genotype (VC)
int	VectSize = 44;

// ------------------------------------
// Genotype-Phenotype Mapping
// ------------------------------------
void GenPhenMapping(TVector<double> &gen, TVector<double> &phen)
{
    // Bias
    for (int i = 1; i <= 7; i++){
        phen(i) = MapSearchParameter(gen(i), -BiasRange, BiasRange);
    }
    // Time Constant
    for (int i = 8; i <= 14; i++){
        phen(i) = MapSearchParameter(gen(i), TauMin, TauMax);
    }
    // Self connections
    for (int i = 15; i <= 21; i++){
        phen(i) = MapSearchParameter(gen(i), -SCRange, SCRange);
    }
    // Chemical synapses
    for (int i = 22; i <=30; i++){
        phen(i) = MapSearchParameter(gen(i), -CSRange, CSRange);
    }
    
    // Gap junctions
    phen(31) = MapSearchParameter(gen(31), 0.0, ESRange);

    // Intersegment synapse tested
    phen(40) = MapSearchParameter(gen(40), -CSRange, CSRange);  // DB to DDnext
    phen(41) = MapSearchParameter(gen(41), -CSRange, CSRange);  // VAnext to DD
    phen(42) = MapSearchParameter(gen(42), 0.0, ESRange);       // AS -- VAnext
    phen(43) = MapSearchParameter(gen(43), 0.0, ESRange);       // DA -- ASnext
    phen(44) = MapSearchParameter(gen(44), 0.0, ESRange);       // VB -- DBnext

    // NMJ Weight
    phen(32) = MapSearchParameter(gen(32), 0.0, NMJmax);       // AS
    phen(33) = MapSearchParameter(gen(33), 0.0, NMJmax);       // DA
    phen(34) = MapSearchParameter(gen(34), NMJmax, NMJmax);       // DB
    phen(35) = MapSearchParameter(gen(35), -NMJmax, 0.0);      // DD
    phen(36) = MapSearchParameter(gen(36), -NMJmax, 0.0);      // VD
    phen(37) = MapSearchParameter(gen(37), NMJmax, NMJmax);      // VB
    phen(38) = MapSearchParameter(gen(38), 0.0, NMJmax);      // VA
    
    phen(39) = MapSearchParameter(gen(39), 0.2, 1.0);       // Used to be 0.4/0.6 XXX NMJ_Gain Mapping
}
void curvRatio(TVector<double> &v, TVector<double> &antposcurv){
    for (int i = 1; i <= N_curvs; i++){
        if (i <= 11)
            antposcurv(1) += fabs(v(i));
        else
            antposcurv(2) += fabs(v(i));
    }
}
//////// evaluation ////////////
/////////////////////////////
double EvaluationFunction(TVector<double> &v, RandomState &rs, int &ConstantInput){
    // Fitness
    double bodyorientation, anglediff;
    double movementorientation, distancetravelled = 0, displacement, temp;
    TVector<double> curvature(1, N_curvs);
    TVector<double> antpostcurv(1, 2);
    antpostcurv.FillContents(0.0);

    // Genotype-Phenotype Mapping
    TVector<double> phenotype(1, VectSize);
    GenPhenMapping(v, phenotype);
    
    // Intersegment synapse tested
    double SynWeight = phenotype(43);
    phenotype(43) = 0.0; // removed DA -- AS

    Worm w(phenotype);
    w.InitializeState(rs);
    
    // Transient XXX
    w.SetAVB(0.0);
    w.SetAVA(0.0);


    for (int seg = 0; seg < 6; seg++){ // set opposite external input to DA and ASnext
        w.n.SetNeuronExternalInput(DA + seg*7,      SynWeight * 0.75*(ConstantInput%21- 10));
        w.n.SetNeuronExternalInput(AS + (seg+1)*7,  SynWeight * 0.75*(ConstantInput/21- 10));
    }

    for (double t = 0.0; t < 20; t += StepSize){    //load input to neurons in segment 3
        w.Step(StepSize);
    }
    
    double xt = w.CoMx(), xtp;
    double yt = w.CoMy(), ytp;
    
    for (double t = 0.0; t <= Duration; t += StepSize){    //load input to neurons in segment 3
        w.Step(StepSize);
        //// Locomotion
        // Current and past centroid position
        xtp = xt; ytp = yt;
        xt = w.CoMx(); yt = w.CoMy();
        
        // Integration error check
        if (isnan(xt) || isnan(yt) || sqrt(pow(xt-xtp,2)+pow(yt-ytp,2)) > 100*AvgSpeed*StepSize){
            return 0.0;
        }
        
        // Fitness
        bodyorientation = w.Orientation();                  // Orientation of the body position
        movementorientation = atan2(yt-ytp,xt-xtp);         // Orientation of the movement
        anglediff = movementorientation - bodyorientation;  // Check how orientations align
        temp = cos(anglediff) > 0.0 ? 1.0 : -1.0;           // Add to fitness only movement forward
        distancetravelled += temp * sqrt(pow(xt-xtp,2)+pow(yt-ytp,2));
    }
    
    return (1 - (fabs(BBCfit-distancetravelled)/BBCfit));
}

// ------------------------------------
// The main program
// ------------------------------------
int main (int argc, const char* argv[]){
    // ConstantInput [0, 200] as an argument from Python script 
    int ConstantInput = atoi(argv[1]);

    RandomState rs;
    long seed = static_cast<long>(time(NULL));
    rs.SetRandomSeed(seed);
    
    ofstream fitness;
    fitness.open("fitness.dat");
    
    // Code to run simulation:
    InitializeBodyConstants();
    
    ifstream BestIndividualFile;
    TVector<double> bestVector(1, VectSize);
    BestIndividualFile.open("./best.gen.dat");
    BestIndividualFile >> bestVector;
    for (int rep= 1; rep<=8; rep++){
        fitness << EvaluationFunction(bestVector, rs, ConstantInput) << " ";
    }
    fitness.close();
    return 0;
}

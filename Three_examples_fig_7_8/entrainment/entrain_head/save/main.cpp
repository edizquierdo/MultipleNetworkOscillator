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
#include <string>
#include <sstream>

int skip_steps = 1;
int ind;
using namespace std;

// Integration parameters
const double Duration = 40.0;       //
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

//////// evaluation ////////////
/////////////////////////////
double EvaluationFunction2(TVector<double> &v, RandomState &rs){
    ofstream actfile(("act_" + to_string(ind) + ".dat").c_str());

    // Genotype-Phenotype Mapping
    TVector<double> phenotype(1, VectSize);
    GenPhenMapping(v, phenotype);
    
    Worm w(phenotype);
    w.InitializeState(rs);
    w.n.RandomizeCircuitOutput(0.5, 0.5, rs);
    w.SetAVB(0.0);
    w.SetAVA(0.0);

    for (double t = 0.0; t < 100; t += StepSize) {
        // Step simulation
        w.Step(StepSize);
    }
    // Time loop
    for (double t = 0.0; t <= 250; t += StepSize) {
        for (int n = 1; n<=7; n++){actfile << w.n.NeuronState(n) << endl;} // save neural traces for first segment
        // Step simulation
        w.Step(StepSize);
    }
    actfile.close();
    return 0.0;
}

// ------------------------------------
// The main program
// ------------------------------------
int main (int argc, const char* argv[]){
    RandomState rs;
    long seed = static_cast<long>(time(NULL));
    rs.SetRandomSeed(seed);
    if (argc > 1){
        ind = atoi(argv[1]);
    }
    cout << ("act_" + to_string(ind) + ".dat").c_str() << endl;
    // Code to run simulation:
    InitializeBodyConstants();
    
    ifstream BestIndividualFile;
    TVector<double> bestVector(1, VectSize);
    BestIndividualFile.open("./best.gen.dat");
    BestIndividualFile >> bestVector;
    EvaluationFunction2(bestVector, rs);
    return 0;
}

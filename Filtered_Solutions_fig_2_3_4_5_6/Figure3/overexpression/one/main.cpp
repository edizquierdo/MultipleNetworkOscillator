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
#include "Segment.h"

//#define EVOLVE
//#define PRINTTOFILE
#define OUTPUT
int skip_steps = 4;
double GJFactor = 1;
using namespace std;

// Integration parameters
const double Duration = 40.0;       //
const double Transient = 10.0;       //
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
    // 5 times Augmented jup junction (UNC-9 overexpresion as un Xu et al, 2018)
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
//////// Stage 1 ////////////
double EvaluationFunction1(TVector<double> &v, RandomState &rs){
    return 0;
}
//////// Stage 2 ////////////
/////////////////////////////
double EvaluationFunction2(TVector<double> &v, RandomState &rs){
    ofstream bodyfile;
    bodyfile.open("body.dat");
    // Genotype-Phenotype Mapping
    TVector<double> phenotype(1, VectSize);
    GenPhenMapping(v, phenotype);
    Worm w(phenotype);
    double wvbdb = phenotype(44);
    w.InitializeState(rs);
    
    GJFactor = 0.0;
    for (int s = 0; s < 6; s++){w.n.SetElectricalSynapseWeight(VB+7*s, DB+7*(s+1), GJFactor + wvbdb);}
    for (int s = 0; s < 6; s++){w.n.SetElectricalSynapseWeight(VB+7*s, VB+7*(s+1), GJFactor);}
    for (int s = 0; s < 6; s++){w.n.SetElectricalSynapseWeight(DB+7*s, DB+7*(s+1), GJFactor);}
    for (double t = 0.0; t <= Transient; t += StepSize){
        w.Step(StepSize);
    }
    for (double t = 0.0; t <= 30; t += StepSize){
        w.Step(StepSize);
        w.DumpBodyState(bodyfile, skip_steps);
    }

    GJFactor = 0.25;
    for (int s = 0; s < 6; s++){w.n.SetElectricalSynapseWeight(VB+7*s, DB+7*(s+1), GJFactor + wvbdb);}
    for (int s = 0; s < 6; s++){w.n.SetElectricalSynapseWeight(VB+7*s, VB+7*(s+1), GJFactor);}
    for (int s = 0; s < 6; s++){w.n.SetElectricalSynapseWeight(DB+7*s, DB+7*(s+1), GJFactor);}
    for (double t = 0.0; t <= Transient; t += StepSize){
        w.Step(StepSize);
    }
    for (double t = 0.0; t <= 30; t += StepSize){
        w.Step(StepSize);
        w.DumpBodyState(bodyfile, skip_steps);
    }

    GJFactor = 0.5;
    for (int s = 0; s < 6; s++){w.n.SetElectricalSynapseWeight(VB+7*s, DB+7*(s+1), GJFactor + wvbdb);}
    for (int s = 0; s < 6; s++){w.n.SetElectricalSynapseWeight(VB+7*s, VB+7*(s+1), GJFactor);}
    for (int s = 0; s < 6; s++){w.n.SetElectricalSynapseWeight(DB+7*s, DB+7*(s+1), GJFactor);}
    for (double t = 0.0; t <= Transient; t += StepSize){
        w.Step(StepSize);
    }
    for (double t = 0.0; t <= 30; t += StepSize){
        w.Step(StepSize);
        w.DumpBodyState(bodyfile, skip_steps);
    }

    GJFactor = 1.0;
    for (int s = 0; s < 6; s++){w.n.SetElectricalSynapseWeight(VB+7*s, DB+7*(s+1), GJFactor + wvbdb);}
    for (int s = 0; s < 6; s++){w.n.SetElectricalSynapseWeight(VB+7*s, VB+7*(s+1), GJFactor);}
    for (int s = 0; s < 6; s++){w.n.SetElectricalSynapseWeight(DB+7*s, DB+7*(s+1), GJFactor);}
    for (double t = 0.0; t <= Transient; t += StepSize){
        w.Step(StepSize);
    }
    for (double t = 0.0; t <= 30; t += StepSize){
        w.Step(StepSize);
        w.DumpBodyState(bodyfile, skip_steps);
    }

    GJFactor = 2.0;
    for (int s = 0; s < 6; s++){w.n.SetElectricalSynapseWeight(VB+7*s, DB+7*(s+1), GJFactor + wvbdb);}
    for (int s = 0; s < 6; s++){w.n.SetElectricalSynapseWeight(VB+7*s, VB+7*(s+1), GJFactor);}
    for (int s = 0; s < 6; s++){w.n.SetElectricalSynapseWeight(DB+7*s, DB+7*(s+1), GJFactor);}
    for (double t = 0.0; t <= Transient; t += StepSize){
        w.Step(StepSize);
    }
    for (double t = 0.0; t <= 30; t += StepSize){
        w.Step(StepSize);
        w.DumpBodyState(bodyfile, skip_steps);
    }

    GJFactor = 4.0;
    for (int s = 0; s < 6; s++){w.n.SetElectricalSynapseWeight(VB+7*s, DB+7*(s+1), GJFactor + wvbdb);}
    for (int s = 0; s < 6; s++){w.n.SetElectricalSynapseWeight(VB+7*s, VB+7*(s+1), GJFactor);}
    for (int s = 0; s < 6; s++){w.n.SetElectricalSynapseWeight(DB+7*s, DB+7*(s+1), GJFactor);}
    for (double t = 0.0; t <= Transient; t += StepSize){
        w.Step(StepSize);
    }
    for (double t = 0.0; t <= 30; t += StepSize){
        w.Step(StepSize);
        w.DumpBodyState(bodyfile, skip_steps);
    }

    GJFactor = 8.0;
    for (int s = 0; s < 6; s++){w.n.SetElectricalSynapseWeight(VB+7*s, DB+7*(s+1), GJFactor + wvbdb);}
    for (int s = 0; s < 6; s++){w.n.SetElectricalSynapseWeight(VB+7*s, VB+7*(s+1), GJFactor);}
    for (int s = 0; s < 6; s++){w.n.SetElectricalSynapseWeight(DB+7*s, DB+7*(s+1), GJFactor);}
    for (double t = 0.0; t <= Transient; t += StepSize){
        w.Step(StepSize);
    }
    for (double t = 0.0; t <= 30; t += StepSize){
        w.Step(StepSize);
        w.DumpBodyState(bodyfile, skip_steps);
    }

    GJFactor = 16.0;
    for (int s = 0; s < 6; s++){w.n.SetElectricalSynapseWeight(VB+7*s, DB+7*(s+1), GJFactor + wvbdb);}
    for (int s = 0; s < 6; s++){w.n.SetElectricalSynapseWeight(VB+7*s, VB+7*(s+1), GJFactor);}
    for (int s = 0; s < 6; s++){w.n.SetElectricalSynapseWeight(DB+7*s, DB+7*(s+1), GJFactor);}
    for (double t = 0.0; t <= Transient; t += StepSize){
        w.Step(StepSize);
    }
    for (double t = 0.0; t <= 30; t += StepSize){
        w.Step(StepSize);
        w.DumpBodyState(bodyfile, skip_steps);
    }


    bodyfile.close();

    return 0.0;
}

// ------------------------------------
// Display functions
// ------------------------------------
void EvolutionaryRunDisplay(int Generation, double BestPerf, double AvgPerf, double PerfVar){
    cout << Generation << " " << BestPerf << " " << AvgPerf << " " << PerfVar << endl;
}

void ResultsDisplay(TSearch &s){
    TVector<double> bestVector;
    ofstream BestIndividualFile;
    
    bestVector = s.BestIndividual();
    BestIndividualFile.open("best.gen.dat");
    BestIndividualFile << setprecision(32);
    BestIndividualFile << bestVector << endl;
    BestIndividualFile.close();
}

int finish_Bosc(int Generation,double BestPerf,double AvgPerf,double PerfVar){
    if (BestPerf > 0.99) return 1;
    else return 0;
}

// ------------------------------------
// The main program
// ------------------------------------
#ifdef EVOLVE
int main (int argc, const char* argv[]){
    std::cout << std::setprecision(10);
    
    long randomseed = static_cast<long>(time(NULL));
    
    if (argc == 2)
        randomseed += atoi(argv[1]);
    
    TSearch s(VectSize);
    TVector<double> phenotype(1, VectSize);
    
    // save the seed to a file
    ofstream seedfile;
    seedfile.open ("seed.dat");
    seedfile << randomseed << endl;
    seedfile.close();
    
    // configure the search
    s.SetRandomSeed(randomseed);
    s.SetPopulationStatisticsDisplayFunction(EvolutionaryRunDisplay);
    s.SetSearchResultsDisplayFunction(ResultsDisplay);
    s.SetSelectionMode(RANK_BASED);             //{FITNESS_PROPORTIONATE,RANK_BASED}
    s.SetReproductionMode(GENETIC_ALGORITHM);	// {HILL_CLIMBING, GENETIC_ALGORITHM}
    s.SetPopulationSize(100);
    s.SetMaxGenerations(2000);
    s.SetMutationVariance(0.1);
    s.SetCrossoverProbability(0.5);
    s.SetCrossoverMode(UNIFORM);              //{UNIFORM, TWO_POINT}
    s.SetMaxExpectedOffspring(1.1);
    s.SetElitistFraction(0.04);
    s.SetSearchConstraint(1);
    s.SetCheckpointInterval(0);
    s.SetReEvaluationFlag(0);

    // Stage 1 //
    s.SetSearchTerminationFunction(finish_Bosc);
    s.SetEvaluationFunction(EvaluationFunction1);
    // redirect standard output to a file
#ifdef PRINTTOFILE
    ofstream evolfile;
    evolfile.open ("fitness.dat");
    cout.rdbuf(evolfile.rdbuf());
#endif
    s.ExecuteSearch();

    // Stage 2 //
    s.SetSearchTerminationFunction(NULL);
    s.SetEvaluationFunction(EvaluationFunction2);
    InitializeBodyConstants();
    s.ExecuteSearch();

#ifdef PRINTTOFILE
    evolfile.close();
#endif
    return 0;
}
#else
int main (int argc, const char* argv[]){
    RandomState rs;
    long seed = static_cast<long>(time(NULL));
    rs.SetRandomSeed(seed);
    
    std::cout << std::setprecision(10);
    
    // Code to run simulation:
    InitializeBodyConstants();
    
    ifstream BestIndividualFile;
    TVector<double> bestVector(1, VectSize);
    BestIndividualFile.open("best.gen.dat");
    BestIndividualFile >> bestVector;
    EvaluationFunction2(bestVector, rs);
    return 0;
}
#endif

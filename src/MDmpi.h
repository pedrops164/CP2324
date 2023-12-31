
/* CONSTANTS */
//  Lennard-Jones parameters in natural units!
#define root_rank 0 // defines the rank of the root
int sigma = 1;
int sigma2 = 1; //sigma2 = sigma1 * sigma1
int epsilon = 1;
int epsilon8 = 8;
int m = 1;
int kB = 1;
double NA = 6.022140857e23;
double kBSI = 1.38064852e-23;  // m^2*kg/(s^2*K)
//  Size of box, which will be specified in natural units
double L;
//  Initial Temperature in Natural Units
double Tinit;  //2;
// atom type
char atype[10];

/* MPI processes local variables */
int start_index=-1;
int end_index=-1;

//void final_prints() {
//    fprintf(afp,"  Total Time (s)      T (K)               P (Pa)      PV/nT (J/(mol K))         Z           V (m^3)              N\n");
//    fprintf(afp," --------------   -----------        ---------------   --------------   ---------------   ------------   -----------\n");
//    fprintf(afp,"  %8.4e  %15.5f       %15.5f     %10.5f       %10.5f        %10.5e         %i\n",i*dt*timefac,Tavg,Pavg,gc,Z,Vol*VolFac,N);
//    
//    printf("\n  TO ANIMATE YOUR SIMULATION, OPEN THE FILE \n  '%s' WITH VMD AFTER THE SIMULATION COMPLETES\n",tfn);
//    printf("\n  TO ANALYZE INSTANTANEOUS DATA ABOUT YOUR MOLECULE, OPEN THE FILE \n  '%s' WITH YOUR FAVORITE TEXT EDITOR OR IMPORT THE DATA INTO EXCEL\n",ofn);
//    printf("\n  THE FOLLOWING THERMODYNAMIC AVERAGES WILL BE COMPUTED AND WRITTEN TO THE FILE  \n  '%s':\n",afn);
//    printf("\n  AVERAGE TEMPERATURE (K):                 %15.5f\n",Tavg);
//    printf("\n  AVERAGE PRESSURE  (Pa):                  %15.5f\n",Pavg);
//    printf("\n  PV/nT (J * mol^-1 K^-1):                 %15.5f\n",gc);
//    printf("\n  PERCENT ERROR of pV/nT AND GAS CONSTANT: %15.5f\n",100*fabs(gc-8.3144598)/8.3144598);
//    printf("\n  THE COMPRESSIBILITY (unitless):          %15.5f \n",Z);
//    printf("\n  TOTAL VOLUME (m^3):                      %10.5e \n",Vol*VolFac);
//    printf("\n  NUMBER OF PARTICLES (unitless):          %i \n", N);
//}

/* Data structures declarations */
//  Position data structure
Position position;
//  Velocity data structure
Velocity velocity;
//  Acceleration data structure
Acceleration acceleration;

/* Function headers */
void initializeVariables();
//  Function prototypes
//  initialize positions on simple cubic lattice, also calls function to initialize velocities
void initialize();  
//  update positions and velocities using Velocity Verlet algorithm 
//  print particle coordinates to file for rendering via VMD or other animation software
//  return 'instantaneous pressure'
double VelocityVerlet(double dt, int iter, FILE *fp, double * PE);  
//  Compute Force using F = -dV/dr
//  solve F = ma for use in Velocity Verlet, and compute Potential energy aswell
double computeAccelerationsAndPotential();
//  Numerical Recipes function for generation gaussian distribution
double gaussdist();
//  Initialize velocities according to user-supplied initial Temperature (Tinit)
void initializeVelocities();
//  Compute mean squared velocity from particle velocities
double MeanSquaredVelocity();
//  Compute total kinetic energy from particle mass and velocities
double Kinetic();
// efficient recursive function that implements exponentiation by squaring
double pow_by_squaring(double x, int n);
// potency calculation between two particles, given the inverse of the square of the distance between the particles
double calculateLocalPotential(double inverseSquaredDistance);
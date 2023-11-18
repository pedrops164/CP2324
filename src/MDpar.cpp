/*
 MD.c - a simple molecular dynamics program for simulating real gas properties of Lennard-Jones particles.
 
 Copyright (C) 2016  Jonathan J. Foley IV, Chelsea Sweet, Oyewumi Akinfenwa
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
 Electronic Contact:  foleyj10@wpunj.edu
 Mail Contact:   Prof. Jonathan Foley
 Department of Chemistry, William Paterson University
 300 Pompton Road
 Wayne NJ 07470
 
 */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<cmath>
#include<string.h>

// Number of particles
//int N;
const int N = 5000;

// Number of threads
const int nworkers = 4;

//  Lennard-Jones parameters in natural units!
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

// Max number of particles
const int MAXPART=5001;

/*
Here we define the Position, Velocity, Acceleration and Force matrices as 3 seperate arrays instead.
We use Structure of Arrays instead of Array of Structures.
With this implementation, since a cache line can hold 8 doubles, every time an element of the array is accessed in the memory, that element and the next 7 doubles
are copied to the cache. This will be helpful because this program demonstrates temporal and spatial locality, so this aspect is very useful
*/
//  Position
double r_x[N];
double r_y[N];
double r_z[N];
//  Velocity
double v_x[N];
double v_y[N];
double v_z[N];
//  Acceleration
double a_x[N];
double a_y[N];
double a_z[N];
//  Force
double F_x[N];
double F_y[N];
double F_z[N];

// atom type
char atype[10];
//  Function prototypes
//  initialize positions on simple cubic lattice, also calls function to initialize velocities
void initialize();  
//  update positions and velocities using Velocity Verlet algorithm 
//  print particle coordinates to file for rendering via VMD or other animation software
//  return 'instantaneous pressure'
double VelocityVerlet(double dt, int iter, FILE *fp);  
//  Compute Force using F = -dV/dr
//  solve F = ma for use in Velocity Verlet
void computeAccelerations();
//  Numerical Recipes function for generation gaussian distribution
double gaussdist();
//  Initialize velocities according to user-supplied initial Temperature (Tinit)
void initializeVelocities();
//  Compute total potential energy from particle coordinates
double Potential();
//  Compute mean squared velocity from particle velocities
double MeanSquaredVelocity();
//  Compute total kinetic energy from particle mass and velocities
double Kinetic();
// efficient recursive function that implements exponentiation by squaring
double pow_by_squaring(double x, int n);

int main()
{
    
    //  variable delcarations
    int i;
    double dt, Vol, Temp, Press, Pavg, Tavg, rho;
    double VolFac, TempFac, PressFac, timefac;
    double KE, PE, mvs;
    char prefix[1000], tfn[1000], ofn[1000], afn[1000];
    FILE *infp, *tfp, *ofp, *afp;
    
    
    printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("                  WELCOME TO WILLY P CHEM MD!\n");
    printf("  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("\n  ENTER A TITLE FOR YOUR CALCULATION!\n");
    scanf("%s",prefix);
    strcpy(tfn,prefix);
    strcat(tfn,"_traj.xyz");
    strcpy(ofn,prefix);
    strcat(ofn,"_output.txt");
    strcpy(afn,prefix);
    strcat(afn,"_average.txt");
    
    printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("                  TITLE ENTERED AS '%s'\n",prefix);
    printf("  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    
    /*     Table of values for Argon relating natural units to SI units:
     *     These are derived from Lennard-Jones parameters from the article
     *     "Liquid argon: Monte carlo and molecular dynamics calculations"
     *     J.A. Barker , R.A. Fisher & R.O. Watts
     *     Mol. Phys., Vol. 21, 657-673 (1971)
     *
     *     mass:     6.633e-26 kg          = one natural unit of mass for argon, by definition
     *     energy:   1.96183e-21 J      = one natural unit of energy for argon, directly from L-J parameters
     *     length:   3.3605e-10  m         = one natural unit of length for argon, directly from L-J parameters
     *     volume:   3.79499-29 m^3        = one natural unit of volume for argon, by length^3
     *     time:     1.951e-12 s           = one natural unit of time for argon, by length*sqrt(mass/energy)
     ***************************************************************************************/
    
    //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //  Edit these factors to be computed in terms of basic properties in natural units of
    //  the gas being simulated
    
    
    printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("  WHICH NOBLE GAS WOULD YOU LIKE TO SIMULATE? (DEFAULT IS ARGON)\n");
    printf("\n  FOR HELIUM,  TYPE 'He' THEN PRESS 'return' TO CONTINUE\n");
    printf("  FOR NEON,    TYPE 'Ne' THEN PRESS 'return' TO CONTINUE\n");
    printf("  FOR ARGON,   TYPE 'Ar' THEN PRESS 'return' TO CONTINUE\n");
    printf("  FOR KRYPTON, TYPE 'Kr' THEN PRESS 'return' TO CONTINUE\n");
    printf("  FOR XENON,   TYPE 'Xe' THEN PRESS 'return' TO CONTINUE\n");
    printf("  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    scanf("%s",atype);
    
    if (strcmp(atype,"He")==0) {
        
        VolFac = 1.8399744000000005e-29;
        PressFac = 8152287.336171632;
        TempFac = 10.864459551225972;
        timefac = 1.7572698825166272e-12;
        
    }
    else if (strcmp(atype,"Ne")==0) {
        
        VolFac = 2.0570823999999997e-29;
        PressFac = 27223022.27659913;
        TempFac = 40.560648991243625;
        timefac = 2.1192341945685407e-12;
        
    }
    else if (strcmp(atype,"Ar")==0) {
        
        VolFac = 3.7949992920124995e-29;
        PressFac = 51695201.06691862;
        TempFac = 142.0950000000000;
        timefac = 2.09618e-12;
        //strcpy(atype,"Ar");
        
    }
    else if (strcmp(atype,"Kr")==0) {
        
        VolFac = 4.5882712000000004e-29;
        PressFac = 59935428.40275003;
        TempFac = 199.1817584391428;
        timefac = 8.051563913585078e-13;
        
    }
    else if (strcmp(atype,"Xe")==0) {
        
        VolFac = 5.4872e-29;
        PressFac = 70527773.72794868;
        TempFac = 280.30305642163006;
        timefac = 9.018957925790732e-13;
        
    }
    else {
        
        VolFac = 3.7949992920124995e-29;
        PressFac = 51695201.06691862;
        TempFac = 142.0950000000000;
        timefac = 2.09618e-12;
        strcpy(atype,"Ar");
        
    }
    printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("\n                     YOU ARE SIMULATING %s GAS! \n",atype);
    printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    
    printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("\n  YOU WILL NOW ENTER A FEW SIMULATION PARAMETERS\n");
    printf("  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("\n\n  ENTER THE INTIAL TEMPERATURE OF YOUR GAS IN KELVIN\n");
    scanf("%lf",&Tinit);
    // Make sure temperature is a positive number!
    if (Tinit<0.) {
        printf("\n  !!!!! ABSOLUTE TEMPERATURE MUST BE A POSITIVE NUMBER!  PLEASE TRY AGAIN WITH A POSITIVE TEMPERATURE!!!\n");
        exit(0);
    }
    // Convert initial temperature from kelvin to natural units
    Tinit /= TempFac;
    
    
    printf("\n\n  ENTER THE NUMBER DENSITY IN moles/m^3\n");
    printf("  FOR REFERENCE, NUMBER DENSITY OF AN IDEAL GAS AT STP IS ABOUT 40 moles/m^3\n");
    printf("  NUMBER DENSITY OF LIQUID ARGON AT 1 ATM AND 87 K IS ABOUT 35000 moles/m^3\n");
    
    scanf("%lf",&rho);
    
    //N = 10*216;
    Vol = N/(rho*NA);
    
    Vol /= VolFac;
    
    //  Limiting N to MAXPART for practical reasons
    if (N>=MAXPART) {
        
        printf("\n\n\n  MAXIMUM NUMBER OF PARTICLES IS %i\n\n  PLEASE ADJUST YOUR INPUT FILE ACCORDINGLY \n\n", MAXPART);
        exit(0);
        
    }
    //  Check to see if the volume makes sense - is it too small?
    //  Remember VDW radius of the particles is 1 natural unit of length
    //  and volume = L*L*L, so if V = N*L*L*L = N, then all the particles
    //  will be initialized with an interparticle separation equal to 2xVDW radius
    if (Vol<N) {
        
        printf("\n\n\n  YOUR DENSITY IS VERY HIGH!\n\n");
        printf("  THE NUMBER OF PARTICLES IS %i AND THE AVAILABLE VOLUME IS %f NATURAL UNITS\n",N,Vol);
        printf("  SIMULATIONS WITH DENSITY GREATER THAN 1 PARTCICLE/(1 Natural Unit of Volume) MAY DIVERGE\n");
        printf("  PLEASE ADJUST YOUR INPUT FILE ACCORDINGLY AND RETRY\n\n");
        exit(0);
    }
    // Vol = L*L*L;
    // Length of the box in natural units:
    L = pow(Vol,(1./3));
    
    //  Files that we can write different quantities to
    tfp = fopen(tfn,"w");     //  The MD trajectory, coordinates of every particle at each timestep
    ofp = fopen(ofn,"w");     //  Output of other quantities (T, P, gc, etc) at every timestep
    afp = fopen(afn,"w");    //  Average T, P, gc, etc from the simulation
    
    int NumTime;
    if (strcmp(atype,"He")==0) {
        
        // dt in natural units of time s.t. in SI it is 5 f.s. for all other gasses
        dt = 0.2e-14/timefac;
        //  We will run the simulation for NumTime timesteps.
        //  The total time will be NumTime*dt in natural units
        //  And NumTime*dt multiplied by the appropriate conversion factor for time in seconds
        NumTime=50000;
    }
    else {
        dt = 0.5e-14/timefac;
        NumTime=200;
        
    }
    
    //  Put all the atoms in simple crystal lattice and give them random velocities
    //  that corresponds to the initial temperature we have specified
    initialize();
    
    //  Based on their positions, calculate the ininial intermolecular forces
    //  The accellerations of each particle will be defined from the forces and their
    //  mass, and this will allow us to update their positions via Newton's law
    computeAccelerations();
    
    
    // Print number of particles to the trajectory file
    fprintf(tfp,"%i\n",N);
    
    //  We want to calculate the average Temperature and Pressure for the simulation
    //  The variables need to be set to zero initially
    Pavg = 0;
    Tavg = 0;
    
    
    int tenp = floor(NumTime/10);
    fprintf(ofp,"  time (s)              T(t) (K)              P(t) (Pa)           Kinetic En. (n.u.)     Potential En. (n.u.) Total En. (n.u.)\n");
    printf("  PERCENTAGE OF CALCULATION COMPLETE:\n  [");
    //#pragma omp parallel num_threads(nworkers)
    for (i=0; i<NumTime+1; i++) {
        
        //  This just prints updates on progress of the calculation for the users convenience
        if (i==tenp) printf(" 10 |");
        else if (i==2*tenp) printf(" 20 |");
        else if (i==3*tenp) printf(" 30 |");
        else if (i==4*tenp) printf(" 40 |");
        else if (i==5*tenp) printf(" 50 |");
        else if (i==6*tenp) printf(" 60 |");
        else if (i==7*tenp) printf(" 70 |");
        else if (i==8*tenp) printf(" 80 |");
        else if (i==9*tenp) printf(" 90 |");
        else if (i==10*tenp) printf(" 100 ]\n");
        fflush(stdout);
        
        
        // This updates the positions and velocities using Newton's Laws
        // Also computes the Pressure as the sum of momentum changes from wall collisions / timestep
        // which is a Kinetic Theory of gasses concept of Pressure
        Press = VelocityVerlet(dt, i+1, tfp);
        Press *= PressFac;
        
        //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //  Now we would like to calculate somethings about the system:
        //  Instantaneous mean velocity squared, Temperature, Pressure
        //  Potential, and Kinetic Energy
        //  We would also like to use the IGL to try to see if we can extract the gas constant
        mvs = MeanSquaredVelocity();
        KE = Kinetic();
        PE = Potential();
        
        // Temperature from Kinetic Theory
        Temp = m*mvs/(3*kB) * TempFac;
        
        // Instantaneous gas constant and compressibility - not well defined because
        // pressure may be zero in some instances because there will be zero wall collisions,
        // pressure may be very high in some instances because there will be a number of collisions
        //gc = NA*Press*(Vol*VolFac)/(N*Temp);
        //Z  = Press*(Vol*VolFac)/(N*kBSI*Temp);
        
        Tavg += Temp;
        Pavg += Press;
        
        fprintf(ofp,"  %8.4e  %20.8f  %20.8f %20.8f  %20.8f  %20.8f \n",i*dt*timefac,Temp,Press,KE, PE, KE+PE);
        
        
    }
    
    // Because we have calculated the instantaneous temperature and pressure,
    // we can take the average over the whole simulation here
    Pavg /= NumTime;
    Tavg /= NumTime;
    double Z = Pavg*(Vol*VolFac)/(N*kBSI*Tavg);
    double gc = NA*Pavg*(Vol*VolFac)/(N*Tavg);
    fprintf(afp,"  Total Time (s)      T (K)               P (Pa)      PV/nT (J/(mol K))         Z           V (m^3)              N\n");
    fprintf(afp," --------------   -----------        ---------------   --------------   ---------------   ------------   -----------\n");
    fprintf(afp,"  %8.4e  %15.5f       %15.5f     %10.5f       %10.5f        %10.5e         %i\n",i*dt*timefac,Tavg,Pavg,gc,Z,Vol*VolFac,N);
    
    printf("\n  TO ANIMATE YOUR SIMULATION, OPEN THE FILE \n  '%s' WITH VMD AFTER THE SIMULATION COMPLETES\n",tfn);
    printf("\n  TO ANALYZE INSTANTANEOUS DATA ABOUT YOUR MOLECULE, OPEN THE FILE \n  '%s' WITH YOUR FAVORITE TEXT EDITOR OR IMPORT THE DATA INTO EXCEL\n",ofn);
    printf("\n  THE FOLLOWING THERMODYNAMIC AVERAGES WILL BE COMPUTED AND WRITTEN TO THE FILE  \n  '%s':\n",afn);
    printf("\n  AVERAGE TEMPERATURE (K):                 %15.5f\n",Tavg);
    printf("\n  AVERAGE PRESSURE  (Pa):                  %15.5f\n",Pavg);
    printf("\n  PV/nT (J * mol^-1 K^-1):                 %15.5f\n",gc);
    printf("\n  PERCENT ERROR of pV/nT AND GAS CONSTANT: %15.5f\n",100*fabs(gc-8.3144598)/8.3144598);
    printf("\n  THE COMPRESSIBILITY (unitless):          %15.5f \n",Z);
    printf("\n  TOTAL VOLUME (m^3):                      %10.5e \n",Vol*VolFac);
    printf("\n  NUMBER OF PARTICLES (unitless):          %i \n", N);
    
    
    
    
    fclose(tfp);
    fclose(ofp);
    fclose(afp);
    
    return 0;
}


void initialize() {
    // Number of atoms in each direction
    int n = int(ceil(pow(N, 1.0/3)));
    
    //  spacing between atoms along a given direction
    double pos = L / n;
    
    //  index for number of particles assigned positions
    int p = 0;
    //  initialize positions
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            for (int k=0; k<n; k++) {
                if (p<N) {
                    r_x[p] = (i + 0.5)*pos;
                    r_y[p] = (j + 0.5)*pos;
                    r_z[p] = (k + 0.5)*pos;
                }
                p++;
            }
        }
    }
    
    // Call function to initialize velocities
    initializeVelocities();
    
    /***********************************************
     *   Uncomment if you want to see what the initial positions and velocities are
     printf("  Printing initial positions!\n");
     for (i=0; i<N; i++) {
     printf("  %6.3e  %6.3e  %6.3e\n",r[i][0],r[i][1],r[i][2]);
     }
     
     printf("  Printing initial velocities!\n");
     for (i=0; i<N; i++) {
     printf("  %6.3e  %6.3e  %6.3e\n",v[i][0],v[i][1],v[i][2]);
     }
     */
    
    
    
}   


//  Function to calculate the averaged velocity squared
double MeanSquaredVelocity() { 
    
    double vx2 = 0;
    double vy2 = 0;
    double vz2 = 0;
    double v2;
    
    for (int i=0; i<N; i++) {
        
        vx2 = vx2 + v_x[i]*v_x[i];
        vy2 = vy2 + v_y[i]*v_y[i];
        vz2 = vz2 + v_z[i]*v_z[i];
        
    }
    v2 = (vx2+vy2+vz2)/N;
    
    
    //printf("  Average of x-component of velocity squared is %f\n",v2);
    return v2;
}

//  Function to calculate the kinetic energy of the system
double Kinetic() { //Write Function here!  
    
    
    double kin = 0.;
    for (int i=0; i<N; i++) {
        
        double vxi = v_x[i] * v_x[i];
        double vyi = v_y[i] * v_y[i];
        double vzi = v_z[i] * v_z[i];
        double v2 = vxi + vyi + vzi;
        kin += m*v2/2.;
        
    }
    
    //printf("  Total Kinetic Energy is %f\n",N*mvs*m/2.);
    return kin;
    
}


// Function to calculate the potential energy of the system
double Potential() {
    /*
    The potential energy between two particles i and j is symmetric, meaning U(i,j) = U(j,i).
    So, we can eliminate almost half of the calculations by only computing the potential energy for unique pairs of particles, without repeating.
    */
    double Pot = 0.;

    for (int i=0; i < N-1; i++) {
        for(int j = i+1; j < N; j++) {

            // we unroll the k loop: computes the difference for each of the coordinates between the two particles
            // we square the differences of each coordinate between the two particles and we sum them up.
            double rij_0 = r_x[i] - r_x[j]; // x distance between the two particles
            double rij_1 = r_y[i] - r_y[j]; // y distance between the two particles
            double rij_2 = r_z[i] - r_z[j]; // z distance between the two particles
            double sq0 = rij_0*rij_0;
            double sq1 = rij_1*rij_1;
            double sq2 = rij_2*rij_2;
            // This is the square of the distance between the two particles
            double distance_sqrd = sq0 + sq1 + sq2;

            /*
            Instead the following calculation:    quot = sigma / sqrt(distance_sqrd)

            We calculate:    quot^2 = (sigma / sqrt(distance_sqrd))^2
                                    = sigma^2 / distance_sqrd
                                    = sigma * sigma / distance_sqrd

            This way, we avoid calculating the square root of distance_sqrd, which is an expensive calculation

            Due to the fact that after calculating quot we need to calculate quot^6 and quot^12, it is actually useful to calculate quot^2
            directly and just calculating the other two powers from quot^2, which is what we do below:
            */
            
            double quot2 = sigma2 / distance_sqrd;
            /*
            Here we calculate quot^6 and quot^12 by calling the function 'pow', but with integer exponents because the function runs faster.
            Calculating these values this way is way faster because 
            */
            double quot6 = quot2*quot2*quot2;  //quot6 = (quot^2)^3
            double quot12 = quot6*quot6; // quot12 = quot6^2

            // We multiply by 8 instead of 4 to account for the symmetry property when calculating the Potential
            // Every pair of particle is only processed once instead of twice this way
            Pot += epsilon8 * (quot12 - quot6);  
        }
    }

    return Pot;
}

//inline double calculatePotential(int i, int j) {
//    double diff1 = r_x[i] - r_x[j];
//    double diff2 = r_y[i] - r_y[j];
//    double diff3 = r_z[i] - r_z[j];
//    double r2 = diff1 * diff1 + diff2 * diff2 + diff3 * diff3;
//    //double rnorm = sqrt(r2);
//    //double quot = sigma / rnorm;
//
//    //like this we avoid calculating the sqrt of r2, which is a costly operation
//    double quot = sigma * sigma / r2;
//    double quot6 = std::pow(quot, 6);
//    double quot12 = quot6 * quot6;
//    return epsilon8 * (quot12 - quot6);
//}
//
//double Potential() {
//    double Pot = 0.;
//
//    for (int i = 0; i < N-1; i++) {
//        // Process pairs of elements in the inner loop
//        int j = i+1;
//        // Process groups of 8 elements in the inner loop
//        for(; j+7 < N; j+=8) {
//            Pot += calculatePotential(i, j);
//            Pot += calculatePotential(i, j+1);
//            Pot += calculatePotential(i, j+2);
//            Pot += calculatePotential(i, j+3);
//            Pot += calculatePotential(i, j+4);
//            Pot += calculatePotential(i, j+5);
//            Pot += calculatePotential(i, j+6);
//            Pot += calculatePotential(i, j+7);
//        }
//
//        // Process any remaining elements
//        for(; j < N; j++) {
//            Pot += calculatePotential(i, j);
//        }
//    }
//
//    return Pot;
//}

//   Uses the derivative of the Lennard-Jones potential to calculate
//   the forces on each atom.  Then uses a = F/m to calculate the
//   accelleration of each atom. 
void computeAccelerations() {

    for (int i = 0; i < N; i++) {  // set all accelerations to zero
        // we unroll the k loop
        a_x[i] = 0;
        a_y[i] = 0;
        a_z[i] = 0;
    }

    // Declare the variables outside the parallel region
    double ai0, ai1, ai2;
    for (int i = 0; i < N-1; i++) {   // loop over all distinct pairs i,j

        // we define temporary variables so that we don't access the memory so often in the j loop
        ai0 = ai1 = ai2 = 0;

        #pragma omp parallel for schedule(static,100) reduction(+:ai0) reduction(+:ai1) reduction(+:ai2) num_threads(nworkers)
        for (int j = i+1; j < N; j++) {
            
            // we unroll the k loop: computes the difference for each of the coordinates between the two particles
            // we square the differences of each coordinate between the two particles and we sum them up.
            double rij_0 = r_x[i] - r_x[j]; // x distance between the two particles
            double rij_1 = r_y[i] - r_y[j]; // y distance between the two particles
            double rij_2 = r_z[i] - r_z[j]; // z distance between the two particles
            double sq0 = rij_0*rij_0;
            double sq1 = rij_1*rij_1;
            double sq2 = rij_2*rij_2;
            // This is the square of the distance between the two particles
            double distance_sqrd = sq0 + sq1 + sq2;

            // d1 is the inverse of the square of the distance
            double d1 = 1 / distance_sqrd;
            // d2 = d1^2
            double d2 = d1*d1;
            // d4 = d2^2 = d1^4
            double d4 = d2*d2;
            
            //  From derivative of Lennard-Jones with sigma and epsilon set equal to 1 in natural units!
            //double f = 24 * (2 * pow(distance_sqrd, -7) - pow(distance_sqrd, -4));

            // pow(distance_sqrd, -7) = d1^7 = d1^4 * d1^2 * d1 = d4*d2*d1
            // pow(distance_sqrd, -4) = d1^4 = d4
            double f = 24 * (2 * d4*d2*d1 - d4);


            // we unroll the k loop
            // Since the value of i doesn't change inside the j loop, we can access the temporary variables instead of accessing the memory each iteration
            //#pragma omp atomic
            ai0 += rij_0 * f;
            //#pragma omp atomic
            ai1 += rij_1 * f;
            //#pragma omp atomic
            ai2 += rij_2 * f;

            //a_x[i] += rij_0 * f;
            //a_y[i] += rij_1 * f;
            //a_z[i] += rij_2 * f;

            //#pragma omp atomic
            a_x[j] -= rij_0 * f;
            //#pragma omp atomic
            a_y[j] -= rij_1 * f;
            //#pragma omp atomic
            a_z[j] -= rij_2 * f;
        }
        // we update the memory with the value of the temporary variables calculated inside the j loop
        //#pragma omp atomic
        a_x[i] += ai0;
        //#pragma omp atomic
        a_y[i] += ai1;
        //#pragma omp atomic
        a_z[i] += ai2;
    }
}

// returns sum of dv/dt*m/A (aka Pressure) from elastic collisions with walls
double VelocityVerlet(double dt, int iter, FILE *fp) {
    int i, j, k;
    
    double psum = 0.;
    
    //  Compute accelerations from forces at current position
    // this call was removed (commented) for predagogical reasons
    //computeAccelerations();
    //  Update positions and velocity with current velocity and acceleration
    //printf("  Updated Positions!\n");
    for (i=0; i<N; i++) {
        r_x[i] += v_x[i]*dt + 0.5*a_x[i]*dt*dt;
        v_x[i] += 0.5*a_x[i]*dt;
        
        r_y[i] += v_y[i]*dt + 0.5*a_y[i]*dt*dt;
        v_y[i] += 0.5*a_y[i]*dt;
        
        r_z[i] += v_z[i]*dt + 0.5*a_z[i]*dt*dt;
        v_z[i] += 0.5*a_z[i]*dt;


        //printf("  %i  %6.4e   %6.4e   %6.4e\n",i,r[i][0],r[i][1],r[i][2]);
    }
    //  Update accellerations from updated positions
    computeAccelerations();
    //  Update velocity with updated acceleration
    for (i=0; i<N; i++) {
        v_x[i] += 0.5*a_x[i]*dt;
        v_y[i] += 0.5*a_y[i]*dt;
        v_z[i] += 0.5*a_z[i]*dt;
    }
    
    // Elastic walls
    for (i=0; i<N; i++) {
        if (r_x[i] < 0 || r_x[i] >= L) {
            v_x[i] = -v_x[i]; //- elastic walls
            psum += 2*m*fabs(v_x[i])/dt;  // contribution to pressure from "left" walls
        }
        
        if (r_y[i] < 0 || r_y[i] >= L) {
            v_y[i] = -v_y[i]; //- elastic walls
            psum += 2*m*fabs(v_y[i])/dt;  // contribution to pressure from "left" walls
        }
        
        if (r_z[i] < 0 || r_z[i] >= L) {
            v_z[i] = -v_z[i]; //- elastic walls
            psum += 2*m*fabs(v_z[i])/dt;  // contribution to pressure from "left" walls
        }
    }
    
    return psum/(6*L*L);
}


void initializeVelocities() {
    
    int i, j;
    
    for (i=0; i<N; i++) {
        
        //  Pull a number from a Gaussian Distribution
        v_x[i] = gaussdist();
        v_y[i] = gaussdist();
        v_z[i] = gaussdist();
    }
    
    // Vcm = sum_i^N  m*v_i/  sum_i^N  M
    // Compute center-of-mas velocity according to the formula above
    double vCM[3] = {0, 0, 0};
    
    for (i=0; i<N; i++) {
        vCM[0] += m*v_x[i];
        vCM[1] += m*v_y[i];
        vCM[2] += m*v_z[i];
    }
    
    
    for (i=0; i<3; i++) vCM[i] /= N*m;
    
    //  Subtract out the center-of-mass velocity from the
    //  velocity of each particle... effectively set the
    //  center of mass velocity to zero so that the system does
    //  not drift in space!
    for (i=0; i<N; i++) {
        v_x[i] -= vCM[0];
        v_y[i] -= vCM[1];
        v_z[i] -= vCM[2];
    }
    
    //  Now we want to scale the average velocity of the system
    //  by a factor which is consistent with our initial temperature, Tinit
    double vSqdSum = 0;
    for (i=0; i<N; i++) {
        vSqdSum += v_x[i]*v_x[i];
        vSqdSum += v_y[i]*v_y[i];
        vSqdSum += v_z[i]*v_z[i];
    }
    
    double lambda = sqrt( 3*(N-1)*Tinit/vSqdSum);
    
    for (i=0; i<N; i++) {
        v_x[i] *= lambda;
        v_y[i] *= lambda;
        v_z[i] *= lambda;
    }
}


//  Numerical recipes Gaussian distribution number generator
double gaussdist() {
    static bool available = false;
    static double gset;
    double fac, rsq, v1, v2;
    if (!available) {
        do {
            v1 = 2.0 * rand() / double(RAND_MAX) - 1.0;
            v2 = 2.0 * rand() / double(RAND_MAX) - 1.0;
            rsq = v1 * v1 + v2 * v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        
        fac = sqrt(-2.0 * log(rsq) / rsq);
        gset = v1 * fac;
        available = true;
        
        return v2*fac;
    } else {
        
        available = false;
        return gset;
        
    }
}

// efficient recursive function that implements exponentiation by squaring
/*
x^n is
if even: x(^(n/2))^2
if odd: x * x^(n-1)

does in log(N)
*/
double pow_by_squaring(double x, int n) {
    if (n == 0) {
        return 1.0;
    } else if (n < 0) {
        return 1.0 / (x * pow_by_squaring(x, -n - 1));
    } else if (n % 2 == 0) {  // n is even
        double y = pow_by_squaring(x, n/2);
        return y * y;
    } else {  // n is odd
        return x * pow_by_squaring(x, n-1);
    }
}
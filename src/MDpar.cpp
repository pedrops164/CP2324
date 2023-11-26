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
#include <omp.h>
#include "data_structures.h"

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
//  Position data structure
Position position;
//  Velocity data structure
Velocity velocity;
//  Acceleration data structure
Acceleration acceleration;

// atom type
char atype[10];
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
    // We don't need the Potential yet, so we don't store it.
    computeAccelerationsAndPotential();
    
    
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
        Press = VelocityVerlet(dt, i+1, tfp, &PE);
        Press *= PressFac;
        
        //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //  Now we would like to calculate somethings about the system:
        //  Instantaneous mean velocity squared, Temperature, Pressure
        //  Potential, and Kinetic Energy
        //  We would also like to use the IGL to try to see if we can extract the gas constant
        mvs = MeanSquaredVelocity();
        KE = Kinetic();
        //PE = Potential();
        
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
                    position.x[p] = (i + 0.5)*pos;
                    position.y[p] = (j + 0.5)*pos;
                    position.z[p] = (k + 0.5)*pos;
                }
                p++;
            }
        }
    }
    
    // Call function to initialize velocities
    initializeVelocities();
}   


//  Function to calculate the averaged velocity squared
double MeanSquaredVelocity() { 
    double vx2 = 0;
    double vy2 = 0;
    double vz2 = 0;
    double v2;
    
    for (int i=0; i<N; i++) {
        vx2 = vx2 + velocity.x[i]*velocity.x[i];
        vy2 = vy2 + velocity.y[i]*velocity.y[i];
        vz2 = vz2 + velocity.z[i]*velocity.z[i];
    }
    v2 = (vx2+vy2+vz2)/N;
    return v2;
}

//  Function to calculate the kinetic energy of the system
double Kinetic() {
    
    
    double kin = 0.;
    for (int i=0; i<N; i++) {
        
        double vxi = velocity.x[i] * velocity.x[i];
        double vyi = velocity.y[i] * velocity.y[i];
        double vzi = velocity.z[i] * velocity.z[i];
        double v2 = vxi + vyi + vzi;
        kin += m*v2/2.;
        
    }
    
    //printf("  Total Kinetic Energy is %f\n",N*mvs*m/2.);
    return kin;
    
}

inline double calculateLocalPotential(double inverseSquaredDistance) {
    /*
    Given the inverse of the square of the distance between two particles, returns the potential between them
    */

    /*
    Instead of the following calculation:    quot = sigma / sqrt(distance_sqrd)
    We calculate:    quot^2 = (sigma / sqrt(distance_sqrd))^2
                            = sigma^2 / distance_sqrd
                            = sigma * sigma / distance_sqrd
    This way, we avoid calculating the square root of distance_sqrd, which is an expensive calculation
    Due to the fact that after calculating quot we need to calculate quot^6 and quot^12, it is actually useful to calculate quot^2
    directly and just calculating the other two powers from quot^2, which is what we do below:
    */
    double quot2 = sigma2 * inverseSquaredDistance;
    /*
    We calculate these values with large exponents manually, to be as quick as possible
    */
    double quot6 = quot2*quot2*quot2;  //quot6 = (quot^2)^3
    double quot12 = quot6*quot6; // quot12 = quot6^

    // We multiply by 8 instead of 4 to account for the symmetry property when calculating the Potential
    // Every pair of particle should only be processed once instead of twice this way, to save computation time!
    return epsilon8 * (quot12 - quot6); 
}

// Uses the derivative of the Lennard-Jones potential to calculate
// the forces on each atom.  Then uses a = F/m to calculate the
// accelleration of each atom. 
// Returns the Potential
// We calculate acceleration and Potential in the same function to speed up computation (two rabbits one stone)
double computeAccelerationsAndPotential() {

    for (int i = 0; i < N; i++) {  // set all accelerations to zero
        // we unroll the k loop
        acceleration.x[i] = 0;
        acceleration.y[i] = 0;
        acceleration.z[i] = 0;
    }

    double Pot = 0.;
    // Declare the variables outside the parallel region
    // we set the arrays inside the acceleration structure as local variables in order to be able to use them in
    // the openmp reduction clause!! Otherwise it wouldn't work
    double * acc_x = acceleration.x;
    double * acc_y = acceleration.y;
    double * acc_z = acceleration.z;
    /*
    Dynamic scheduling is used because the iterations have different sizes (the range of the inner loop depends on
    the iteration of the outer loop). We chose a chunksize of 50 because it is a good tradeoff between having good
    fine-grain distribution, while minimizing the time spent assigning chunks to workers.
    The reduction clause is used to avoid data races in the calculated Potential value, and also in the acceleration
    data structure.
    */
    #pragma omp parallel for schedule(dynamic,50) reduction(+:Pot,acc_x[:N],acc_y[:N],acc_z[:N])
    for (int i = 0; i < N-1; i++) {   // loop over all distinct pairs i,j
        // we define temporary variables so that we don't access the memory so often in the j loop
        double ai0=0, ai1=0, ai2=0;

        for (int j = i+1; j < N; j++) {
            
            // we unroll the k loop: computes the difference for each of the coordinates between the two particles
            // we square the differences of each coordinate between the two particles and we sum them up.
            double rij_0 = position.x[i] - position.x[j]; // x distance between the two particles
            double rij_1 = position.y[i] - position.y[j]; // y distance between the two particles
            double rij_2 = position.z[i] - position.z[j]; // z distance between the two particles
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

            // We pass as argument the inverse of the distance to avoid calculating the same division
            Pot += calculateLocalPotential(d1);
            
            //  From derivative of Lennard-Jones with sigma and epsilon set equal to 1 in natural units!
            //double f = 24 * (2 * pow(distance_sqrd, -7) - pow(distance_sqrd, -4));
            // -> pow(distance_sqrd, -7) = d1^7 = d1^4 * d1^2 * d1 = d4*d2*d1
            // -> pow(distance_sqrd, -4) = d1^4 = d4
            double f = 24 * (2 * d4*d2*d1 - d4);
            double rij0f = rij_0 * f;
            double rij1f = rij_1 * f;
            double rij2f = rij_2 * f;

            // Since the value of i doesn't change inside the j loop, we can access the temporary variables
            // instead of accessing the memory each iteration, avoiding memory access overhead
            ai0 += rij0f;
            ai1 += rij1f;
            ai2 += rij2f;

            acc_x[j] -= rij0f;
            acc_y[j] -= rij1f;
            acc_z[j] -= rij2f;
        }
        // we update the memory with the value of the temporary variables calculated inside the j loop
        acc_x[i] += ai0;
        acc_y[i] += ai1;
        acc_z[i] += ai2;
    }
    return Pot;
}

// returns sum of dv/dt*m/A (aka Pressure) from elastic collisions with walls
double VelocityVerlet(double dt, int iter, FILE *fp, double * PE) {
    int i, j, k;
    
    double psum = 0.;
    
    //  Compute accelerations from forces at current position
    //  Update positions and velocity with current velocity and acceleration
    for (i=0; i<N; i++) {
        position.x[i] += velocity.x[i]*dt + 0.5*acceleration.x[i]*dt*dt;
        velocity.x[i] += 0.5*acceleration.x[i]*dt;
        
        position.y[i] += velocity.y[i]*dt + 0.5*acceleration.y[i]*dt*dt;
        velocity.y[i] += 0.5*acceleration.y[i]*dt;
        
        position.z[i] += velocity.z[i]*dt + 0.5*acceleration.z[i]*dt*dt;
        velocity.z[i] += 0.5*acceleration.z[i]*dt;
    }
    //  Update accellerations from updated positions, and calculate potential
    *PE = computeAccelerationsAndPotential();
    //  Update velocity with updated acceleration
    for (i=0; i<N; i++) {
        velocity.x[i] += 0.5*acceleration.x[i]*dt;
        velocity.y[i] += 0.5*acceleration.y[i]*dt;
        velocity.z[i] += 0.5*acceleration.z[i]*dt;
    }
    
    // Elastic walls
    for (i=0; i<N; i++) {
        if (position.x[i] < 0 || position.x[i] >= L) {
            velocity.x[i] = -velocity.x[i]; //- elastic walls
            psum += 2*m*fabs(velocity.x[i])/dt;  // contribution to pressure from "left" walls
        }
        
        if (position.y[i] < 0 || position.y[i] >= L) {
            velocity.y[i] = -velocity.y[i]; //- elastic walls
            psum += 2*m*fabs(velocity.y[i])/dt;  // contribution to pressure from "left" walls
        }
        
        if (position.z[i] < 0 || position.z[i] >= L) {
            velocity.z[i] = -velocity.z[i]; //- elastic walls
            psum += 2*m*fabs(velocity.z[i])/dt;  // contribution to pressure from "left" walls
        }
    }
    
    return psum/(6*L*L);
}


void initializeVelocities() {
    
    int i, j;
    
    for (i=0; i<N; i++) {
        
        //  Pull a number from a Gaussian Distribution
        velocity.x[i] = gaussdist();
        velocity.y[i] = gaussdist();
        velocity.z[i] = gaussdist();
    }
    
    // Vcm = sum_i^N  m*v_i/  sum_i^N  M
    // Compute center-of-mas velocity according to the formula above
    double vCM[3] = {0, 0, 0};
    
    for (i=0; i<N; i++) {
        vCM[0] += m*velocity.x[i];
        vCM[1] += m*velocity.y[i];
        vCM[2] += m*velocity.z[i];
    }
    
    
    for (i=0; i<3; i++) vCM[i] /= N*m;
    
    //  Subtract out the center-of-mass velocity from the
    //  velocity of each particle... effectively set the
    //  center of mass velocity to zero so that the system does
    //  not drift in space!
    for (i=0; i<N; i++) {
        velocity.x[i] -= vCM[0];
        velocity.y[i] -= vCM[1];
        velocity.z[i] -= vCM[2];
    }
    
    //  Now we want to scale the average velocity of the system
    //  by a factor which is consistent with our initial temperature, Tinit
    double vSqdSum = 0;
    for (i=0; i<N; i++) {
        vSqdSum += velocity.x[i]*velocity.x[i];
        vSqdSum += velocity.y[i]*velocity.y[i];
        vSqdSum += velocity.z[i]*velocity.z[i];
    }
    
    double lambda = sqrt( 3*(N-1)*Tinit/vSqdSum);
    
    for (i=0; i<N; i++) {
        velocity.x[i] *= lambda;
        velocity.y[i] *= lambda;
        velocity.z[i] *= lambda;
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
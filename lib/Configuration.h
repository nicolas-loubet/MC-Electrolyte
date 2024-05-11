#pragma once
#include "Particle.h"
#include <iostream>
#include <ctime>
#include <random>
#include <chrono>
#include <algorithm>
#include <vector>

/**
 * This class creates a Configuration object, with an array of Water objects
 */
class Configuration {
	protected:
		Particle** ps; //The pointer to the pointer to the first element of the array
		int N_PARTICLE= 0; //The number of Water objects in the array
		float TEMPERATURE= 0.; //Temperature of the system
		Vector* bounds; //The bounds of the system
        float relative_permitivity= 78.;
        float R_SPHERE= 0.;

        /**
         * Creates a random float between bounds. srand(time(NULL)) should be used before using this function.
         * @param LO Minimum number (included)
         * @param HI Maximum number (not included)
         * @return A random number, of type float
         */
        float randomBetween(const float LO, const float HI) {
            return LO + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(HI-LO)));
        }

        /**
         * Checks if the particle is in the box, and edit the coordinates if neccesary
         * @param v Vector* to the vector to check
         */
        void checkPBC(Vector* v) {
            if(v->x >= bounds->x)   v->x-= bounds->x;
            else if(v->x < 0.)      v->x+= bounds->x;
            if(v->y >= bounds->y)   v->y-= bounds->y;
            else if(v->y < 0.)      v->y+= bounds->y;
            if(v->z >= bounds->z)   v->z-= bounds->z;
            else if(v->z < 0.)      v->z+= bounds->z;
        }

        /**
         * Computes the proportion of `true` values in a boolean array.
         * @param arr Pointer to the boolean array.
         * @param N Number of elements in the array.
         * @return The proportion of `true` values in the array, represented as a float.
         */
        float getProportion(bool* arr, const int N) {
            float f= 0.;
            for(int i= 0; i < N; i++)
                if(arr[i]) f++;
            return f/N;
        }

        /**
         * Performs a single iteration of a Monte Carlo simulation for particle movement.
         * It randomly selects particles and attempts to move them within the system, accepting or rejecting the move based on the change in energy and a Metropolis criterion.
         *
         * @param gen Reference to a Mersenne Twister random number generator.
         * @param r_max_adjust Maximum adjustment factor for the random displacement of particles.
         * @return The number of particles that successfully changed their positions during the iteration.
         */
        int iteration(mt19937& gen, const float r_max_adjust) {
            const float k_B= 1.380649E-23;

            vector<int> list_shuffle;
            for(int i= 1; i <= N_PARTICLE; i++)
                list_shuffle.push_back(i);
            shuffle(list_shuffle.begin(), list_shuffle.end(), gen);

            int changed= 0;
            for(int i: list_shuffle) {
                Vector* dq_not_norm= new Vector(randomBetween(-.5,.5),randomBetween(-.5,.5),randomBetween(-.5,.5));
                Vector* dq= *dq_not_norm * r_max_adjust;
                delete(dq_not_norm);

                float E_ref= getEnergyParticle(i);

                Vector* pos_new= *ps[i-1]->getPosition() + *dq;
                checkPBC(pos_new);
                delete(dq);
                Vector* pos_prev= ps[i-1]->getPosition();
                ps[i-1]->setPosition(pos_new);

                float E_ens= getEnergyParticle(i);

                if(randomBetween(0,1) > exp(-(E_ens-E_ref)/(k_B*TEMPERATURE))) { //Decline
                    ps[i-1]->setPosition(pos_prev);
                    delete(pos_new);
                } else { //Accept
                    changed++;
                    delete(pos_prev);
                }
            }
            return changed;
        }

	public:
        //Getters
		int getNParticle() { return N_PARTICLE; }
		Particle* getParticle(int id) { return ps[id-1]; }
		Vector* getBounds() { return bounds; }

        /**
         * Constructs a Configuration object with specified parameters.
         * This constructor initializes a Configuration object with the specified number of particles, spherical boundary radius, temperature, periodic boundary conditions, and relative permittivity.
         * It creates an array of Particle pointers, initializes each Particle with a unique identifier, random position within the given periodic boundary conditions, and assigns a charge based on its position in the particle array.
         * Additionally, it ensures that there are no particle superpositions by randomly repositioning particles if their initial positions are too close.
         * @param n_part The number of particles in the configuration.
         * @param r_sphere The radius of the spherical boundary that confines the particles.
         * @param temp The temperature of the system.
         * @param pbc Pointer to a Vector representing the periodic boundary conditions.
         * @param rel_perm The relative permittivity of the medium.
         */
		Configuration(const int n_part, const float r_sphere, const float temp, Vector* pbc, const float rel_perm) {
			bounds= pbc;
			N_PARTICLE= n_part;
			TEMPERATURE= temp;
			ps= new Particle*[n_part];
			relative_permitivity= rel_perm;
			R_SPHERE= r_sphere;

			for(int i= 0; i < n_part; i++)
                ps[i]= new Particle(i+1,new Vector(randomBetween(0.,pbc->x),randomBetween(0.,pbc->y),randomBetween(0.,pbc->z)),(i < n_part/2 ? +1:-1));

            //Check for superpositions
            for(int i= 0; i < n_part; i++)
                for(int j= 0; j < n_part; j++) {
                    if(i == j) continue;
                    if(ps[i]->distanceTo(ps[j],bounds) < R_SPHERE) {
                        ps[i]->getPosition()->x= randomBetween(0.,pbc->x);
                        ps[i]->getPosition()->y= randomBetween(0.,pbc->y);
                        ps[i]->getPosition()->z= randomBetween(0.,pbc->z);
                        j= 0;
                        i--;
                        if(i < 0) i= 0;
                    }
                }
		}

        /**
         * Destructor. It destroys the particle array and bounds.
         */
		~Configuration() {
			for(int i= 0; i < N_PARTICLE; i++)
				delete(ps[i]);
			delete(ps);
			delete(bounds);
		}

		/**
         * Prints configuration to an XYZ file. Uses append, so you should check if the file already exists.
         * @param file_name The name of the file to which the particle positions will be written.
         * @param comment A comment string that will be written as a header in the XYZ file.
         */
        void printXYZ(string file_name, const string comment) {
            const string SEP= "\t";
            ofstream f(file_name,ios::app);

            f << N_PARTICLE << endl;
            f << comment << endl;
            for(int i= 0; i < N_PARTICLE; i++)
                f << (ps[i]->getCharge() > 0 ? "Na":"Cl") << SEP << ps[i]->getPosition()->x << SEP << ps[i]->getPosition()->y << SEP << ps[i]->getPosition()->z << endl;

            f.flush();
            f.close();
		}

		/**
		 * Computes the energy for the ID particle.
		 * @param ID The identifier (1-N) of the particle.
		 * @return sum of all the pair-waise energies involving particle ID, or 10E100 if overlap.
		 */
		float getEnergyParticle(const int ID) {
            const float const_energy= 2.3066E-19/relative_permitivity;
            float output= 0.;
            for(int j= 0; j < N_PARTICLE; j++) {
                if(ID-1 == j) continue;
                if(ps[ID-1]->distanceTo(ps[j],bounds) < R_SPHERE) return 10E100;
                output+= const_energy*(ps[ID-1]->getCharge()*ps[j]->getCharge())/ps[ID-1]->distanceTo(ps[j],bounds);
            }
            return output;
		}

		/**
		 * Computes the energy of the whole system (configuration)
		 */
		float getEnergy() {
            const float const_energy= 2.3066E-19/relative_permitivity;
            float output= 0.;
            for(int i= 0; i < N_PARTICLE; i++)
                for(int j= i+1; j < N_PARTICLE; j++)
                    output+= const_energy*(ps[i]->getCharge()*ps[j]->getCharge())/ps[i]->distanceTo(ps[j],bounds);
            return output;
		}

		/**
         * Performs equilibration of the system using Monte Carlo simulation for a specified number of iterations. It adjusts the maximum displacement factor for particle movement dynamically based on the acceptance rate of particle movements.
         * @param N_ITER The number of iterations for equilibration.
         * @param R_MAX The maximum adjustment factor for particle displacement (default: 1.0).
         * @param print_every The interval for printing system energy during equilibration (default: 200).
         * @return The adjusted maximum displacement factor after equilibration.
         */
        float equilibrate(const int N_ITER, const float R_MAX= 1., const int print_every= 200) {
            const float R_ADJUST_RATIO= .01;
            const float MIN_PROPORTION_INCREASE_R= .2;
            const float MAX_PROPORTION_DECREASE_R= .2;
            const int CHECK_RMAX_EVERY= 100;

            random_device rd;
            mt19937 gen(rd());
            uniform_int_distribution<int> dis(1, N_PARTICLE);

            float r_max_adjust= R_MAX;
            int n_ens= 0;
            int n_ens_ac= 0;
            float f= 0;

            chrono::high_resolution_clock::time_point start_time= chrono::high_resolution_clock::now();

            for(int iter= 0; iter <= N_ITER; iter++) {
                n_ens_ac+= iteration(gen,r_max_adjust);
                n_ens+= N_PARTICLE;

                if(iter%CHECK_RMAX_EVERY == 0 && iter > 0) {
                    f= float(n_ens_ac)/n_ens;
                    if(f > MIN_PROPORTION_INCREASE_R) r_max_adjust*= 1+R_ADJUST_RATIO;
                    else if(f <= MAX_PROPORTION_DECREASE_R) r_max_adjust*= 1-R_ADJUST_RATIO;
                    if(r_max_adjust > bounds->magnitude()/2) r_max_adjust*= 1-R_ADJUST_RATIO*2;
                    n_ens= 0;
                    n_ens_ac= 0;
                }

                if(iter % print_every == 0) cout << "Energy = " << getEnergy() << " for step " << iter << " using R=" << r_max_adjust << " with acceptance rate f=" << f << endl;
            }

            cout << endl << "------------------------------------------------" << endl;
            chrono::high_resolution_clock::time_point end_time= chrono::high_resolution_clock::now();
            chrono::duration<double> elapsed_seconds= end_time - start_time;
            cout << "Time for equilibration: " << elapsed_seconds.count() << " s" << endl;
            cout << "Final energy for equilibration: " << getEnergy() << " J" << endl;
            cout << "------------------------------------------------" << endl;

            return r_max_adjust;
        }

		/**
         * Performs production runs of the MC simulation for a specified number of iterations, performing particle movements and periodically printing the configuration to an XYZ file.
         * @param file_name The name of the file to which configurations will be written.
         * @param N_ITER The number of iterations for the production run.
         * @param print_every The interval for printing configurations to the file (default: 200).
         * @param R_MAX The maximum adjustment factor for particle displacement (default: 1.0).
         */
        void production(const string file_name, const int N_ITER, const int print_every= 200, const float R_MAX= 1.) {
            random_device rd;
            mt19937 gen(rd());
            uniform_int_distribution<int> dis(1, N_PARTICLE);

            int counter= 0;

            chrono::high_resolution_clock::time_point start_time= chrono::high_resolution_clock::now();

            for(int iter= 0; iter <= N_ITER; iter++) {
                iteration(gen,R_MAX);
                if(iter % print_every == 0) {
                    cout << "Printing step " << iter << endl;
                    printXYZ(file_name,"Step "+to_string(iter));
                    counter++;
                }
            }

            cout << endl << "------------------------------------------------" << endl;
            chrono::high_resolution_clock::time_point end_time= chrono::high_resolution_clock::now();
            chrono::duration<double> elapsed_seconds= end_time - start_time;
            cout << "Time for production: " << elapsed_seconds.count() << " s" << endl;
            cout << "Final energy for production: " << getEnergy() << " J" << endl;
            cout << "Number of configurations: " << counter << endl;
            cout << "------------------------------------------------" << endl;

        }

};

#pragma once
#include "Util.h"

/**
 * This class creates a Particle object, with an ID to identify it and a position in a 3D system
 */
class Particle {
	protected:
		Vector* pos;
		int ID= -1;
		float q= 0;

	public:
        //Setters and getters
		void setPosition(Vector* position) { pos= position; }
		Vector* getPosition() { return pos; }
		void setID(int i) { ID= i; }
		int getID() { return ID; }
		void setCharge(float charge) { q= charge; }
		int getCharge() { return q; }

        /**
         * Basic constructor
         * @param id Is the ID to identify the Particle
         * @param position Is the *Vector of the position
         * @param charge Is the charge of the particle for Coulombic potential
         */
		Particle(int id, Vector* position, float charge= 0) {
			pos= position;
			ID= id;
			q= charge;
		}

        /**
         * Constructor for assigment
         * @param p Other Particle already defined
         */
        Particle(const Particle& p): ID(p.ID), pos(p.pos), q(p.q) {}

        /**
         * Destructor. It destroys the position if not already destroyed.
         */
		virtual ~Particle() {
            if(pos != NULL) {
                delete(pos);
                pos= NULL;
            }
		}

		/**
         * It measures the distance to another Particle
         * @param p *Particle to other particle
         * @param bounds The coordinate of the last point, so the components are the width, height and length
         */
		float distanceTo(Particle* p, Vector* bounds) {
			return dist(this->pos, p->pos, bounds);
		}

};

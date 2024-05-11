#pragma once
#include <ostream>
#include <math.h>
using namespace std;

/**
 * This class creates a mathematical Vector, with components x, y and z (3D system)
 */
class Vector {
	public:
		float x, y, z;

        /**
         * Basic constructor
         * @param v_x Is the x component of the vector (default=0)
         * @param v_y Is the y component of the vector (default=0)
         * @param v_z Is the z component of the vector (default=0)
         */
		Vector(float v_x= 0, float v_y= 0, float v_z= 0) {
			x= v_x;
			y= v_y;
			z= v_z;
		}

        /**
         * Constructor for assigment
         * @param v Other vector already defined
         */
		Vector(const Vector& v): x(v.x), y(v.y), z(v.z) {}

		/**
         * @return module of the vector
         */
		float magnitude() const {
			return sqrt(x*x+y*y+z*z);
		}

		/**
         * To compare two vectors in magnitude. Otherwise, use ==
         * @param v Other vector to compare
         * @return true if both vectors have the same magnitude (be careful with the precision managed)
         */
		bool equalMagnitud(Vector& v) {
			return magnitude()==v.magnitude();
		}

		/**
         * OVERLOAD of operator + for two vectors
         * @param v Other vector to sum
         * @return The sum of each component in a new *Vector
         */
		Vector* operator +(const Vector& v) {
			return new Vector(x+v.x, y+v.y, z+v.z);
		}

		/**
         * OVERLOAD of operator - for two vectors
         * @param v Other vector to substract
         * @return The substraction of each component in a new *Vector
         */
		Vector* operator -(const Vector& v) {
			return new Vector(x-v.x, y-v.y, z-v.z);
		}

		/**
         * OVERLOAD of operator * for two vectors
         * @param v Other vector to operate
         * @return The scalar product of both vectors
         */
		float operator *(const Vector& v) {
			return x*v.x + y*v.y + z*v.z;
		}

		/**
         * OVERLOAD of operator * for one vector and a scalar
         * @param k A constant (scalar)
         * @return The product of each component with the scalar k in a new *Vector
         */
		Vector* operator *(const float k) {
			return new Vector(x*k, y*k, z*k);
		}

		/**
         * OVERLOAD of operator / for one vector and a scalar
         * @param k A constant (scalar)
         * @return The division of each component with the scalar k in a new *Vector
         */
		Vector* operator /(const float k) {
			return new Vector(x/k, y/k, z/k);
		}

        /**
         * OVERLOAD of operator % for two vectors
         * @param v Other vector to operate
         * @return The cross (vector) product of both vectors in a new *Vector
         */
		Vector* operator %(const Vector& v) {
			return new Vector(y*v.z-z*v.y, z*v.x-x*v.z, x*v.y-y*v.x);
		}

		/**
		 * OVERLOAD of operator == for two vectors
         * To compare two vectors in component. For comparing magnitud use equalsMagnitud()
         * @param v Other vector to compare
         * @return true if both vectors have the same components in x, y, and z (be careful with the precision managed)
         */
		bool operator ==(const Vector& v) {
			return x==v.x && y==v.y && z==v.z;
		}

		/**
		 * OVERLOAD of operator != for two vectors
         * To compare two vectors in component. For comparing magnitud use !equalsMagnitud()
         * @param v Other vector to compare
         * @return true if both vectors DO NOT have the same components in x, y, and z
         */
		bool operator !=(const Vector& v) {
			return !(*this==v);
		}

		/**
		 * OVERLOAD of operator < for two vectors
         * To compare two vectors in magnitud
         * @param v Other vector to compare
         * @return true if this vector has a minor magnitud than v
         */
		bool operator <(const Vector& v) {
			return magnitude()<v.magnitude();
		}

		/**
		 * OVERLOAD of operator > for two vectors
         * To compare two vectors in magnitud
         * @param v Other vector to compare
         * @return true if this vector has a major magnitud than v
         */
		bool operator >(const Vector& v) {
			return magnitude()>v.magnitude();
		}

		/**
		 * OVERLOAD of operator <= for two vectors
         * To compare two vectors in magnitud
         * @param v Other vector to compare
         * @return true if this vector has a minor or equal magnitud than v
         */
		bool operator <=(const Vector& v) {
			return magnitude()<=v.magnitude();
		}

		/**
		 * OVERLOAD of operator >= for two vectors
         * To compare two vectors in magnitud
         * @param v Other vector to compare
         * @return true if this vector has a major or equal magnitud than v
         */
		bool operator >=(const Vector& v) {
			return magnitude()>=v.magnitude();
		}

};

/**
 * OVERLOAD of operator << for a Vector
 * Writes the components of the vector in the format: -->(x,y,z)
 * @param o Ostream to write in
 * @param v Vector to write
 * @return Ostream with the writed string
 */
ostream& operator <<(ostream& o, Vector& v) {
    o << "-->(" << v.x << "," << v.y << "," << v.z << ")";
    return o;
}

/**
 * OVERLOAD of operator << for a *Vector
 * Writes the components of the vector in the format: -->(x,y,z)
 * @param o Ostream to write in
 * @param v Vector to write
 * @return Ostream with the writed string
 */
ostream& operator <<(ostream& o, Vector* v) {
    o << "-->(" << v->x << "," << v->y << "," << v->z << ")";
    return o;
}

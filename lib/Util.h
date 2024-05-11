#pragma once
#include "Vector.h"
#include <math.h>
#include <string>
#include <fstream>
using namespace std;

#define ERROR_VALUE -3

/**
 * STRUCT that contains the pointer to an array of int and the size of that array
 */
struct ArrInt {
    int *arr;
    int size= 0;
};

/**
 * Calculates the distance between two coordinates (vectors) given the bounds to use periodic boundary conditions
 * @param c1 One of the coordinates
 * @param c2 The other coordinate
 * @param bounds The coordinate of the last point, so the components are the width, height and length
 * @return The distance between the two coordinates
 */
float dist(Vector* c1, Vector* c2, Vector* bounds) {
    Vector* dif= *c1-*c2;
    Vector* b_1= *bounds/2;
    Vector* b_2= *bounds/(-2);

    if(dif->x > b_1->x) dif->x-= bounds->x;
    if(dif->x < b_2->x) dif->x+= bounds->x;
    if(dif->y > b_1->y) dif->y-= bounds->y;
    if(dif->y < b_2->y) dif->y+= bounds->y;
    if(dif->z > b_1->z) dif->z-= bounds->z;
    if(dif->z < b_2->z) dif->z+= bounds->z;

    float output= dif->magnitude();

    delete(dif);
    delete(b_1);
    delete(b_2);

    return output;
}

/**
 * Calculates the angle that 3 points form in a 3D system
 * @param c1 One of the point in the edges
 * @param c2 The center point
 * @param c3 The other point in an edge
 * @param bounds The coordinate of the last point, so the components are the width, height and length
 * @return The angle in radians formed by the 3 points
 */
float getAngle(Vector* c1, Vector* c2, Vector* c3, Vector* bounds) {
    float a= dist(c1,c3,bounds); //Opposite to the angle
    float b= dist(c1,c2,bounds);
    float c= dist(c2,c3,bounds);
    //Law of cosines
    return abs(acos((pow(b,2)+pow(c,2)-pow(a,2))/(2*b*c)));
}

/**
 * Function that returns the position in a set of bins of a value (for a histogram)
 * @param value The value to analyze
 * @param LIMIT_MIN The value of the first bin
 * @param LIMIT_MAX The value of the last bin
 * @param N_BINS The number of bins
 * @param exact If false, when the position if less than 0 or more than N_BINS, it returns the nearest position. If true, in the case of "out of bounds" the result would be ERROR_VALUE
 * @return position of the bin
 */
int getBinPosition(float value, const float LIMIT_MIN, const float LIMIT_MAX, const int N_BINS, bool exact) {
    int pos= int((value-LIMIT_MIN)*N_BINS/(LIMIT_MAX-LIMIT_MIN));
    if(!exact) {
        if(pos < 0)
            pos= 0;
        if(pos >= N_BINS)
            pos= N_BINS-1;
    } else
        if(pos < 0 || pos >= N_BINS)
            pos= ERROR_VALUE;
    return pos;
}

/**
 * Function that returns the position in a set of bins of a value (for a histogram)
 * @param value The value to analyze
 * @param LIMIT_MIN The value of the first bin
 * @param LIMIT_MAX The value of the last bin
 * @param N_BINS The number of bins
 * @return position of the bin
 */
int getBinPosition(float value, const float LIMIT_MIN, const float LIMIT_MAX, const int N_BINS) {
	return int((value-LIMIT_MIN)*N_BINS/(LIMIT_MAX-LIMIT_MIN));
}


#ifndef MACRO_H
#define MACRO_H

#define M_Velocity 30.0 // Unit: KM/H => Kilometer per hour
#define PER_MINUTE 1.0

#define M  1000.0 // TODO: Adapt this big-M constant !!!!! It is too big...
#define epz  0.0001
#define Rconst 6370
#define Positive 1
#define Negative -1
#define M_EPSILON 2.0
#define MAIN_STEP 10000
//Rejection
#define RMAX 10
#define WMAX 10
#define M_RejectionCost 1
#define GENERATOR_NO 10000
#define UNUSED(expr) (void)(expr);

#include <unordered_set>

#define M_1 0.8
#define M_2 -0.8

enum REQUESTTYPE
{
    PD = 1,
    PND = 2,
    NPD = 3,
    NPND = 4
};

#endif // MACRO_H

// Number of particles
#define N 5000

typedef struct {
    double x[N];
    double y[N];
    double z[N];
} Vector3NArray;

// Typedefs for specific physical attributes
typedef Vector3NArray Position;
typedef Vector3NArray Velocity;
typedef Vector3NArray Acceleration;
typedef Vector3NArray Force;
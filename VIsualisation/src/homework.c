#include "fem.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
void chk(int ierr) {
    if (ierr != 0) {
        printf("Error %d\n",ierr);
        exit(1);
    }
}
// double hermiteInterpolation(double d, double dStar, double h0, double hStar) {
//     if (d <= 0) return h0;
//     if (d >= dStar) return hStar;
//     double t = d / dStar;
//     double t2 = t * t;
//     double t3 = t2 * t;
//     return h0 * (2 * t3 - 3 * t2 + 1) + hStar * (3 * t2 - 2 * t3);
// }

double geoSize(double x, double y) {
    femGeo* theGeometry = geoGetGeometry();
    
    double h_max = theGeometry->h;       
    double h_min = h_max * 0.1;  // Ajuste ce facteur pour un effet plus visible
    double y_min = -theGeometry->MiddleY + theGeometry->MiddleY ;  
    double y_max = theGeometry->MiddleY  * theGeometry->hexRadius;  

    // Éviter une division par zéro
    if (y_max == y_min) return h_max;

    // Progression linéaire de la taille de maille entre h_min et h_max
    double t = (y - y_min) / (y_max - y_min);
    t = fmax(0.0, fmin(1.0, t)); // S'assurer que t reste dans [0,1]
    
    double h = h_min + (h_max - h_min) * t; // Progression linéaire
    
    return h;


}






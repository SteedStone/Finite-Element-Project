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
    double h_min;
    double h_max;
    double y_min;
    double y_max;
    if (theGeometry->hexa_triangles == 1) {
        double hexRadius = theGeometry->hexRadius;
        double numHexX = theGeometry->NumberOfHexagonsInX;
        double numHexY = theGeometry->NumberOfHexagonsInY;
        
        double distance = hexRadius * 1.1; // Distance entre les lignes
        h_max = theGeometry->h;       
        h_min = h_max * 0.2;  // Ajuste ce facteur pour un effet plus visible
        y_min =0;
        y_max = numHexY * distance - hexRadius;  
    } else {
        double hexRadius = theGeometry->hexRadius;
        double numHexX = theGeometry->NumberOfHexagonsInX;
        double numHexY = theGeometry->NumberOfHexagonsInY;
        h_max = theGeometry->h;       
        h_min = h_max * 0.1;  // Ajuste ce facteur pour un effet plus visible
        y_min =-hexRadius*sqrt(3)/2 - hexRadius/3;  
        y_max = numHexY * sqrt(3) * hexRadius+ hexRadius/3;  
    }

    
    
    
    
    // Éviter une division par zéro

    if (y_max == y_min) {
        return h_max;
    }

    // Progression linéaire de la taille de maille entre h_min et h_max
    double t = (y - y_min) / (y_max - y_min);
    t = fmax(0.0, fmin(1.0, t)); // S'assurer que t reste dans [0,1]
    
    // double h = h_min + (h_max - h_min) * pow(1 - t, 3); 
    double h = h_min + (h_max - h_min) * (1 - t);      
    //  printf("geoSize called: x=%.2f, y=%.2f -> h=%.5f\n", x, y, h); // 🔴 DEBUG
    return h;


}






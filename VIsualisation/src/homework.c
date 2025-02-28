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
double hermiteInterpolation(double d, double dStar, double h0, double hStar) {
    if (d <= 0) return h0;
    if (d >= dStar) return hStar;
    double t = d / dStar;
    double t2 = t * t;
    double t3 = t2 * t;
    return h0 * (2 * t3 - 3 * t2 + 1) + hStar * (3 * t2 - 2 * t3);
}

double geoSize(double x, double y){

    femGeo* theGeometry = geoGetGeometry();
    
    double h = theGeometry->h;
    double x0 = theGeometry->xNotch;
    double y0 = theGeometry->yNotch;
    double r0 = theGeometry->rNotch;
    double h0 = theGeometry->hNotch;
    double d0 = theGeometry->dNotch;
  
    
    double x1 = theGeometry->xHole;
    double y1 = theGeometry->yHole;
    double r1 = theGeometry->rHole;
    double h1 = theGeometry->hHole;
    double d1 = theGeometry->dHole;


//
//     A modifier !
//     
// Your contribution starts here ....
//
    
    double distNotch = sqrt((x - x0) * (x - x0) + (y - y0) * (y - y0)) - r0;
    double distHole = sqrt((x - x1) * (x - x1) + (y - y1) * (y - y1)) - r1;

    // Interpolation Hermite pour chaque zone
    double sizeNotch = hermiteInterpolation(distNotch, d0, h0, h);
    double sizeHole = hermiteInterpolation(distHole, d1, h1, h);

    // Prendre la taille minimale pour un meilleur raffinement
    return fmin(sizeNotch, sizeHole);
    
//   
// Your contribution ends here :-)
//

}


#define ___ 0

void geoMeshGenerate() {


    
    trianglePlot(0,0) ;
    HexagonPlot(0,0);
//
//  -1- Construction de la g�om�trie avec OpenCascade
//      On cr�e le rectangle
//      On cr�e les deux cercles
//      On soustrait les cercles du rectangle :-)
//
    
   

        
    // gmshModelOccSynchronize(&ierr);
    // gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    // gmshModelMeshGenerate(2, &ierr);
//
//  -2- D�finition de la fonction callback pour la taille de r�f�rence
//      Synchronisation de OpenCascade avec gmsh
//      G�n�ration du maillage (avec l'option Mesh.SaveAll :-)
                  
   
       
//
//  Generation de quads :-)
//
//    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
//    gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
//    gmshOptionSetNumber("Mesh.Algorithm", 8, &ierr);  chk(ierr);
//    gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);  chk(ierr);
//    gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  chk(ierr);
//    gmshModelMeshGenerate(2, &ierr);  
   
 
// //
// //  Plot of Fltk
// //
//   gmshFltkInitialize(&ierr);
//   gmshFltkRun(&ierr);  chk(ierr);

    
}




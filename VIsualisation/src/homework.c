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


    
 
//
//  -1- Construction de la g�om�trie avec OpenCascade
//      On cr�e le rectangle
//      On cr�e les deux cercles
//      On soustrait les cercles du rectangle :-)
//
    
    femGeo* theGeometry = geoGetGeometry();

    double hexRadius = theGeometry->hexRadius;  // Rayon de l'hexagone
    int numHexX = 10;  // Nombre d'hexagones en largeur
    int numHexY = 10;  // Nombre d'hexagones en hauteur
    int ierr;
    double meshSize = 0.1;

    int mainWireTags[100]; // Contient tous les hexagones extérieurs
    int innerWireTags[100]; // Contient les petits hexagones à soustraire
    int wireCount = 0, innerWireCount = 0;

    for (int i = 0; i < numHexX; i++) {
        for (int j = 0; j < numHexY; j++) {
            double x = i * 1.5 * hexRadius;
            double y = j * sqrt(3) * hexRadius + (i % 2) * sqrt(3) * hexRadius / 2;

            int points[6], innerPoints[6];
            for (int k = 0; k < 6; k++) {
                double angle = M_PI / 3 * k;
                double px = x + hexRadius * cos(angle);
                double py = y + hexRadius * sin(angle);
                points[k] = gmshModelOccAddPoint(px, py, 0, meshSize, -1, &ierr);
                ErrorGmsh(ierr);

                double innerPx = x + (hexRadius * 0.8) * cos(angle);
                double innerPy = y + (hexRadius * 0.8) * sin(angle);
                innerPoints[k] = gmshModelOccAddPoint(innerPx, innerPy, 0, meshSize, -1, &ierr);
                ErrorGmsh(ierr);
            }

            int lines[6], innerLines[6];
            for (int k = 0; k < 6; k++) {
                lines[k] = gmshModelOccAddLine(points[k], points[(k + 1) % 6], -1, &ierr);
                ErrorGmsh(ierr);
                innerLines[k] = gmshModelOccAddLine(innerPoints[k], innerPoints[(k + 1) % 6], -1, &ierr);
                ErrorGmsh(ierr);
            }

            int wire = gmshModelOccAddWire(lines, 6, -1, 1, &ierr);
            ErrorGmsh(ierr);
            int innerWire = gmshModelOccAddWire(innerLines, 6, -1, 1, &ierr);
            ErrorGmsh(ierr);

            mainWireTags[wireCount++] = wire;
            innerWireTags[innerWireCount++] = innerWire;
        }
    }

    // Création d'un grand rectangle englobant la structure
    int rectPoints[4];
    rectPoints[0] = gmshModelOccAddPoint(-hexRadius, -hexRadius, 0, meshSize, -1, &ierr);
    rectPoints[1] = gmshModelOccAddPoint((numHexX -1)* 1.5* hexRadius + hexRadius, -hexRadius, 0, meshSize, -1, &ierr);
    rectPoints[2] = gmshModelOccAddPoint((numHexX -1 ) * 1.5 * hexRadius + hexRadius, numHexY * sqrt(3) * hexRadius, 0, meshSize, -1, &ierr);
    rectPoints[3] = gmshModelOccAddPoint(-hexRadius, numHexY * sqrt(3) * hexRadius, 0, meshSize, -1, &ierr);

    int rectLines[4];
    for (int k = 0; k < 4; k++) {
        rectLines[k] = gmshModelOccAddLine(rectPoints[k], rectPoints[(k + 1) % 4], -1, &ierr);
        ErrorGmsh(ierr);
    }

    int rectWire = gmshModelOccAddWire(rectLines, 4, -1, 1, &ierr);
    ErrorGmsh(ierr);

    // Création de la surface principale
    int plateSurface = gmshModelOccAddPlaneSurface(&rectWire, 1, -1, &ierr);
    ErrorGmsh(ierr);

    // Création des surfaces hexagonales
    int hexSurfaces[100];
    for (int i = 0; i < innerWireCount; i++) {
        hexSurfaces[i] = gmshModelOccAddPlaneSurface(&innerWireTags[i], 1, -1, &ierr);
        ErrorGmsh(ierr);
    }
    

    // Découper les hexagones de la plaque
    int notch[innerWireCount][2];
    for (int i = 0; i < innerWireCount; i++) {
        notch[i][0] = 2;
        notch[i][1] = hexSurfaces[i];
    }
    int plate[2] = {2, plateSurface};

    gmshModelOccCut(plate, 2, (int *)notch, 2*innerWireCount, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    ErrorGmsh(ierr);

    // Synchronisation après ajout des surfaces
    gmshModelOccSynchronize(&ierr);
    ErrorGmsh(ierr);

    // Génération du maillage
    gmshModelMeshGenerate(2, &ierr);
    ErrorGmsh(ierr);

    // Affichage dans Gmsh
    gmshFltkRun(&ierr);
    ErrorGmsh(ierr);



        
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




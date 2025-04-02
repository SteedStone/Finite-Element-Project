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
        h_min = h_max * 0.1;  // Ajuste ce facteur pour un effet plus visible
        y_min =0;
        y_max = numHexY * distance - hexRadius;  
    } else {
        double hexRadius = theGeometry->hexRadius;
        double numHexX = theGeometry->NumberOfHexagonsInX;
        double numHexY = theGeometry->NumberOfHexagonsInY;
        h_max = theGeometry->h;       
        h_min = h_max * 0.01;  // Ajuste ce facteur pour un effet plus visible
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
    
    double h = h_min + (h_max - h_min) * pow(1 - t, 3); 
    // double h = h_min + (h_max - h_min) * (1 - t);      
    //  printf("geoSize called: x=%.2f, y=%.2f -> h=%.5f\n", x, y, h); // 🔴 DEBUG
    return h;


}


void trianglePlot() {

    int ierr;

    // gmshModelAdd("TriangleModel",  &ierr); // Créer un modèle pour les triangles
    // ErrorGmsh(ierr);


    femGeo* theGeometry = geoGetGeometry();
    double hexRadius = theGeometry->hexRadius;
    int numHexX = theGeometry->NumberOfTrianglesInX;
    int numHexY = theGeometry->NumberOfTrianglesInY;

    // Définition des tailles de maillage variables
    double meshSizeMin = 0.07;  // Taille de maille en haut (fine)
    double meshSizeMax = 0.3;   // Taille de maille en bas (grossière)
    double distance = hexRadius * 1.1; // Distance entre les lignes

    int innerWireTags[1000];
    int innerWireCount = 0;

    for (int i = 0; i < numHexX; i++) {
        for (int j = 0; j < numHexY; j++) {
            double x = i * 1.3 * hexRadius;
            double y = j * distance - hexRadius;

            // Calcul dynamique du meshSize en fonction de la hauteur
            double meshSize = meshSizeMax - (meshSizeMax - meshSizeMin) * ((double)j / numHexY);

            int innerPoints[3];

            if (i % 2 == 0) {
                innerPoints[0] = gmshModelOccAddPoint(x, y, 0, meshSize, -1, &ierr);
                innerPoints[1] = gmshModelOccAddPoint(x + hexRadius, y + sqrt(3) * hexRadius / 2, 0, meshSize, -1, &ierr);
                innerPoints[2] = gmshModelOccAddPoint(x - hexRadius, y + sqrt(3) * hexRadius / 2, 0, meshSize, -1, &ierr);
            } else {
                innerPoints[0] = gmshModelOccAddPoint(x, y + sqrt(3) * hexRadius / 2, 0, meshSize, -1, &ierr);
                innerPoints[1] = gmshModelOccAddPoint(x + hexRadius, y, 0, meshSize, -1, &ierr);
                innerPoints[2] = gmshModelOccAddPoint(x - hexRadius, y, 0, meshSize, -1, &ierr);
            }

            ErrorGmsh(ierr);

            int innerLines[3];
            for (int k = 0; k < 3; k++) {
                innerLines[k] = gmshModelOccAddLine(innerPoints[k], innerPoints[(k + 1) % 3], -1, &ierr);
                ErrorGmsh(ierr);
            }

            int innerWire = gmshModelOccAddWire(innerLines, 3, -1, 1, &ierr);
            ErrorGmsh(ierr);

            innerWireTags[innerWireCount++] = innerWire;
        }
    }

    int rectPoints[4];
    rectPoints[0] = gmshModelOccAddPoint(0, 0, 0, meshSizeMax, -1, &ierr);
    rectPoints[1] = gmshModelOccAddPoint(numHexX * 1.3 * hexRadius - (1.5) * hexRadius, 0, 0, meshSizeMax, -1, &ierr);
    rectPoints[2] = gmshModelOccAddPoint(numHexX * 1.3 * hexRadius - (1.5) * hexRadius, numHexY * distance - hexRadius, 0, meshSizeMin, -1, &ierr);
    rectPoints[3] = gmshModelOccAddPoint(0, numHexY * distance - hexRadius, 0, meshSizeMin, -1, &ierr);

    int rectLines[4];
    for (int k = 0; k < 4; k++) {
        rectLines[k] = gmshModelOccAddLine(rectPoints[k], rectPoints[(k + 1) % 4], -1, &ierr);
        ErrorGmsh(ierr);
    }

    int rectWire = gmshModelOccAddWire(rectLines, 4, -1, 1, &ierr);
    ErrorGmsh(ierr);

    int plateSurface = gmshModelOccAddPlaneSurface(&rectWire, 1, -1, &ierr);
    ErrorGmsh(ierr);

    int triSurfaces[100];
    for (int i = 0; i < innerWireCount; i++) {
        triSurfaces[i] = gmshModelOccAddPlaneSurface(&innerWireTags[i], 1, -1, &ierr);
        ErrorGmsh(ierr);
    }

    int notch[innerWireCount][2];
    for (int i = 0; i < innerWireCount; i++) {
        notch[i][0] = 2;
        notch[i][1] = triSurfaces[i];
    }

    int plate[2] = {2, plateSurface};
    gmshModelOccCut(plate, 2, (int *)notch, 2 * innerWireCount, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    ErrorGmsh(ierr);

    gmshModelOccSynchronize(&ierr);
    ErrorGmsh(ierr);


    

    // gmshModelMeshGenerate(2, &ierr);
    // ErrorGmsh(ierr);

    // gmshFltkRun(&ierr);
    // ErrorGmsh(ierr);
}



void HexagonPlot(){

//
//  -1- Construction de la g�om�trie avec OpenCascade
//      On cr�e le rectangle
//      On cr�e les deux cercles
//      On soustrait les cercles du rectangle :-)
//

    int ierr;

    // gmshModelAdd("HexagonModel",  &ierr);
    // ErrorGmsh(ierr);
    
    femGeo* theGeometry = geoGetGeometry();

    double hexRadius = theGeometry->hexRadius;  // Rayon de l'hexagone
    int numHexX = theGeometry->NumberOfHexagonsInX;  // Nombre d'hexagones en largeur
    int numHexY = theGeometry->NumberOfHexagonsInY;  // Nombre d'hexagones en hauteur
    double meshSizeMin = 0.07;  // Taille de maille en haut (fine)
    double meshSizeMax = 0.3;   // Taille de maille en bas (grossière)
    double f = 0.81;

    int mainWireTags[1000]; // Contient tous les hexagones extérieurs
    int innerWireTags[1000]; // Contient les petits hexagones à soustraire
    int trianglecoteTags[1000]; // Contient les triangles à soustraire
    int quadrilatereTags[1000]; // Contient les quadrilateres à soustraire
    int losangeTags[1000]; // Contient les losanges à soustraire
    

    int wireCount = 0, innerWireCount = 0, triangleWireCount = 0, quadrilatereWireCount = 0, losangeWireCount = 0;

    for (int i = 0; i < numHexX; i++) {
        for (int j = 0; j < numHexY; j++) {
            double x = i * 1.5 * hexRadius;
            double y = j * sqrt(3) * hexRadius + (i % 2) * sqrt(3) * hexRadius / 2;
            int trianglecotePoints[3], trianglecoteLines[3], previoustrianglepoints[numHexY];
            double meshSize = meshSizeMax - (meshSizeMax - meshSizeMin) * ((double)j / numHexY);
            if (j == 0 && i == 0) {
                for (int k = 0; k < 2; k++)
                {
                    double angle = M_PI / 3 * (k+3);
                    double px = x + hexRadius * cos(angle);
                    double py = y + hexRadius * sin(angle);
                    trianglecotePoints[k] = gmshModelOccAddPoint(px, py, 0, meshSize, -1, &ierr);
                    ErrorGmsh(ierr);
                }
                trianglecotePoints[2] = gmshModelOccAddPoint(-hexRadius, -hexRadius*sqrt(3)/2, 0, meshSize, -1, &ierr);
                previoustrianglepoints[0] = gmshModelOccAddPoint(hexRadius * cos(M_PI), hexRadius*sin(M_PI), 0, meshSize, -1, &ierr);
                ErrorGmsh(ierr);
                for (int k = 0; k < 3; k++)
                {
                    trianglecoteLines[k] = gmshModelOccAddLine(trianglecotePoints[k], trianglecotePoints[(k+1)%3], -1, &ierr);
                    ErrorGmsh(ierr);
                }
                
            }
            if(i == 0 && j > 0 ) {
                for (int k = 0; k < 2; k++)
                {
                    double angle = M_PI / 3 * (k+3);
                    double px = x + hexRadius * cos(angle);
                    double py = y + hexRadius * sin(angle);
                    trianglecotePoints[k] = gmshModelOccAddPoint(px, py, 0, meshSize, -1, &ierr);
                    ErrorGmsh(ierr);
                }
                trianglecotePoints[2] = previoustrianglepoints[j-1];
                previoustrianglepoints[j] = gmshModelOccAddPoint(hexRadius * cos(M_PI) + x, hexRadius*sin(M_PI)+ y, 0, meshSize, -1, &ierr);
                ErrorGmsh(ierr);
                for (int k = 0; k < 3; k++)
                {
                    trianglecoteLines[k] = gmshModelOccAddLine(trianglecotePoints[k], trianglecotePoints[(k+1)%3], -1, &ierr);
                    ErrorGmsh(ierr);
                }
            }
            if(i == numHexX-1 && j == 0) {
                for (int k = 0; k < 2; k++)
                {
                    double angle = -M_PI / 3 * (k);
                    double px = x + hexRadius * cos(angle);
                    double py = y + hexRadius * sin(angle);
                    trianglecotePoints[k] = gmshModelOccAddPoint(px, py, 0, meshSize, -1, &ierr);
                    ErrorGmsh(ierr);
                }
                trianglecotePoints[2] = gmshModelOccAddPoint((numHexX -1)* 1.5* hexRadius + hexRadius, -hexRadius*sqrt(3)/2, 0, meshSize, -1, &ierr);
                previoustrianglepoints[0] = gmshModelOccAddPoint(hexRadius * cos(0)+ x, hexRadius*sin(0), 0, meshSize, -1, &ierr);
                ErrorGmsh(ierr);
                for (int k = 0; k < 3; k++)
                {
                    trianglecoteLines[k] = gmshModelOccAddLine(trianglecotePoints[k], trianglecotePoints[(k+1)%3], -1, &ierr);
                    ErrorGmsh(ierr);
                }
            }
            if(i == numHexX-1 && j > 0 ) {
                for (int k = 0; k < 2; k++)
                {
                    double angle = -M_PI / 3 * (k);
                    double px = x + hexRadius * cos(angle);
                    double py = y + hexRadius * sin(angle);
                    trianglecotePoints[k] = gmshModelOccAddPoint(px, py, 0, meshSize, -1, &ierr);
                    ErrorGmsh(ierr);
                }
                trianglecotePoints[2] = previoustrianglepoints[j-1];
                previoustrianglepoints[j] = gmshModelOccAddPoint(hexRadius * cos(0) + x, hexRadius*sin(0)+ y, 0, meshSize, -1, &ierr);
                ErrorGmsh(ierr);
                for (int k = 0; k < 3; k++)
                {
                    trianglecoteLines[k] = gmshModelOccAddLine(trianglecotePoints[k], trianglecotePoints[(k+1)%3], -1, &ierr);
                    ErrorGmsh(ierr);
                }
            }
            int quadrilaterePoints[6], quadrilatereLines[6];
            if( i == 0 && j == numHexY - 1) {
                for (int k = 0; k < 3; k++)
                {
                    double angle = M_PI / 3 * (k+1);
                    double px = x + hexRadius * cos(angle);
                    double py = y + hexRadius * sin(angle);
                    quadrilaterePoints[k] = gmshModelOccAddPoint(px, py, 0, meshSize, -1, &ierr);
                    ErrorGmsh(ierr);
                }
                quadrilaterePoints[3] = gmshModelOccAddPoint(-hexRadius, numHexY * sqrt(3) * hexRadius+ hexRadius/3, 0, meshSize, -1, &ierr);
                quadrilaterePoints[5] = gmshModelOccAddPoint(1 * 1.5 * hexRadius + hexRadius * cos(M_PI/3 * 2), j * sqrt(3) * hexRadius + (1 % 2) * sqrt(3) * hexRadius / 2 + hexRadius * sin(M_PI/3 * 2), 0, meshSize, -1, &ierr);
                quadrilaterePoints[4] = gmshModelOccAddPoint(1 * 1.5 * hexRadius + hexRadius * cos(M_PI/3 * 2), numHexY * sqrt(3) * hexRadius+ hexRadius/3, 0, meshSize, -1, &ierr);

                ErrorGmsh(ierr);
                for (int k = 0; k < 6; k++)
                {
                    quadrilatereLines[k] = gmshModelOccAddLine(quadrilaterePoints[k], quadrilaterePoints[(k+1)%6], -1, &ierr);
                    ErrorGmsh(ierr);
                }
            }
            if( i == numHexX - 1 && j == numHexY - 1 ){
                for (int k = 0; k < 3; k++)
                {
                    double angle = M_PI / 3 * (k);
                    double px = x + hexRadius * cos(angle);
                    double py = y + hexRadius * sin(angle);
                    quadrilaterePoints[k] = gmshModelOccAddPoint(px, py, 0, meshSize, -1, &ierr);
                    ErrorGmsh(ierr);
                }

                quadrilaterePoints[5] = gmshModelOccAddPoint((numHexX -1 ) * 1.5 * hexRadius + hexRadius, numHexY * sqrt(3) * hexRadius+ hexRadius/3, 0, meshSize, -1, &ierr);
                quadrilaterePoints[3] = gmshModelOccAddPoint((numHexX - 2) * 1.5 * hexRadius + hexRadius * cos(M_PI/3 ), j * sqrt(3) * hexRadius + ((numHexX - 2) % 2) * sqrt(3) * hexRadius / 2 + hexRadius * sin(M_PI/3 ), 0, meshSize, -1, &ierr);
                quadrilaterePoints[4] = gmshModelOccAddPoint((numHexX - 2) * 1.5 * hexRadius + hexRadius * cos(M_PI/3 ), numHexY * sqrt(3) * hexRadius+ hexRadius/3, 0, meshSize, -1, &ierr);

                ErrorGmsh(ierr);
                for (int k = 0; k < 6; k++)
                {
                    quadrilatereLines[k] = gmshModelOccAddLine(quadrilaterePoints[k], quadrilaterePoints[(k+1)%6], -1, &ierr);
                    ErrorGmsh(ierr);
                }
            }
            int losangePoints[4], losangeLines[4];
            if(i % 2 != 0 && j == 0) {
                for (int k = 1; k < 3; k++) {
                    double angle = -M_PI / 3 * (k);
                    double px = x + hexRadius * cos(angle);
                    double py = y + hexRadius * sin(angle);
                    losangePoints[k] = gmshModelOccAddPoint(px, py, 0, meshSize, -1, &ierr);
                    ErrorGmsh(ierr);
                }
                losangePoints[3] = gmshModelOccAddPoint((i - 1) * 1.5 * hexRadius + hexRadius * cos(-M_PI/3 ), j * sqrt(3) * hexRadius + ((i- 1) % 2) * sqrt(3) * hexRadius / 2 + hexRadius * sin(-M_PI/3 ), 0, meshSize, -1, &ierr);
                losangePoints[0] = gmshModelOccAddPoint((i + 1) * 1.5 * hexRadius + hexRadius * cos(-M_PI/3 *2), j * sqrt(3) * hexRadius + ((i+ 1) % 2) * sqrt(3) * hexRadius / 2 + hexRadius * sin(-M_PI/3 *2), 0, meshSize, -1, &ierr);
                for (int k = 0; k < 4; k++) {
                    losangeLines[k] = gmshModelOccAddLine(losangePoints[k], losangePoints[(k + 1) % 4], -1, &ierr);
                    ErrorGmsh(ierr);
                }
                int losangeWire = gmshModelOccAddWire(losangeLines, 4, -1, 1, &ierr);
                ErrorGmsh(ierr);
                losangeTags[losangeWireCount++] = losangeWire;
            }
            if( i % 2 == 0 && j == numHexY - 1 && i != 0 && i != numHexX - 1) {
                for (int k = 1; k < 3; k++) {
                    double angle = M_PI / 3 * (k);
                    double px = x + hexRadius * cos(angle);
                    double py = y + hexRadius * sin(angle);
                    losangePoints[k] = gmshModelOccAddPoint(px, py, 0, meshSize, -1, &ierr);
                    ErrorGmsh(ierr);
                }
                losangePoints[3] = gmshModelOccAddPoint((i - 1) * 1.5 * hexRadius + hexRadius * cos(M_PI/3 ), j * sqrt(3) * hexRadius + ((i- 1) % 2) * sqrt(3) * hexRadius / 2 + hexRadius * sin(M_PI/3 ), 0, meshSize, -1, &ierr);
                losangePoints[0] = gmshModelOccAddPoint((i + 1) * 1.5 * hexRadius + hexRadius * cos(M_PI/3 *2), j * sqrt(3) * hexRadius + ((i+ 1) % 2) * sqrt(3) * hexRadius / 2 + hexRadius * sin(M_PI/3 *2), 0, meshSize, -1, &ierr);
                for (int k = 0; k < 4; k++) {
                    losangeLines[k] = gmshModelOccAddLine(losangePoints[k], losangePoints[(k + 1) % 4], -1, &ierr);
                    ErrorGmsh(ierr);
                }
                int losangeWire = gmshModelOccAddWire(losangeLines, 4, -1, 1, &ierr);
                ErrorGmsh(ierr);
                losangeTags[losangeWireCount++] = losangeWire;
            }
        
           
            // triangle interne 
            int points[6], innerPoints[6];
            for (int k = 0; k < 6; k++) {
                double angle = M_PI / 3 * k;
                double px = x + hexRadius * cos(angle);
                double py = y + hexRadius * sin(angle);
                points[k] = gmshModelOccAddPoint(px, py, 0, meshSize, -1, &ierr);
                ErrorGmsh(ierr);

                double innerPx = x + (hexRadius * f) * cos(angle);
                double innerPy = y + (hexRadius * f) * sin(angle);
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
            if (j == 0 && i==0) {
                int trianglecoteWire = gmshModelOccAddWire(trianglecoteLines, 3, -1, 1, &ierr);
                ErrorGmsh(ierr);
                trianglecoteTags[triangleWireCount++] = trianglecoteWire;
            }
            if(i == 0 && j > 0 ) {
                int trianglecoteWire = gmshModelOccAddWire(trianglecoteLines, 3, -1, 1, &ierr);
                ErrorGmsh(ierr);
                trianglecoteTags[triangleWireCount++] = trianglecoteWire;
            }
            if(i == numHexX-1 && j == 0) {
                int trianglecoteWire = gmshModelOccAddWire(trianglecoteLines, 3, -1, 1, &ierr);
                ErrorGmsh(ierr);
                trianglecoteTags[triangleWireCount++] = trianglecoteWire;
            }
            if(i == numHexX-1 && j > 0) {
                int trianglecoteWire = gmshModelOccAddWire(trianglecoteLines, 3, -1, 1, &ierr);
                ErrorGmsh(ierr);
                trianglecoteTags[triangleWireCount++] = trianglecoteWire;
            }
            if( i == 0 && j == numHexY - 1) {
                int quadrilatereWire = gmshModelOccAddWire(quadrilatereLines, 6, -1, 1, &ierr);
                ErrorGmsh(ierr);
                quadrilatereTags[0] = quadrilatereWire;
            }
            if( i == numHexX - 1 && j == numHexY - 1) {
                int quadrilatereWire = gmshModelOccAddWire(quadrilatereLines, 6, -1, 1, &ierr);
                ErrorGmsh(ierr);
                quadrilatereTags[1] = quadrilatereWire;
            }
            
            
        }
    }

    // Création d'un grand rectangle englobant la structure
    int rectPoints[4];
    rectPoints[0] = gmshModelOccAddPoint(-hexRadius, -hexRadius*sqrt(3)/2 - hexRadius/3, 0, meshSizeMax, -1, &ierr);
    rectPoints[1] = gmshModelOccAddPoint((numHexX -1)* 1.5* hexRadius + hexRadius, -hexRadius*sqrt(3)/2 - hexRadius/3, 0, meshSizeMax, -1, &ierr);
    rectPoints[2] = gmshModelOccAddPoint((numHexX -1)* 1.5* hexRadius + hexRadius, numHexY * sqrt(3) * hexRadius+ hexRadius/3, 0, meshSizeMin, -1, &ierr);
    rectPoints[3] = gmshModelOccAddPoint(-hexRadius, numHexY * sqrt(3) * hexRadius+ hexRadius/3, 0, meshSizeMin, -1, &ierr);

    // int newRectPoints[4];
    // newRectPoints[0] = gmshModelOccAddPoint(-hexRadius, -hexRadius*sqrt(3)/2 , 0, meshSizeMax, -1, &ierr);
    // newRectPoints[1] = gmshModelOccAddPoint((numHexX -1)* 1.5* hexRadius + hexRadius, -hexRadius*sqrt(3)/2 , 0, meshSizeMax, -1, &ierr);
    // newRectPoints[2] = gmshModelOccAddPoint((numHexX -1 ) * 1.5 * hexRadius + hexRadius, numHexY * sqrt(3) * hexRadius, 0, meshSizeMin, -1, &ierr);
    // newRectPoints[3] = gmshModelOccAddPoint(-hexRadius, numHexY * sqrt(3) * hexRadius, 0, meshSizeMin, -1, &ierr);
    

    int rectLines[4];
    for (int k = 0; k < 4; k++) {
        rectLines[k] = gmshModelOccAddLine(rectPoints[k], rectPoints[(k + 1) % 4], -1, &ierr);
        ErrorGmsh(ierr);
    }

    // int newRectLines[4];
    // for (int k = 0; k < 4; k++) {
    //     newRectLines[k] = gmshModelOccAddLine(newRectPoints[k], newRectPoints[(k + 1) % 4], -1, &ierr);
    //     ErrorGmsh(ierr);
    // }


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

    // Création des surfaces triangulaires
    int triangleSurfaces[100];
    for (int i = 0; i < triangleWireCount; i++) {
        triangleSurfaces[i] = gmshModelOccAddPlaneSurface(&trianglecoteTags[i], 1, -1, &ierr);
        ErrorGmsh(ierr);


    }
    // Creation des surfaces quadrilateres
    int quadrilatereSurfaces[2];
    for (int i = 0; i < 2; i++) {
        quadrilatereSurfaces[i] = gmshModelOccAddPlaneSurface(&quadrilatereTags[i], 1, -1, &ierr);
        ErrorGmsh(ierr);
    }
    int plate[2] = {2, plateSurface};
    // // Découper les quadrilateres
    int notchquadrilatere[2][2];
    for (int i = 0; i < 2; i++) {
        notchquadrilatere[i][0] = 2 ;
        notchquadrilatere[i][1] = quadrilatereSurfaces[i] ;
    }
    gmshModelOccCut(plate, 2, (int *)notchquadrilatere, 2*2  ,NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccSynchronize(&ierr);

    // Creation des surfaces losanges
    int losangeSurfaces[100];
    for (int i = 0; i < losangeWireCount; i++) {
        losangeSurfaces[i] = gmshModelOccAddPlaneSurface(&losangeTags[i], 1, -1, &ierr);
        ErrorGmsh(ierr);
    }
    // Découper les losanges
    int notchlosange[losangeWireCount][2];
    for (int i = 0; i < losangeWireCount; i++) {
        notchlosange[i][0] = 2 ;
        notchlosange[i][1] = losangeSurfaces[i] ;
    }
    gmshModelOccCut(plate, 2, (int *)notchlosange, 2*losangeWireCount  ,NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccSynchronize(&ierr);

    // Découper les triangles 
    int notchtriangle[triangleWireCount][2];
    for (int i = 0; i < triangleWireCount; i++) {
        notchtriangle[i][0] = 2 ;
        notchtriangle[i][1] = triangleSurfaces[i] ;
    }
    gmshModelOccCut(plate, 2, (int *)notchtriangle, 2*triangleWireCount  ,NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccSynchronize(&ierr);

    // Découper les hexagones de la plaque
    int notch[innerWireCount][2];
    for (int i = 0; i < innerWireCount; i++) {
        notch[i][0] = 2;
        notch[i][1] = hexSurfaces[i];
    }

    gmshModelOccCut(plate, 2, (int *)notch, 2*innerWireCount, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    ErrorGmsh(ierr);

    // Synchronisation apres ajout des surfaces
    gmshModelOccSynchronize(&ierr);
    ErrorGmsh(ierr);

  
    // // Affichage dans Gmsh
    // gmshFltkRun(&ierr);
    // ErrorGmsh(ierr);



}

void femFindBoundaryNodes(femGeo *theProblem, double targetY, double epsilon , char *name)
{
    femMesh *theEdges = theProblem->theEdges;

    int size = theProblem->theElements->nElem;
    int *map = malloc(size * sizeof(int));
    if (map == NULL) { Error("Memory allocation failed !"); exit(EXIT_FAILURE); return; }
    for (int i = 0; i < size; i++) { map[i] = 0; }

    int nDomains = theProblem->nDomains;
    for (int iDomain = 0; iDomain < nDomains; iDomain++)
    {
        femDomain *theDomain = theProblem->theDomains[iDomain];
        for (int i = 0; i < theDomain->nElem; i++)
        {
           int iEdge = theDomain->elem[i];
            for (int j = 0; j < 2; j++)
            {
               int iNode = theEdges->elem[iEdge * 2 + j];
               map[iNode] = 1;
            }
        }
    }

    // Compter les nœuds de la frontière qui sont à targetY
    int nFilteredBoundary = 0;
    for (int i = 0; i < size; i++)
    {
        
        if (map[i] == 1 && fabs(theProblem->theNodes->Y[i] - targetY) < epsilon) { 
            nFilteredBoundary++; 
        }
    }

    // Création du domaine pour les nœuds filtrés
    femDomain *theFilteredBoundary = malloc(sizeof(femDomain));
    if (theFilteredBoundary == NULL) { Error("Memory allocation failed !"); exit(EXIT_FAILURE); return; }
    theProblem->nDomains++;
    theProblem->theDomains = realloc(theProblem->theDomains, theProblem->nDomains * sizeof(femDomain *));
    if (theProblem->theDomains == NULL) { Error("Memory allocation failed !"); exit(EXIT_FAILURE); return; }
    theProblem->theDomains[theProblem->nDomains - 1] = theFilteredBoundary;
    theFilteredBoundary->nElem = nFilteredBoundary;
    theFilteredBoundary->elem = malloc(nFilteredBoundary * sizeof(int));
    if (theFilteredBoundary->elem == NULL) { Error("Memory allocation failed !"); exit(EXIT_FAILURE); return; }

    theFilteredBoundary->mesh = NULL;
    sprintf(theFilteredBoundary->name, name);

    // Stocker les nœuds filtrés
    int index = 0;
    for (int i = 0; i < size; i++)
    {
        if (map[i] == 1 && fabs(theProblem->theNodes->Y[i] - targetY) < epsilon ) {
            // printf("TArgetY: %f\n", targetY);
            // printf("Node %d: %f\n", i, theProblem->theNodes->Y[i]);
            theFilteredBoundary->elem[index++] = i; 
        }
    }
    
    free(map);
}



void femSolverSet(femSolver *mySolver, double **newA, double *newB)
{
    if (mySolver->type == FEM_FULL)
    {
        femFullSystem *mySystem = mySolver->solver;
        mySystem->A = newA;
        mySystem->B = newB;
    }
    else if (mySolver->type == FEM_BAND)
    {
        femBandSystem *mySystem = mySolver->solver;
        mySystem->A = newA;
        mySystem->B = newB;
    }
    else { Error("Unexpected solver type"); }
}
void femElasticitySigma(femProblem *theProblem, double *sigmaXX, double *sigmaYY, double *sigmaXY)
{
    femIntegration *theRule  = theProblem->rule;
    femDiscrete *theSpace    = theProblem->space;
    femGeo *theGeometry = theProblem->geometry;
    femNodes *theNodes       = theGeometry->theNodes;
    femMesh *theMesh         = theGeometry->theElements;

    int nLocal = theSpace->n;

    femSolver *theSolver;
    double **M1, **M2, **M3, *B1, *B2, *B3, *epsilonXX, *epsilonYY, *epsilonXY, *theSoluce;
    double x[nLocal], y[nLocal], phi[nLocal], dphidxsi[nLocal], dphideta[nLocal], dphidx[nLocal], dphidy[nLocal];
    double u[nLocal], v[nLocal], xsi, eta, weight, jac, weightedJac, dxdxsi, dxdeta, dydxsi, dydeta, a, b, c;
    int iElem, iInteg, iEdge, i, j, map[nLocal], nNodes;

    nNodes    = theNodes->nNodes;
    int* number = theNodes->number;
    theSoluce = theProblem->soluce;    
    theSolver = femSolverFullCreate(nNodes , theProblem->sizeLoc);

    if (theSolver->type == FEM_FULL) {
        M1 = ((femFullSystem *) theSolver->solver)->A  ;
        B1 = ((femFullSystem *) theSolver->solver)->B ;
    } else if(theSolver->type == FEM_BAND) 
    {
        M1 = ((femBandSystem *) theSolver->solver)->A  ;
        B1 = ((femBandSystem *) theSolver->solver)->B ;
    } 
    
    

    B2 = (double *) malloc(nNodes * sizeof(double));
    B3 = (double *) malloc(nNodes * sizeof(double));

    if (B2 == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return; }
    if (B3 == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return; }

    for (i = 0; i < nNodes; i++) { B2[i] = 0.0; B3[i] = 0.0; }
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++)
    {
        for (i = 0; i < theSpace->n; i++)
        {
            map[i] = theMesh->elem[iElem * nLocal + i];
            x[i] = theNodes->X[map[i]];
            y[i] = theNodes->Y[map[i]];
            u[i] = theSoluce[2 * map[number[i]]];
            v[i] = theSoluce[2 * map[number[i]] + 1];

        }

        for (iInteg = 0; iInteg < theRule->n; iInteg++)
        {
            xsi    = theRule->xsi[iInteg];
            eta    = theRule->eta[iInteg];
            weight = theRule->weight[iInteg];

            femDiscretePhi2(theSpace, xsi, eta, phi);
            femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);

            dxdxsi = 0.0; dydxsi = 0.0;
            dxdeta = 0.0; dydeta = 0.0;
            for (i = 0; i < theSpace->n; i++)
            {
                dxdxsi += x[i] * dphidxsi[i];
                dxdeta += x[i] * dphideta[i];
                dydxsi += y[i] * dphidxsi[i];
                dydeta += y[i] * dphideta[i];
            }

            jac = dxdxsi * dydeta - dxdeta * dydxsi;
            if (jac < 0.0) { printf("Negative jacobian! Your mesh is oriented in reverse. The normals will be wrong\n"); }
            jac = fabs(jac);

            for (i = 0; i < theSpace->n; i++)
            {
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
            }

            weightedJac = jac * weight;

            for (i = 0; i < theSpace->n; i++)
            {
                for (j = 0; j < theSpace->n; j++) { M1[map[i]][map[j]] += phi[i] * phi[j] * weightedJac; }
                B1[map[i]] += phi[i] * dphidx[i] * u[i] * weightedJac;
                B2[map[i]] += phi[i] * dphidy[i] * v[i] * weightedJac;
                B3[map[i]] += phi[i] * (dphidy[i] * u[i] + dphidx[i] * v[i]) * weightedJac / 2.0;
            }
        }
    }

    M2 = (double **) malloc(sizeof(double *) * nNodes);
    M3 = (double **) malloc(sizeof(double *) * nNodes);

    if (M2 == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    if (M3 == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }

    for (i = 0; i < nNodes; i++)
    {
        M2[i] = (double *) malloc(nNodes * sizeof(double));
        M3[i] = (double *) malloc(nNodes * sizeof(double));

        if (M2[i] == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
        if (M3[i] == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }

        memcpy(M2[i], M1[i], nNodes * sizeof(double));
        memcpy(M3[i], M1[i], nNodes * sizeof(double));
    }

    epsilonXX = (double *) malloc(nNodes * sizeof(double));
    epsilonYY = (double *) malloc(nNodes * sizeof(double));
    epsilonXY = (double *) malloc(nNodes * sizeof(double));

    if (epsilonXX == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return; }
    if (epsilonYY == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return; }
    if (epsilonXY == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return; }

    if (theSolver->type == FEM_FULL) {

        femSolverEliminate(theSolver);
        memcpy(epsilonXX, ((femFullSystem *) theSolver->solver)->B, nNodes * sizeof(double));
        
        femSolverSet(theSolver, M2, B2);
        femSolverEliminate(theSolver);
        memcpy(epsilonYY, ((femFullSystem *) theSolver->solver)->B, nNodes * sizeof(double));

        femSolverSet(theSolver, M3, B3);
        femSolverEliminate(theSolver);
        memcpy(epsilonXY, ((femFullSystem *) theSolver->solver)->B, nNodes * sizeof(double));
    } else  if(theSolver->type == FEM_BAND ) {
        femSolverEliminate(theSolver);
        memcpy(epsilonXX, ((femBandSystem *) theSolver->solver)->B, nNodes * sizeof(double));

        femSolverSet(theSolver, M2, B2);
        femSolverEliminate(theSolver);
        memcpy(epsilonYY, ((femBandSystem *) theSolver->solver)->B, nNodes * sizeof(double));

        femSolverSet(theSolver, M3, B3);
        femSolverEliminate(theSolver);
        memcpy(epsilonXY, ((femBandSystem *) theSolver->solver)->B, nNodes * sizeof(double));
    }

    a = theProblem->A;
    b = theProblem->B;
    c = theProblem->C;

    for (i = 0; i < nNodes; i++)
    {
        sigmaXX[i] = a * epsilonXX[i] + b * epsilonYY[i];
        sigmaXY[i] = 2 * c * epsilonXY[i];
        sigmaYY[i] = b * epsilonXX[i] + a * epsilonYY[i];
    }

    free(epsilonXX); epsilonXX = NULL;
    free(epsilonYY); epsilonYY = NULL;
    free(epsilonXY); epsilonXY = NULL;
    for (i = 0; i < nNodes; i++) { free(M2[i]); M2[i] = NULL; free(M3[i]); M3[i] = NULL; }
    free(M2); M2 = NULL;
    free(M3); M3 = NULL;
    free(B2); B2 = NULL;
    free(B3); B3 = NULL;
}


 




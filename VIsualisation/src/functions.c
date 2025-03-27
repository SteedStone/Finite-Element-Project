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






/*
 *  main.c
 *  Library for EPL1110 : Finite Elements for dummies
 *  Utilisation de l'API de GMSH pour cr�er un maillage
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */
 
 #include "glfem.h"
 #define TRIANGLE 1 
 #define HEXAGON 2
 double fun(double x, double y) 
 {
     return 1;
 }
 int main(int argc, char *argv[]) {  
    if (argc < 5) {
        printf("Usage: %s <type> <show_visualisation>\n", argv[0]);
        printf("type: 1 (TRIANGLE), 2 (HEXAGON)\n");
        printf("show_visualisation: 2 (non), 1 (oui)\n");
        printf("Type de maillage: 0 (Triangle) , 1 (Quadrilatère)\n");
        printf("Grand problème: 0 (Non) , 1 (Oui)\n");
        return EXIT_FAILURE;
    }
    
    int type = atoi(argv[1]);
    int show_visualisation_maillage = atoi(argv[2]);
    int type1 = atoi(argv[3]);
    int simple = 0;
    int big = atoi(argv[4]);
     
 

  
    
       
     int ierr;
     femGeo* theGeometry ;
     if(!simple){
         geoInitialize();
         theGeometry = geoGetGeometry();
         if (big) {
         theGeometry->NumberOfHexagonsInX = 9;
         theGeometry->NumberOfHexagonsInY = 5;
         theGeometry->NumberOfTrianglesInX = 11;
         theGeometry->NumberOfTrianglesInY = 7;
         
         }
            else {
                theGeometry->NumberOfHexagonsInX = 3;
                theGeometry->NumberOfHexagonsInY = 2;
                theGeometry->NumberOfTrianglesInX = 5;
                theGeometry->NumberOfTrianglesInY = 3;
            }
         theGeometry->hexRadius = 0.01;
         theGeometry->MiddleX = (theGeometry->NumberOfHexagonsInX -1 ) * 1.5 * theGeometry->hexRadius + theGeometry->hexRadius - (-theGeometry->hexRadius) ; 
         theGeometry->MiddleY = theGeometry->NumberOfHexagonsInY * sqrt(3) * theGeometry->hexRadius - (-theGeometry->hexRadius*sqrt(3)/2) ;
         theGeometry->elementType = type1;
         theGeometry->hexa_triangles = type ;
         
         geoMeshGenerate();
         geoMeshImport();

         if (theGeometry->hexa_triangles == 1) {
            if (big) {
                geoSetDomainName(181 , "Top"); // 0
                geoSetDomainName(201 , "Bottom"); // 1
            } else {
                geoSetDomainName(25 , "Top"); 
                geoSetDomainName(33 , "Bottom"); 
            }

         
         } else {
            if (big) {
            geoSetDomainName(547 , "Top"); // hexagone (9;5) = 549 , Hexagone(3;3) = 149 115
            geoSetDomainName(564 , "Bottom"); // hexagone (9;5) = 564 , Hexagone(3;3) = 158 124
            } else {
                geoSetDomainName(79 , "Top"); // hexagone (3;3) = 149 , Hexagone(3;3) = 149 115
                geoSetDomainName(84 , "Bottom"); // hexagone (3;3) = 158 , Hexagone(3;3) = 158 124
            }
        }
         
     }
     else {
         double Lx = 1.0;
         double Ly = 1.0;
         
         geoInitialize();
         theGeometry = geoGetGeometry();
         
         theGeometry->LxPlate     =  Lx;
         theGeometry->LyPlate     =  Ly;     
         theGeometry->h           =  Lx * 0.075;    
         theGeometry->elementType = FEM_QUAD;
         theGeometry->hexa_triangles = 3;
 
         geoMeshGenerate();
         geoMeshImport();
         geoSetDomainName(0,"Symmetry");
         geoSetDomainName(7,"Bottom");
         geoSetDomainName(1,"Top");
 
         geoMeshWrite("data/elasticity.txt");
     }
     
     
     
 
     // EXemple plus simple 
     
 
     
 
 //
 //  -2- Creation du fichier du maillage
 //
     if(!simple){
     char filename[] = "data/mesh.txt";
     geoMeshWrite(filename);
     }
     
 
 //
 //  -3- Visualisation du maillage
 //  
     if(show_visualisation_maillage) {
         double *meshSizeField = malloc(theGeometry->theNodes->nNodes*sizeof(double));
         femNodes *theNodes = theGeometry->theNodes;
         for(int i=0; i < theNodes->nNodes; ++i)
             meshSizeField[i] = geoSize(theNodes->X[i], theNodes->Y[i]);
         double hMin = femMin(meshSizeField,theNodes->nNodes);  
         double hMax = femMax(meshSizeField,theNodes->nNodes);  
         printf(" ==== Global requested h : %14.7e \n",theGeometry->h);
         printf(" ==== Minimum h          : %14.7e \n",hMin);
         printf(" ==== Maximum h          : %14.7e \n",hMax);
 
 
         int mode = 1; // Change mode by pressing "j", "k", "l"
         int domain = 530;
         int freezingButton = FALSE;
         double t, told = 0;
         char theMessage[256];
         double pos[2] = {20,460};
 
 
         GLFWwindow* window = glfemInit("EPL1110 : Mesh generation ");
         glfwMakeContextCurrent(window);
 
         do {
             int w,h;
 
 
             glfwGetFramebufferSize(window,&w,&h);
             glfemReshapeWindows(theGeometry->theNodes,w,h);
 
             t = glfwGetTime();  
             if (glfwGetKey(window,'D') == GLFW_PRESS) { mode = 0;}
             if (glfwGetKey(window,'V') == GLFW_PRESS) { mode = 1;}
             if (glfwGetKey(window,'N') == GLFW_PRESS && freezingButton == FALSE) { domain++; freezingButton = TRUE; told = t;}
 
             
             if (t-told > 0.5) {freezingButton = FALSE; }
                 
             
             
                 
             if (mode == 1) {
                 glfemPlotField(theGeometry->theElements, meshSizeField);
                 glfemPlotMesh(theGeometry->theElements); 
                 sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
 
                 
                 glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
 
                 
                 
                 }
             if (mode == 0) {
                 domain = domain % theGeometry->nDomains;
                 glfemPlotDomain( theGeometry->theDomains[domain]); 
                 
                 
                 
                 sprintf(theMessage, "%s : %d ",theGeometry->theDomains[domain]->name,domain);
 
                 
                 glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
                 }
                 
                 glfwSwapBuffers(window);
                 glfwPollEvents();
         } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
                     glfwWindowShouldClose(window) != 1 );
                 
 
         free(meshSizeField);  
     }
             
     //
     //  -4- Creation probleme 
     //
      
     geoFinalize();
     glfwTerminate(); 
     
 
 
     
     exit(EXIT_SUCCESS);
     return 0;  
 }
 
  
 
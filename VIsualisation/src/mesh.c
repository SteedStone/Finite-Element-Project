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
    if (argc < 4) {
        printf("Usage: %s <type> <show_visualisation>\n", argv[0]);
        printf("type: 1 (TRIANGLE), 2 (HEXAGON)\n");
        printf("show_visualisation: 0 (non), 1 (oui)\n");
        printf("Type de maillage: 0 (Triangle) , 1 (Quadrilatère)\n");
        return EXIT_FAILURE;
    }
    
    int type = atoi(argv[1]);
    int show_visualisation_maillage = atoi(argv[2]);
    int type1 = atoi(argv[3]);
    int simple = 0;
     
 
 
     // On va devoir avoir un gros rectangle de lx = 2 et ly = 0.5
     // Puis deux petits cercles au bout de chaque côté du rectangle de rayon 0.1 et 0.2 
     // et enfin deux rectangle dans les cercles. 
  
    
       
     int ierr;
     femGeo* theGeometry ;
     if(!simple){
         geoInitialize();
         theGeometry = geoGetGeometry();
         theGeometry->h = 11;
         theGeometry->NumberOfHexagonsInX = 9;
         theGeometry->NumberOfHexagonsInY = 5;
         theGeometry->NumberOfTrianglesInX = 5;
         theGeometry->NumberOfTrianglesInY = 3;
         theGeometry->hexRadius = 9.0;
         theGeometry->MiddleX = (theGeometry->NumberOfHexagonsInX -1 ) * 1.5 * theGeometry->hexRadius + theGeometry->hexRadius - (-theGeometry->hexRadius) ; 
         theGeometry->MiddleY = theGeometry->NumberOfHexagonsInY * sqrt(3) * theGeometry->hexRadius - (-theGeometry->hexRadius*sqrt(3)/2) ;
         theGeometry->elementType = type1;
         theGeometry->hexa_triangles = type ;
         
         // Appel nos fonctions pour créer notre maillage avec soit les hexagones soit les triangles avec la librairie gsmh 
         geoMeshGenerate();
         // On remplit la structure TheGeometry avec les noeuds, edges ect depuis gsmh vers notre programme 
         geoMeshImport();
         // On trouve créer les différents domaines ou on va venir appliquer nos forces celon si on apllique sur les hexagones ou les triangles.
        //  if (theGeometry->hexa_triangles == 1) {
        //      femFindBoundaryNodes(theGeometry, 0 , 10e-6 , "Bottom"); // numéro 
        //      femFindBoundaryNodes(theGeometry, 90 , 10e-6 , "Top");
        //  } else {
        //     geoSetDomainName(117 , "Top");
        //     geoSetDomainName(124 , "Bottom");
        //  }
         if (theGeometry->hexa_triangles == 1) {

        geoSetDomainName(25 , "Top"); // 298
            geoSetDomainName(33 , "Bottom"); // 327
         
         } else {
            geoSetDomainName(549 , "Top"); // 549
            geoSetDomainName(564 , "Bottom"); // 564
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
     // On a rajouté une option pour voir si on veut l'afficher ou pas 
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
         int domain = 0;
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
         //    glfemChangeState(&mode, theMeshes->nMesh);
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
                 
         // Check if the ESC key was pressed or the window was closed
 
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
 
  
 
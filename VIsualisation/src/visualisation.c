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
 #define STB_IMAGE_WRITE_IMPLEMENTATION
 #include "stb_image_write.h"
 double fun(double x, double y) 
 {
     return 1;
 }
 
 
 
 void captureFrame(int frame, GLFWwindow* window) {
     int width, height;
     glfwGetFramebufferSize(window, &width, &height);
 
     // Créer un tableau pour stocker les pixels
     unsigned char* pixels = (unsigned char*)malloc(3 * width * height);
 
     printf("Framebuffer size: %d x %d\n", width, height);
 
     GLint viewport[4];
     glGetIntegerv(GL_VIEWPORT, viewport);
     printf("Viewport: x=%d, y=%d, width=%d, height=%d\n", viewport[0], viewport[1], viewport[2], viewport[3]);
 
     int winWidth, winHeight;
     glfwGetWindowSize(window, &winWidth, &winHeight);
     printf("Window size: %d x %d\n", winWidth, winHeight);
 
 
 
     // Lire le contenu du framebuffer (RGB)
     glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels);
 
 
 
     // Inverser l'image verticalement (OpenGL utilise une origine en bas à gauche)
     unsigned char* flippedPixels = (unsigned char*)malloc(3 * width * height);
     for (int y = 0; y < height; y++) {
         for (int x = 0; x < width; x++) {
             int srcIndex = 3 * (y * width + x);
             int destIndex = 3 * ((height - y - 1) * width + x);
             flippedPixels[destIndex] = pixels[srcIndex];
             flippedPixels[destIndex + 1] = pixels[srcIndex + 1];
             flippedPixels[destIndex + 2] = pixels[srcIndex + 2];
         }
     }
 
     // Sauvegarder l'image sous un fichier
     char filename[100];
     
     sprintf(filename, "frame/frame_%d.png", frame/10);
     stbi_write_png(filename, width, height, 3, flippedPixels, width * 3);
 
     // Libérer la mémoire
     free(pixels);
     free(flippedPixels);
 }
 
 int main(int argc, char *argv[])
 {  
     printf("\n\n    V : Mesh and size mesh field \n");
     printf("    D : Domains \n");
     printf("    N : Next domain highlighted\n");
 
     
 
     if (argc < 2) {
         fprintf(stderr, "Usage: %s <path_to_file>\n", argv[0]);
         return EXIT_FAILURE;
     }
 
     char *file_path = argv[1];
     
 
 
     // On va devoir avoir un gros rectangle de lx = 2 et ly = 0.5
     // Puis deux petits cercles au bout de chaque côté du rectangle de rayon 0.1 et 0.2 
     // et enfin deux rectangle dans les cercles. 
  
    
       
     int ierr;
     femGeo* theGeometry ;
     geoInitialize();
     theGeometry = geoGetGeometry();
     
     
     
 
 
     geoMeshRead(file_path);
     
 
     
     double E   = 211.e9;
     double nu  = 0.3;
     double rho = 7.85e3; 
     double g   = 9.81;
    
     // Initialisation du problème avec les conditions aux bords
     // ON remplit juste la structure theProblem avec les valeurs de E, nu, rho, g et le type de problème qu'on veut résoudre
     femProblem* theProblem = femElasticityCreate(theGeometry,E,nu,rho,g,PLANAR_STRESS , FEM_FULL , FEM_NO);
     
     femElasticityAddBoundaryCondition(theProblem,"Bottom",DIRICHLET_X,0.0);
     femElasticityAddBoundaryCondition(theProblem,"Bottom",DIRICHLET_Y,0.0);
 
     femElasticityAddBoundaryCondition(theProblem,"Top",NEUMANN_Y,-1e6);
     
     femElasticityPrint(theProblem);
     // femFullSystemPrint(theProblem->solver->solver);
     // printf("Taille system %d" , theProblem->size);
     
     //  -3- Resolution du probleme et calcul des forces
     
 
     double *theSoluce = femElasticitySolve(theProblem);
 
     double *theForces = femElasticityForces(theProblem);
     
     
     double area = femElasticityIntegrate(theProblem, fun);   
 
     //
     //  -4- Deformation du maillage pour le plot final
     //      Creation du champ de la norme du deplacement
     //
 
 
     femNodes *theNodes = theGeometry->theNodes;
     double *normDisplacement = malloc(theNodes->nNodes * sizeof(double));
     double *forcesX = malloc(theNodes->nNodes * sizeof(double));
     double *forcesY = malloc(theNodes->nNodes * sizeof(double));
 
 
     
    GLFWwindow* window = glfemInit("EPL1110 : Recovering forces on constrained nodes");
    glfwMakeContextCurrent(window);
 
 
 
     double deformationFactor_start = 0.0;  // Commence à 0 (pas de déformation initiale)
     double deformationFactor_end = 1000.0;  // Incrément du facteur de déformation à chaque itération
     double deformationFactor = 0.0;
     
     double *X0 = malloc(sizeof(double) * theNodes->nNodes);
     double *Y0 = malloc(sizeof(double) * theNodes->nNodes);
 
     for (int j = 0; j < theNodes->nNodes; j++) {
         X0[j] = theNodes->X[j];
         Y0[j] = theNodes->Y[j];
     }
 
 
     // Boucle d'animation
     for (double i = deformationFactor_start; i < deformationFactor_end; i += 10.0) {
         double deformationFactor = i;
     
         // Repartir systématiquement des coordonnées initiales
         for (int j = 0; j < theNodes->nNodes; j++) {
             theNodes->X[j] = X0[j] + deformationFactor * theSoluce[2*j + 0];
             theNodes->Y[j] = Y0[j] + deformationFactor * theSoluce[2*j + 1];
         }
         
         int w, h;
         glfwGetFramebufferSize(window, &w, &h);
         glfemReshapeWindows(theGeometry->theNodes, w, h);
     
         // 3) Effacer l’écran avant de redessiner (optionnel mais souvent nécessaire)
         glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
     
         // Redessiner et capturer l'image
         glfemPlotField(theGeometry->theElements, normDisplacement);
         glfemPlotMesh(theGeometry->theElements);
         glFlush();
         glFinish();
         captureFrame(i, window);
     
         // Swap des buffers et poll events
         glfwSwapBuffers(window);
         glfwPollEvents();
     
        }
     // 3) À la fin, libérer la mémoire
     free(X0);
     free(Y0);
 
     double hMin = femMin(normDisplacement,theNodes->nNodes);  
     double hMax = femMax(normDisplacement,theNodes->nNodes);  
     printf(" ==== Minimum displacement          : %14.7e [m] \n",hMin);
     printf(" ==== Maximum displacement          : %14.7e [m] \n",hMax);
 
     //
     //  -5- Calcul de la force globaleresultante
     //
 
     double theGlobalForce[2] = {0, 0};
     for (int i=0; i<theProblem->geometry->theNodes->nNodes; i++) {
         theGlobalForce[0] += theForces[2*i+0];
         theGlobalForce[1] += theForces[2*i+1]; }
     printf(" ==== Global horizontal force       : %14.7e [N] \n",theGlobalForce[0]);
     printf(" ==== Global vertical force         : %14.7e [N] \n",theGlobalForce[1]);
     printf(" ==== Weight                        : %14.7e [N] \n", area * rho * g);
 
     //
     //  -6- Visualisation du maillage
     //  
 
     
 
 
 
     // Check if the ESC key was pressed or the window was closed
 
     free(normDisplacement);
     free(forcesX);
     free(forcesY);
     femElasticityFree(theProblem) ; 
     geoFinalize();
     glfwTerminate(); 
     
 
 
     
     exit(EXIT_SUCCESS);
     return 0;  
 }
 
  
 
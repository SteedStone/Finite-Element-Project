/*
 *  main_visualization.c
 *  Visualisation avec itération sur la force appliquée et affichage de la rupture
 *  Utilisation de l'API de GMSH pour créer un maillage et de la fonction glfemPlotFailureNodes pour visualiser les ruptures
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 */

 #include <stdio.h>
 #include <stdlib.h>
 #include <string.h>
 #include <math.h>
 #include <time.h>
 #include "glfem.h"
 #define STB_IMAGE_WRITE_IMPLEMENTATION
 #include "stb_image_write.h"
 
 // Fonction de base pour une intégration (ici constante)
 double fun(double x, double y)
 {
     return 1;
 }
 
 void captureFrame(int frame, GLFWwindow* window) {
     int width, height;
     glfwGetFramebufferSize(window, &width, &height);
 
     printf("Framebuffer size: %d x %d\n", width, height);
     GLint viewport[4];
     glGetIntegerv(GL_VIEWPORT, viewport);
     printf("Viewport: x=%d, y=%d, width=%d, height=%d\n", viewport[0], viewport[1], viewport[2], viewport[3]);
     int winWidth, winHeight;
     glfwGetWindowSize(window, &winWidth, &winHeight);
     printf("Window size: %d x %d\n", winWidth, winHeight);
 
     unsigned char* pixels = (unsigned char*)malloc(3 * width * height);
     glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels);
 
     unsigned char* flippedPixels = (unsigned char*)malloc(3 * width * height);
     for (int y = 0; y < height; y++) {
         for (int x = 0; x < width; x++) {
             int srcIndex = 3 * (y * width + x);
             int destIndex = 3 * ((height - y - 1) * width + x);
             flippedPixels[destIndex]     = pixels[srcIndex];
             flippedPixels[destIndex + 1] = pixels[srcIndex + 1];
             flippedPixels[destIndex + 2] = pixels[srcIndex + 2];
         }
     }
 
     char filename[100];
     sprintf(filename, "frame/frame_%d.png", frame);
     stbi_write_png(filename, width, height, 3, flippedPixels, width * 3);
 
     free(pixels);
     free(flippedPixels);
 }
 
 int main(int argc, char *argv[])
 {
     if (argc < 2) {
         fprintf(stderr, "Usage: %s <path_to_mesh_file>\n", argv[0]);
         return EXIT_FAILURE;
     }
     char *file_path = argv[1];
 
     geoInitialize();
     femGeo* theGeometry = geoGetGeometry();
     geoMeshRead(file_path);
 
     FILE *fp = fopen("data/problem.txt", "r");
     if(fp == NULL) {
         perror("Erreur lors de l'ouverture du fichier problem.txt");
         return EXIT_FAILURE;
     }
     char ligne[256];
     double E = 211.e9, nu = 0.3, rho = 7.85e3, g = 9.81, deformation_factor = 1e3;
     double SigmaMax = 230e6; 
     char problemType[50] = {0};
 
     while(fgets(ligne, sizeof(ligne), fp) != NULL) {
         if(strstr(ligne, "Type of problem") != NULL) {
             sscanf(ligne, "Type of problem    : %49[^\n]", problemType);
         }
         else if(strstr(ligne, "Young modulus") != NULL) {
             sscanf(ligne, "Young modulus      : %le", &E);
         }
         else if(strstr(ligne, "Poisson ratio") != NULL) {
             sscanf(ligne, "Poisson ratio      : %le", &nu);
         }
         else if(strstr(ligne, "Mass density") != NULL) {
             sscanf(ligne, "Mass density       : %le", &rho);
         }
         else if(strstr(ligne, "Gravity") != NULL) {
             sscanf(ligne, "Gravity            : %le", &g);
         }
         else if(strstr(ligne, "Deformation Factor") != NULL) {
             sscanf(ligne, "Deformation Factor : %le", &deformation_factor);
         }
         else if(strstr(ligne, "Sigma Max (rupture)") != NULL) {
             sscanf(ligne, "Sigma Max (rupture)          : %le", &SigmaMax);
         }
     }
     fclose(fp);
 
     printf("Paramètres lus depuis problem.txt :\n");
     printf("  Type of problem    : %s\n", problemType);
     printf("  Young modulus      : %le\n", E);
     printf("  Poisson ratio      : %le\n", nu);
     printf("  Mass density       : %le\n", rho);
     printf("  Gravity            : %le\n", g);
     printf("  Deformation Factor : %le\n", deformation_factor);
     printf("  Sigma max (rupture): %le\n", SigmaMax);
 
     femNodes *theNodes = theGeometry->theNodes;
     double *X0 = malloc(theNodes->nNodes * sizeof(double));
     double *Y0 = malloc(theNodes->nNodes * sizeof(double));
     for (int i = 0; i < theNodes->nNodes; i++){
         X0[i] = theNodes->X[i];
         Y0[i] = theNodes->Y[i];
     }
 
     double *normDisplacement = malloc(theNodes->nNodes * sizeof(double));
     double *forcesX = malloc(theNodes->nNodes * sizeof(double));
     double *forcesY = malloc(theNodes->nNodes * sizeof(double));
     double *sigmaXX = malloc(theNodes->nNodes * sizeof(double));
     double *sigmaYY = malloc(theNodes->nNodes * sizeof(double));
     double *sigmaXY = malloc(theNodes->nNodes * sizeof(double));
     if(sigmaXX == NULL || sigmaYY == NULL || sigmaXY == NULL) {
         fprintf(stderr, "Erreur d'allocation des tableaux de contraintes\n");
         exit(EXIT_FAILURE);
     }
 
     glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);
     GLFWwindow* window = glfemInit("Visualisation et rupture (FEM)");
     glfwSetWindowSize(window, 1024, 768);
     glfwMakeContextCurrent(window);
 
     double F_top_start = (0.1e4)*4*4, F_top_end = (5e4)*4*4, F_top_step = (0.1e4)*4*4;
     int frameCounter = 0;
     for (double F_top = F_top_start; F_top <= F_top_end; F_top += F_top_step) {
         printf("Simulation avec Force Top = %.1e N\n", F_top);
         
         femProblem* theProblem = femElasticityCreate(theGeometry, E, nu, rho, g, PLANAR_STRESS, FEM_BAND, FEM_YNUM);
         femElasticityAddBoundaryCondition(theProblem, "Bottom", DIRICHLET_X, 0.0);
         femElasticityAddBoundaryCondition(theProblem, "Bottom", DIRICHLET_Y, 0.0);
         femElasticityAddBoundaryCondition(theProblem, "Top", NEUMANN_Y, -F_top);
 
         double *theSoluce = femElasticitySolve(theProblem);
         double *theForces = femElasticityForces(theProblem);
 
         femMesh *theMesh = theProblem->geometry->theElements;
        int *number = theMesh->nodes->number;
        for (int i=0; i<theNodes->nNodes; i++){
        theNodes->X[i] += theSoluce[2*number[i]+0]*deformation_factor;
        theNodes->Y[i] += theSoluce[2*number[i]+1]*deformation_factor;
        normDisplacement[i] = sqrt(theSoluce[2*number[i]+0]*theSoluce[2*number[i]+0] + 
                                theSoluce[2*number[i]+1]*theSoluce[2*number[i]+1]);
        forcesX[i] = theForces[2 * number[i] + 0];
        forcesY[i] = theForces[2 * number[i] + 1];
         }
 
         femElasticitySigma(theProblem, sigmaXX, sigmaYY, sigmaXY);
 
         int w, h;
         glfwGetFramebufferSize(window, &w, &h);
         glfemReshapeWindows(theNodes, w, h);
 
         glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
 
         glfemPlotField(theGeometry->theElements, normDisplacement);
         glfemPlotMesh(theGeometry->theElements);
 
         glfemPlotFailureNodes(theNodes, sigmaXX, sigmaYY, sigmaXY, SigmaMax);
 
         char theMessage[256];
         sprintf(theMessage, "Force Top: %.1e N", F_top);
         glColor3f(1.0, 0.0, 0.0);
         glfemDrawMessage(20, 375, theMessage);

 
         glfwSwapBuffers(window);
         glfwPollEvents();
         captureFrame(frameCounter, window);
         frameCounter++;
 
         femElasticityFree(theProblem);
     }
 
     free(normDisplacement);
     free(forcesX);
     free(forcesY);
     free(X0);
     free(Y0);
     free(sigmaXX);
     free(sigmaYY);
     free(sigmaXY);
     geoFinalize();
     glfwTerminate();
 
     return EXIT_SUCCESS;
 }
 
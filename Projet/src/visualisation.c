/*
 *  main.c
 *  Library for EPL1110 : Finite Elements for dummies
 *  Utilisation de l'API de GMSH pour créer un maillage
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 */

 #include "glfem.h"
 #define STB_IMAGE_WRITE_IMPLEMENTATION
 #include "stb_image_write.h"
 #define TRIANGLE 1
 #define HEXAGON 2
 
 double fun(double x, double y)
 {
     return 1;
 }
 
 // Capture le contenu du framebuffer et sauvegarde une image PNG
 void captureFrame(int frame, GLFWwindow* window) {
     int width, height;
     glfwGetFramebufferSize(window, &width, &height);
 
     // Debug
     printf("Framebuffer size: %d x %d\n", width, height);
     GLint viewport[4];
     glGetIntegerv(GL_VIEWPORT, viewport);
     printf("Viewport: x=%d, y=%d, width=%d, height=%d\n", viewport[0], viewport[1], viewport[2], viewport[3]);
     int winWidth, winHeight;
     glfwGetWindowSize(window, &winWidth, &winHeight);
     printf("Window size: %d x %d\n", winWidth, winHeight);
 
     // Allouer la mémoire pour récupérer les pixels
     unsigned char* pixels = (unsigned char*)malloc(3 * width * height);
     glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels);
 
     // Inverser l'image verticalement
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
 
     // Sauvegarder l'image dans le dossier "frame" (assure-toi que ce dossier existe)
     char filename[100];
     sprintf(filename, "frame/frame_%d.png", frame);
     stbi_write_png(filename, width, height, 3, flippedPixels, width * 3);
 
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
 
     // Initialisation de la géométrie et lecture du maillage
     geoInitialize();
     femGeo* theGeometry = geoGetGeometry();
     geoMeshRead(file_path);
 
     // Propriétés du matériau et gravité
     double E   = 211.e9;
     double nu  = 0.3;
     double rho = 7.85e3;
     double g   = 9.81;
     double deformation_factor = 1e3;  // Ajuste cette valeur selon la visualisation désirée

     // Stocker les coordonnées initiales des nœuds
     femNodes *theNodes = theGeometry->theNodes;
     double *X0 = malloc(sizeof(double) * theNodes->nNodes);
     double *Y0 = malloc(sizeof(double) * theNodes->nNodes);
     for (int j = 0; j < theNodes->nNodes; j++) {
         X0[j] = theNodes->X[j];
         Y0[j] = theNodes->Y[j];
     }
 
     // Allocation pour les champs de déformation et forces (affichage)
     double *normDisplacement = malloc(theNodes->nNodes * sizeof(double));
     double *forcesX = malloc(theNodes->nNodes * sizeof(double));
     double *forcesY = malloc(theNodes->nNodes * sizeof(double));
 
     // Création de la fenêtre GLFW avec une résolution plus grande
     glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);
     GLFWwindow* window = glfemInit("EPL1110 : Recovering forces on constrained nodes");
     glfwSetWindowSize(window, 1024, 768);
     glfwMakeContextCurrent(window);
 
     // Variation de la force sur le bord supérieur de 1e3 N à 5e3 N
     double F_top_start = 1e5;
     double F_top_end   = 1e6;
     double F_top_step  = 2e4; // incréments de 200 N
     int frameCounter = 0;
     
    double *stressValues = malloc(theNodes->nNodes * sizeof(double));
    for (int j = 0; j < theNodes->nNodes; j++) {
        stressValues[j] = 0.0; // Initialisation
    }

        // Boucle d'animation pour chaque valeur de force
    for (double F_top = F_top_start; F_top <= F_top_end; F_top += F_top_step) {
        printf("Simulation avec force Top = %.1e N\n", F_top);

        // Créer le problème FEM pour cette valeur de force
        femProblem* theProblem = femElasticityCreate(theGeometry, E, nu, rho, g, PLANAR_STRESS, FEM_FULL, FEM_NO);
        // Conditions aux limites :
        femElasticityAddBoundaryCondition(theProblem, "Bottom", DIRICHLET_X, 0.0 ,0.0);
        femElasticityAddBoundaryCondition(theProblem, "Bottom", DIRICHLET_Y, 0.0,0.0);
        // Force appliquée sur le bord supérieur (négative pour une force vers le bas)
        femElasticityAddBoundaryCondition(theProblem, "Top", NEUMANN_Y, -F_top , -F_top);

        femElasticityPrint(theProblem);

        // Résoudre le problème pour obtenir la déformation
        double *theSoluce = femElasticitySolve(theProblem);
        double *theForces = femElasticityForces(theProblem);

        // Mettre à jour les positions des nœuds en repartant des positions initiales
        for (int j = 0; j < theNodes->nNodes; j++) {
            theNodes->X[j] = X0[j] + deformation_factor * theSoluce[2*j + 0];
            theNodes->Y[j] = Y0[j] + deformation_factor * theSoluce[2*j + 1];
            normDisplacement[j] = sqrt(theSoluce[2*j + 0] * theSoluce[2*j + 0] +
                                        theSoluce[2*j + 1] * theSoluce[2*j + 1]);
            forcesX[j] = theForces[2*j + 0];
            forcesY[j] = theForces[2*j + 1];
        }

        // -------- Calcul de la déformation et contrainte pour toutes les lignes --------
        // On parcourt chaque élément, et pour chaque côté (ligne) de l'élément, on calcule
        // L0 (longueur initiale) et L (longueur déformée)
        int nElem = theGeometry->theElements->nElem;
        int nLocal = theGeometry->theElements->nLocalNode;
        // Pour chaque élément, pour chaque côté (de k à (k+1)%nLocal)
        for (int iElem = 0; iElem < nElem; iElem++) {
            int *elem = &(theGeometry->theElements->elem[iElem * nLocal]);
            for (int k = 0; k < nLocal; k++) {
                int node1 = elem[k];
                int node2 = elem[(k + 1) % nLocal];
                double L0 = sqrt((X0[node2] - X0[node1]) * (X0[node2] - X0[node1]) +
                                 (Y0[node2] - Y0[node1]) * (Y0[node2] - Y0[node1]));
                double L = sqrt((theNodes->X[node2] - theNodes->X[node1]) * (theNodes->X[node2] - theNodes->X[node1]) +
                                (theNodes->Y[node2] - theNodes->Y[node1]) * (theNodes->Y[node2] - theNodes->Y[node1]));
                double strain = (L - L0) / L0;
                double stress = E * strain;  // loi de Hooke
                // Ajouter la contrainte aux nœuds (moyenne des contributions des éléments voisins)
                stressValues[node1] += stress;
                stressValues[node2] += stress;


                // Afficher l'information pour cette ligne
                printf("Element %d, ligne entre noeuds %d et %d : L0 = %e, L = %e, strain = %e, stress = %e Pa\n",
                       iElem, node1, node2, L0, L, strain, stress);
                if (stress > 300e20) {
                    printf("  ⚠️  Risque de rupture sur cette ligne !\n");
                }
            }
        }
        // -------- Fin calcul contraintes --------
        for (int j = 0; j < theNodes->nNodes; j++) {
            stressValues[j] /= theNodes->nNodes;  // Normalisation par le nombre de contributions
        }
        
        int w, h;
        glfwGetFramebufferSize(window, &w, &h);
        glfemReshapeWindows(theGeometry->theNodes, w, h);

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glfemPlotField(theGeometry->theElements, stressValues);

        glfemPlotMesh(theGeometry->theElements);
        glFlush();
        glFinish();

        captureFrame(frameCounter, window);
        frameCounter++;

        glfwSwapBuffers(window);
        glfwPollEvents();

        // Ne pas appeler free sur theSoluce et theForces si c'est géré par femElasticityFree
        femElasticityFree(theProblem);
    }

    // Affichage final (optionnel)
    double hMin = femMin(normDisplacement, theNodes->nNodes);
    double hMax = femMax(normDisplacement, theNodes->nNodes);
    printf(" ==== Minimum displacement          : %14.7e [m] \n", hMin);
    printf(" ==== Maximum displacement          : %14.7e [m] \n", hMax);

    double theGlobalForce[2] = {0, 0};
    for (int i = 0; i < theGeometry->theNodes->nNodes; i++) {
        theGlobalForce[0] += forcesX[i];
        theGlobalForce[1] += forcesY[i];
    }
    printf(" ==== Global horizontal force       : %14.7e [N] \n", theGlobalForce[0]);
    printf(" ==== Global vertical force         : %14.7e [N] \n", theGlobalForce[1]);

    free(normDisplacement);
    free(forcesX);
    free(forcesY);
    free(X0);
    free(Y0);
    geoFinalize();
    glfwTerminate();

    exit(EXIT_SUCCESS);
    return 0;
}
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
#include <time.h>
double fun(double x, double y) 
{
    return 1;
}




int main(int argc, char *argv[])
{  
    printf("\n\n    V : Mesh and size mesh field \n");
    printf("    D : Domains \n");
    printf("    N : Next domain highlighted\n");

    

    if (argc < 3) {
        fprintf(stderr, "Usage: %s <path_to_mesh> <path_to_problem>\n", argv[0]);
        return EXIT_FAILURE;
    }

    char *file_mesh = argv[1];
    char *file_problem = argv[2];

    


    // On va devoir avoir un gros rectangle de lx = 2 et ly = 0.5
    // Puis deux petits cercles au bout de chaque côté du rectangle de rayon 0.1 et 0.2 
    // et enfin deux rectangle dans les cercles. 
 
   
      
    int ierr;
    femGeo* theGeometry ;
    geoInitialize();
    theGeometry = geoGetGeometry();
    
    
    

    geoMeshRead(file_mesh);
    

    
    FILE *fp = fopen(file_problem, "r");
    if(fp == NULL) {
        perror("Erreur lors de l'ouverture du fichier problem.txt");
        return EXIT_FAILURE;
    }
    
    char ligne[256];
    double E = 0.0, nu = 0.0, rho = 0.0, g = 0.0, deformationFactor = 0.0, SigmaMax = 0.0;
    int problemCase = PLANAR_STRESS; // Valeur par défaut
    char problemType[50] = {0};      // Pour stocker la chaîne lue

    // Lecture du fichier ligne par ligne
    while(fgets(ligne, sizeof(ligne), fp) != NULL) {
        if(strstr(ligne, "Type of problem") != NULL) {
            // Exemple de ligne : "Type of problem    :  Planar stresses"
            // La directive %[^\n] permet de récupérer toute la chaîne jusqu'au retour à la ligne.
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
            sscanf(ligne, "Deformation Factor : %le", &deformationFactor);
        }
        else if(strstr(ligne, "Sigma Max (rupture)") != NULL) {
            sscanf(ligne, "Sigma Max (rupture): %le", &SigmaMax);
        }
    }
    fclose(fp);
    
    // Traitement du type de problème pour définir la variable problemCase
    if(strstr(problemType, "stresses") != NULL) {
        problemCase = PLANAR_STRESS;
    }
    else if(strstr(problemType, "strains") != NULL) {
        problemCase = PLANAR_STRAIN;
    }

    else {
        // Si le type n'est pas reconnu, on garde la valeur par défaut
        problemCase = PLANAR_STRESS;
    }
    
    // Affichage des valeurs lues
    printf("Type of problem    : %s\n", problemType);
    printf("Young modulus      : %le\n", E);
    printf("Poisson ratio      : %le\n", nu);
    printf("Mass density       : %le\n", rho);
    printf("Gravity            : %le\n", g);
    printf("Deformation Factor : %le\n", deformationFactor);
    printf("Sigma Max (rupture): %le\n", SigmaMax);
    
    // Exemple d'utilisation avec la fonction femElasticityCreate :
    // femProblem* theProblem = femElasticityCreate(theGeometry, E, nu, rho, g, problemCase, FEM_BAND, FEM_NO);
    



    
    // Initialisation du problème avec les conditions aux bords
    // ON remplit juste la structure theProblem avec les valeurs de E, nu, rho, g et le type de problème qu'on veut résoudre
    femProblem* theProblem = femElasticityCreate(theGeometry,E,nu,rho,g,PLANAR_STRESS , FEM_BAND , FEM_YNUM);
    // femElasticityAddBoundaryCondition(theProblem,"Symmetry",DIRICHLET_X,0.0);
    // femElasticityAddBoundaryCondition(theProblem,"Bottom",DIRICHLET_Y,0.0);
    // femElasticityAddBoundaryCondition(theProblem,"Bottom",DIRICHLET_X,0.0);


    // femElasticityAddBoundaryCondition(theProblem,"Top",NEUMANN_Y,-1e6);
    fp = fopen(file_problem, "r");
    while(fgets(ligne, sizeof(ligne), fp) != NULL) {
        if(strstr(ligne, "Boundary condition") != NULL) {
            char typeBC[20];
            double valeur;
            char nomDomaine[20];
            // Exemple de format : "Boundary condition :  Dirichlet-X      =  0.0000000e+00 : Bottom"
            sscanf(ligne, "Boundary condition : %s = %le : %s", typeBC, &valeur, nomDomaine);
            // Selon le type, vous appellerez la fonction correspondante
            if(strcmp(typeBC, "Dirichlet-X") == 0)
                femElasticityAddBoundaryCondition(theProblem, nomDomaine, DIRICHLET_X, valeur);
            else if(strcmp(typeBC, "Dirichlet-Y") == 0)
                femElasticityAddBoundaryCondition(theProblem, nomDomaine, DIRICHLET_Y, valeur);
            else if(strcmp(typeBC, "Neumann-Y") == 0)
                femElasticityAddBoundaryCondition(theProblem, nomDomaine, NEUMANN_Y, valeur);
            // Ajoutez les autres cas si nécessaire.
        }
    }
    fclose(fp);

    
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
    int normal = 0;


    femNodes *theNodes = theGeometry->theNodes;
    double *normDisplacement = malloc(theNodes->nNodes * sizeof(double));
    double *forcesX = malloc(theNodes->nNodes * sizeof(double));
    double *forcesY = malloc(theNodes->nNodes * sizeof(double));


   

    
    femMesh *theMesh = theProblem->geometry->theElements;
    int *number = theMesh->nodes->number;
    for (int i=0; i<theNodes->nNodes; i++){
        theNodes->X[i] += theSoluce[2*number[i]+0]*deformationFactor;
        theNodes->Y[i] += theSoluce[2*number[i]+1]*deformationFactor;
        normDisplacement[i] = sqrt(theSoluce[2*number[i]+0]*theSoluce[2*number[i]+0] + 
                                theSoluce[2*number[i]+1]*theSoluce[2*number[i]+1]);
        forcesX[i] = theForces[2 * number[i] + 0];
        forcesY[i] = theForces[2 * number[i] + 1];
         }
    int nNodes = theNodes->nNodes;

    
    
    // N'oublie pas de libérer le tableau inverse

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
    
    double *sigmaXX = (double *) malloc(nNodes * sizeof(double));
    double *sigmaYY = (double *) malloc(nNodes * sizeof(double));
    double *sigmaXY = (double *) malloc(nNodes * sizeof(double));

    if (sigmaXX == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return EXIT_FAILURE; }
    if (sigmaYY == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return EXIT_FAILURE; }
    if (sigmaXY == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return EXIT_FAILURE; }

    femElasticitySigma(theProblem, sigmaXX, sigmaYY, sigmaXY);
    
    

    femSolverType solverType = FEM_FULL;
    femRenumType  renumType  = FEM_NO;
    int option = 1;    
    int mode = 8;
    femSolverType newSolverType = solverType;
    femRenumType  newRenumType  = renumType;

    GLFWwindow* window = glfemInit("LEPL1110 : Band Solver ");
    glfwMakeContextCurrent(window);
    glPointSize(4.0f);                 // Par exemple, un point de 4 pixels
    glEnable(GL_POINT_SMOOTH);         // Active l’anticrénelage pour les points
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST); 
    
    do 
    {
        int testConvergence;
        char theMessage[256];
        // sprintf(theMessage, "Max : %.4f ",femMax(theProblem->soluce,theProblem->size));
        int w,h;
        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindows(theGeometry->theNodes,w,h);
    
        if (option == 1) {
            glfemPlotField(theGeometry->theElements,normDisplacement);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); 
        }
        else {
            glColor3f(1.0,0.0,0.0);
            glfemPlotSolver(theProblem->solver,theProblem->size,w,h); }
        glColor3f(1.0,0.0,0.0); glfemDrawMessage(20,460,theMessage);              
    
        if (solverType != newSolverType || renumType != newRenumType && option == 0) { 
            solverType = newSolverType;
            renumType = newRenumType;
            femElasticityFree(theProblem);
            theProblem = femElasticityCreate(theGeometry,E,nu,rho,g,PLANAR_STRESS , solverType , renumType);
            femElasticityAddBoundaryCondition(theProblem,"Bottom",DIRICHLET_X,0.0);
            femElasticityAddBoundaryCondition(theProblem,"Bottom",DIRICHLET_Y,0.0);
            femElasticityAddBoundaryCondition(theProblem,"Top",NEUMANN_Y,-1e6);
            theMesh = theProblem->geometry->theElements;
            clock_t tic = clock();
            do {
                femElasticitySolve(theProblem);  
                femSolverPrintInfos(theProblem->solver); 
                testConvergence = femSolverConverged(theProblem->solver); }
            while ( testConvergence == 0);
            if (testConvergence == -1)  printf("    Iterative solver stopped afer a maximum number of iterations\n");
            printf("    CPU time : %.2f [sec] \n", (clock() - tic) * 1.0 /CLOCKS_PER_SEC);
            switch (renumType) {
                case FEM_XNUM : printf("    Renumbering along the x-direction\n"); break;
                case FEM_YNUM : printf("    Renumbering along the y-direction\n"); break;
                default : break; }
            // printf("    Maximum value : %.4f\n", femMax(theProblem->soluce,theProblem->size));
            fflush(stdout); }
        if (glfwGetKey(window,'V') == GLFW_PRESS)   {option = 1;mode = 8;}
        if (glfwGetKey(window,'S') == GLFW_PRESS)   option = 0;
        if (glfwGetKey(window,'F') == GLFW_PRESS)   newSolverType = FEM_FULL; 
        if (glfwGetKey(window,'B') == GLFW_PRESS)   newSolverType = FEM_BAND; 
        // if (glfwGetKey(window,'I') == GLFW_PRESS)   newSolverType = FEM_ITER; 
        if (glfwGetKey(window,'X') == GLFW_PRESS)   newRenumType  = FEM_XNUM; 
        if (glfwGetKey(window,'Y') == GLFW_PRESS)   newRenumType  = FEM_YNUM; 
        if (glfwGetKey(window,'N') == GLFW_PRESS)   newRenumType  = FEM_NO; 
        if (glfwGetKey(window, 'X') == GLFW_PRESS) { mode = 3; }
        if (glfwGetKey(window, 'Y') == GLFW_PRESS) { mode = 4; }
        if (glfwGetKey(window, 'K') == GLFW_PRESS) { mode = 5; }
        if (glfwGetKey(window, 'L') == GLFW_PRESS) { mode = 6; }
        if (glfwGetKey(window, 'M') == GLFW_PRESS) { mode = 7; }
        if (glfwGetKey(window, 'T') == GLFW_PRESS) {glfemPlotFailureNodes(theGeometry->theNodes, sigmaXX, sigmaYY, sigmaXY, SigmaMax); }

        if(option == 1 && mode == 3) {
            glfemPlotField(theGeometry->theElements, forcesX);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
        }
        else if(option == 1 && mode == 4) {
            glfemPlotField(theGeometry->theElements, forcesY);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
        }
        else if(option == 1 && mode == 5) {
            glfemPlotField(theGeometry->theElements, sigmaYY);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
        }
        else if(option == 1 && mode == 6) {
            glfemPlotField(theGeometry->theElements, sigmaXX);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
        }
        else if(option == 1 && mode == 7) {
            glfemPlotField(theGeometry->theElements, sigmaXY);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
        
        } else {
            glfemPlotField(theGeometry->theElements,normDisplacement);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); 
        }
        glfwSwapBuffers(window);
        glfwPollEvents();
    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
            glfwWindowShouldClose(window) != 1 );
            
    // Check if the ESC key was pressed or the window was closed
            
    glfwTerminate(); 




    


    
    exit(EXIT_SUCCESS);
    return 0;  
}

 

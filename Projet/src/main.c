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
    // femElasticityAddBoundaryCondition(theProblem,"Symmetry",DIRICHLET_X,0.0);
    femElasticityAddBoundaryCondition(theProblem,"Bottom",DIRICHLET_X,0.0);
    femElasticityAddBoundaryCondition(theProblem,"Bottom",DIRICHLET_Y,0.0);

    femElasticityAddBoundaryCondition(theProblem,"Top",NEUMANN_Y,-1e6);
    
    femElasticityPrint(theProblem);
    
    // femFullSystemPrint(theProblem->solver->solver);
    // printf("Taille system %d" , theProblem->size);
    
    //  -3- Resolution du probleme et calcul des forces
    

    double *theSoluce = femElasticitySolve(theProblem);
    printf("fin de force") ;

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


   
    double deformationFactor = 300.0; // 5000.0 pour hexa et 300 pour triangle
    
    femMesh *theMesh = theProblem->geometry->theElements;
    int *number = theMesh->nodes->number;
    for (int i=0; i<theNodes->nNodes; i++){
        theNodes->X[i] += theSoluce[2*i+0]*deformationFactor;
        theNodes->Y[i] += theSoluce[2*i+1]*deformationFactor;
        normDisplacement[i] = sqrt(theSoluce[2*i+0]*theSoluce[2*i+0] + 
                                theSoluce[2*i+1]*theSoluce[2*i+1]);
        forcesX[i] = theForces[2*i+0];
        forcesY[i] = theForces[2*i+1]; }
    
    

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

    // int mode = 1; 
    // int domain = 0;
    // int freezingButton = FALSE;
    // double t, told = 0;
    // char theMessage[MAXNAME];


    // GLFWwindow* window = glfemInit("EPL1110 : Recovering forces on constrained nodes");
    // glfwMakeContextCurrent(window);

    // do {
    //     int w,h;
    //     glfwGetFramebufferSize(window,&w,&h);
    //     glfemReshapeWindows(theGeometry->theNodes,w,h);

    //     t = glfwGetTime();  
    //     if (glfwGetKey(window,'D') == GLFW_PRESS) { mode = 0;}
    //     if (glfwGetKey(window,'V') == GLFW_PRESS) { mode = 1;}
    //     if (glfwGetKey(window,'X') == GLFW_PRESS) { mode = 2;}
    //     if (glfwGetKey(window,'Y') == GLFW_PRESS) { mode = 3;}
    //     if (glfwGetKey(window,'N') == GLFW_PRESS && freezingButton == FALSE) { domain++; freezingButton = TRUE; told = t;}
    //     if (t-told > 0.5) {freezingButton = FALSE; }
        
    //     if (mode == 0) {
    //         domain = domain % theGeometry->nDomains;
    //         glfemPlotDomain( theGeometry->theDomains[domain]); 
    //         sprintf(theMessage, "%s : %d ",theGeometry->theDomains[domain]->name,domain);
    //         glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
    //     if (mode == 1) {
    //         glfemPlotField(theGeometry->theElements,normDisplacement);
    //         glfemPlotMesh(theGeometry->theElements); 
    //         sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
    //         glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
    //     if (mode == 2) {
    //         glfemPlotField(theGeometry->theElements,forcesX);
    //         glfemPlotMesh(theGeometry->theElements); 
    //         sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
    //         glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
    //     if (mode == 3) {
    //         glfemPlotField(theGeometry->theElements,forcesY);
    //         glfemPlotMesh(theGeometry->theElements); 
    //         sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
    //         glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
    //     glfwSwapBuffers(window);
    //     glfwPollEvents();
    // } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
    //         glfwWindowShouldClose(window) != 1 );
            
    // // Check if the ESC key was pressed or the window was closed

    // free(normDisplacement);
    // free(forcesX);
    // free(forcesY);
    // femElasticityFree(theProblem) ; 
    // geoFinalize();
    // glfwTerminate(); 

    femSolverType solverType = FEM_FULL;
    femRenumType  renumType  = FEM_NO;
    int option = 1;    
    int mode = 5;
    femSolverType newSolverType = solverType;
    femRenumType  newRenumType  = renumType;

    GLFWwindow* window = glfemInit("LEPL1110 : Band Solver ");
    glfwMakeContextCurrent(window);
    
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
        if (glfwGetKey(window,'V') == GLFW_PRESS)   {option = 1;mode = 5;}
        if (glfwGetKey(window,'S') == GLFW_PRESS)   option = 0;
        if (glfwGetKey(window,'F') == GLFW_PRESS)   newSolverType = FEM_FULL; 
        if (glfwGetKey(window,'B') == GLFW_PRESS)   newSolverType = FEM_BAND; 
        // if (glfwGetKey(window,'I') == GLFW_PRESS)   newSolverType = FEM_ITER; 
        if (glfwGetKey(window,'X') == GLFW_PRESS)   newRenumType  = FEM_XNUM; 
        if (glfwGetKey(window,'Y') == GLFW_PRESS)   newRenumType  = FEM_YNUM; 
        if (glfwGetKey(window,'N') == GLFW_PRESS)   newRenumType  = FEM_NO; 
        if (glfwGetKey(window, 'X') == GLFW_PRESS) { mode = 3; }
        if (glfwGetKey(window, 'Y') == GLFW_PRESS) { mode = 4; }
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

//     int mode = 1;
//     int domain = 0;
//     int freezingButton = FALSE;
//     double t, told = 0;
//     char theMessage[MAXNAME + 5];

//     GLFWwindow *window = glfemInit("EPL1110 : Project 2022-23 ");
//     glfwMakeContextCurrent(window);
//     // glfwSetScrollCallback(window, scroll_callback);
//     // glfwSetMouseButtonCallback(window, mouse_button_callback);

//     do {
//         int w, h;
//         glfwGetFramebufferSize(window, &w, &h);
//         glfemReshapeWindows(theGeometry->theNodes, w, h);

//         t = glfwGetTime();
//         if (glfwGetKey(window, 'D') == GLFW_PRESS) { mode = 0; }
//         if (glfwGetKey(window, 'V') == GLFW_PRESS) { mode = 1; }
//         if (glfwGetKey(window, 'S') == GLFW_PRESS) { mode = 2; }
//         if (glfwGetKey(window, 'X') == GLFW_PRESS) { mode = 3; }
//         if (glfwGetKey(window, 'Y') == GLFW_PRESS) { mode = 4; }
//         if (glfwGetKey(window, 'U') == GLFW_PRESS) { mode = 5; }
//         if (glfwGetKey(window, 'I') == GLFW_PRESS) { mode = 6; }
//         if (glfwGetKey(window, 'O') == GLFW_PRESS) { mode = 7; }
//         if (glfwGetKey(window, 'P') == GLFW_PRESS) { mode = 8; }
//         if (glfwGetKey(window, 'N') == GLFW_PRESS && freezingButton == FALSE)
//         {
//             domain++;
//             freezingButton = TRUE;
//             told = t;
//         }

//         if (t - told > 0.5) { freezingButton = FALSE; }
//         if (mode == 1)
//         {
//             glfemPlotField(theGeometry->theElements, normDisplacement);
//             glfemPlotMesh(theGeometry->theElements);
//             sprintf(theMessage, "Number of elements : %d ", theGeometry->theElements->nElem);
//             glColor3f(1.0, 0.0, 0.0);
//             glfemMessage(theMessage);
//         }
//         else if (mode == 0)
//         {
//             domain = domain % theGeometry->nDomains;
//             glfemPlotDomain(theGeometry->theDomains[domain]);
//             sprintf(theMessage, "%s : %d ", theGeometry->theDomains[domain]->name, domain);
//             glColor3f(1.0, 0.0, 0.0);
//             glfemMessage(theMessage);
//         } 
//         else if (mode == 2)
//         {
//             glColor3f(1.0, 0.0, 0.0);
//             glfemPlotSolver(theProblem->solver, theProblem->size, w, h);
//         }
//         else if (mode == 3)
//         {
//             glfemPlotField(theGeometry->theElements, forcesX);
//             glfemPlotMesh(theGeometry->theElements); 
//             sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
//             glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
//         }
//         else if (mode == 4)
//         {
//             glfemPlotField(theGeometry->theElements, forcesY);
//             glfemPlotMesh(theGeometry->theElements); 
//             sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
//             glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
//         }
//         // else if (mode == 5)
//         // {
//         //     glfemPlotField(theGeometry->theElements, sigmaXX);
//         //     glfemPlotMesh(theGeometry->theElements); 
//         //     sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
//         //     glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
//         // }
//         // else if (mode == 6)
//         // {
//         //     glfemPlotField(theGeometry->theElements, sigmaYY);
//         //     glfemPlotMesh(theGeometry->theElements); 
//         //     sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
//         //     glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
//         // }
//         // else if (mode == 7)
//         // {
//         //     glfemPlotField(theGeometry->theElements, sigmaXY);
//         //     glfemPlotMesh(theGeometry->theElements); 
//         //     sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
//         //     glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
//         // }
//         // else if (mode == 8)
//         // {
//         //     glfemPlotField(theGeometry->theElements, vonMises);
//         //     glfemPlotMesh(theGeometry->theElements); 
//         //     sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
//         //     glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
//         // }
//         else { printf("Mode %d not implemented\n", mode); }  

//         glfwSwapBuffers(window);
//         glfwPollEvents();
//     } while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS && glfwWindowShouldClose(window) != 1);
//     // Check if the ESC key was pressed or the window was closed

//    free(normDisplacement);
//     free(forcesX);
//     free(forcesY);
//     femElasticityFree(theProblem) ; 
//     geoFinalize();
//     glfwTerminate(); 


    


    
    exit(EXIT_SUCCESS);
    return 0;  
}

 

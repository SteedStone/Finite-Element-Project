/*
 *  fem.h
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

 #ifndef _FEM_H_
 #define _FEM_H_
 
 #include <stdio.h>
 #include <stdlib.h>
 #include <math.h>
 #include <string.h>
 #include <limits.h>
 #include "gmshc.h"
 
 #ifndef M_PI
 #define M_PI 3.14159265358979323846
 #endif
 
 // Macros d'erreur et de gestion
 #define ErrorScan(a)   femErrorScan(a,__LINE__,__FILE__)
 #define ErrorGmsh(a)   femErrorGmsh(a,__LINE__,__FILE__)
 #define Error(a)       femError(a,__LINE__,__FILE__)
 #define Warning(a)     femWarning(a,  __LINE__, __FILE__)
 #define FALSE 0 
 #define TRUE  1
 #define MAXNAME 256
 
 // ======================================================================
 // ->Enums
 // ======================================================================
 typedef enum { FEM_TRIANGLE, FEM_QUAD, FEM_EDGE } femElementType;
 typedef enum { DIRICHLET_X, DIRICHLET_Y, NEUMANN_X, NEUMANN_Y } femBoundaryType;
 typedef enum { PLANAR_STRESS, PLANAR_STRAIN, AXISYM } femElasticCase;
 typedef enum { FEM_FULL, FEM_BAND, FEM_ITER } femSolverType;
 typedef enum { FEM_NO, FEM_XNUM, FEM_YNUM } femRenumType;
 
 // ======================================================================
 // ->Structures de données
 // ======================================================================
 
 // 1. Mesh et Domaines
 typedef struct {
     int nNodes;
     double *X;
     double *Y;
     int *number;
 } femNodes;
 
 typedef struct {
     int nLocalNode;
     int nElem;
     int *elem;
     femNodes *nodes;
 } femMesh;
 
 typedef struct {
     femMesh *mesh;
     int nElem;
     int *elem;
     char name[MAXNAME];
 } femDomain;
 
 // 2. Géométrie
 typedef struct {
     int NumberOfHexagonsInX, NumberOfHexagonsInY;
     int NumberOfTrianglesInX, NumberOfTrianglesInY;
     double hexRadius;
     double h;
     double MiddleX;
     double MiddleY;
     femElementType elementType;
     double (*geoSize)(double x, double y);
     femNodes *theNodes;
     femMesh  *theElements;
     femMesh  *theEdges;
     int nDomains;
     femDomain **theDomains;
     double triBase;   // Base des triangles
     double triHeight;
     int hexa_triangles;
     double LxPlate;
     double LyPlate;
 } femGeo;
 
 // 3. Discrétisation et Intégration
 typedef struct {
     int n;
     femElementType type;
     void (*x2)(double *xsi, double *eta);
     void (*phi2)(double xsi, double eta, double *phi);
     void (*dphi2dx)(double xsi, double eta, double *dphidxsi, double *dphideta);
     void (*x)(double *xsi);
     void (*phi)(double xsi, double *phi);
     void (*dphidx)(double xsi, double *dphidxsi);
 } femDiscrete;
 
 typedef struct {
     int n;
     const double *xsi;
     const double *eta;
     const double *weight;
 } femIntegration;
 
 // 4. Systèmes linéaires
 typedef struct {
     double *B;
     double **A;
     int size;
 } femFullSystem;
 
 typedef struct {
     double *B;
     double **A;        
     int size;
     int band;        
 } femBandSystem;
 
 typedef struct {
     double *R;
     double *D;
     double *S;
     double *X; 
     double error;      
     int size;
     int iter;        
 } femIterativeSolver;
 
 typedef struct {
     femSolverType type;
     femFullSystem *local;
     void *solver;
     int size;
 } femSolver;
 
 // 5. Conditions aux limites
 typedef struct {
     femDomain* domain;
     femBoundaryType type; 
     double value;
 } femBoundaryCondition;
 
 // 6. Problème de calcul
 typedef struct {
     double E, nu, rho, g;
     double A, B, C;
     int planarStrainStress;
     int nBoundaryConditions;
     femBoundaryCondition **conditions;  
     int *constrainedNodes; 
     double *soluce;
     double *residuals;
     femGeo *geometry;
     femDiscrete *space;
     femIntegration *rule;
     femDiscrete *spaceEdge;
     femIntegration *ruleEdge;
     femSolver *solver;
     int size;
     int sizeLoc;
 } femProblem;
 
 // ======================================================================
 // ->Prototypes des fonctions
 // ======================================================================
 
 // -------------------------
 // 1. Gestion de la Géométrie
 // -------------------------
 void                geoInitialize();
 femGeo*             geoGetGeometry();
 
 double              geoSize(double x, double y);
 double              geoSizeDefault(double x, double y);
 
 void                geoSetSizeCallback(double (*geoSize)(double x, double y));
 
 void                geoMeshGenerate();
 void                geoMeshImport();
 void                geoMeshPrint();
 void                geoMeshWrite(const char *filename);
 void                geoMeshRead(const char *filename);
 
 void                geoSetDomainName(int iDomain, char *name);
 int                 geoGetDomain(char *name);
 
 void                geoFinalize();
 
 // -------------------------
 // 2. Problème d'Élasticité
 // -------------------------
 femProblem*         femElasticityCreate(femGeo* theGeometry, 
                                          double E, double nu, double rho, double g, 
                                          femElasticCase iCase, femSolverType solverType, femRenumType renumType);
 void                femElasticityFree(femProblem *theProblem);
 void                femElasticityPrint(femProblem *theProblem);
 void                femElasticityAddBoundaryCondition(femProblem *theProblem, char *nameDomain, femBoundaryType type, double value);
 void                femElasticityAssembleElements(femProblem *theProblem);
 void                femElasticityAssembleNeumann(femProblem *theProblem);
 double*             femElasticitySolve(femProblem *theProblem);
 double*             femElasticityForces(femProblem *theProblem);
 double              femElasticityIntegrate(femProblem *theProblem, double (*f)(double x, double y));
 void femSolverGetResidual(femSolver *mySolver, double *residuals, double *theSoluce) ;
 void femFullSystemGetResidual(femFullSystem *mySystem, int size, double *residuals, double *theSoluce) ;
 void femBandSystemGetResidual(femBandSystem *mySystem, int size, double *residuals, double *theSoluce) ;



 // -------------------------
 // 2. Discrétisation et Intégration
 // -------------------------
 femIntegration*     femIntegrationCreate(int n, femElementType type);
 void                femIntegrationFree(femIntegration *theRule);
 
 femDiscrete*        femDiscreteCreate(int n, femElementType type);
 void                femDiscreteFree(femDiscrete* mySpace);
 void                femDiscretePrint(femDiscrete* mySpace);
 
 void                femDiscreteXsi2(femDiscrete* mySpace, double *xsi, double *eta);
 void                femDiscretePhi2(femDiscrete* mySpace, double xsi, double eta, double *phi);
 void                femDiscreteDphi2(femDiscrete* mySpace, double xsi, double eta, double *dphidxsi, double *dphideta);
 
 void                femDiscreteXsi(femDiscrete* mySpace, double *xsi);
 void                femDiscretePhi(femDiscrete* mySpace, double xsi, double *phi);
 void                femDiscreteDphi(femDiscrete* mySpace, double xsi, double *dphidxsi);
 
 // -------------------------
 // 4. Systèmes Linéaires (Matrices et Assemblage)
 // -------------------------
 /* 4.1 Système complet (Full System) */
 femFullSystem*      femFullSystemCreate(int size);
 void                femFullSystemFree(femFullSystem* mySystem);
 void                femFullSystemPrint(femFullSystem* mySystem);
 void                femFullSystemInit(femFullSystem* mySystem);
 void                femFullSystemAlloc(femFullSystem* mySystem, int size);
 double*             femFullSystemEliminate(femFullSystem* mySystem);
 void                femFullSystemConstrain(femFullSystem* mySystem, int myNode, double value, int size);
 double              femFullSystemGet(femFullSystem* myFullSystem, int myRow, int myCol);
 void                femFullSystemAssemble(femFullSystem *mySystem, double *Aloc, double *Bloc, int *mapX, int *mapY, int nLocal);
 
 /* 4.2 Système bande (Band System) */
 femBandSystem*       femBandSystemCreate(int size, int band);
 void                 femBandSystemFree(femBandSystem* myBandSystem);
 void                 femBandSystemInit(femBandSystem *myBand);
 void                 femBandSystemPrint(femBandSystem *myBand);
 void                 femBandSystemPrintInfos(femBandSystem *myBand);
 double*              femBandSystemEliminate(femBandSystem *myBand);
 void                 femBandSystemAssemble(femBandSystem *myBandSystem, double *Aloc, double *Bloc, int *mapX, int *mapY, int nLoc);
 double               femBandSystemGet(femBandSystem* myBandSystem, int i, int j);
 int                  femMeshComputeBand(femMesh *theMesh);
 void femBandSystemConstrain(femBandSystem *mySystem, int myNode, double myValue, int size) ;
 double femBandSystemGetA_Entry(femBandSystem *mySystem, int myRow, int myCol) ;
 int isInBand(int band, int myRow, int myCol) ;



 
 // -------------------------
 // 5. Solvers
 // -------------------------
 /* 5.1 Solvers généraux */
 femSolver*           femSolverFullCreate(int size, int sizeLoc);
 femSolver*           femSolverBandCreate(int size, int sizeLoc, int band);
 femSolver*           femSolverIterativeCreate(int size, int sizeLoc);
 
 /* 5.2 Fonctions communes aux solvers */
 void                 femSolverPrintInfos(femSolver* mySolver);
 double*              femSolverEliminate(femSolver* mySolver);
 double               femSolverGet(femSolver* mySolver, int i, int j);
 int                  femSolverConverged(femSolver *mySolver);
 void                 femSolverAssemble(femSolver* mySolver, double *Aloc, double *Bloc, double *Uloc, int *mapX, int *mapY, int nLoc);
 void femSolverSystemConstrain(femSolver *mySolver, int node, double value) ;

 /* 5.3 Solvers itératifs */
 femIterativeSolver*  femIterativeSolverCreate(int size);
 void                 femIterativeSolverFree(femIterativeSolver* mySolver);
 void                 femIterativeSolverInit(femIterativeSolver* mySolver);
 void                 femIterativeSolverPrint(femIterativeSolver* mySolver);
 void                 femIterativeSolverPrintInfos(femIterativeSolver* mySolver);
 double*              femIterativeSolverEliminate(femIterativeSolver* mySolver);
 void                 femIterativeSolverAssemble(femIterativeSolver* mySolver, double *Aloc, double *Bloc, double *Uloc, int *mapX, int *mapY, int nLoc);
 double               femIterativeSolverGet(femIterativeSolver* mySolver, int i, int j);
 int                  femIterativeSolverConverged(femIterativeSolver *mySolver);
 
 // -------------------------
 // 6. Visualisation et Post-traitement
 // -------------------------
 void               trianglePlot();
 void               HexagonPlot();
 void               femFindBoundaryNodes(femGeo *theProblem, double targetY, double epsilon, char *name);
 void               femSolverSet(femSolver *mySolver, double **newA, double *newB);
 void               femElasticitySigma(femProblem *theProblem, double *sigmaXX, double *sigmaYY, double *sigmaXY);

 // -------------------------
 // 7. Renumérotation
 // -------------------------
 void                femMeshRenumber(femMesh *theMesh, femRenumType renumType);
 
 
 // -------------------------
 // 8. Gestion des Erreurs et Utilitaires
 // -------------------------
 double              femMin(double *x, int n);
 double              femMax(double *x, int n);
 void                femError(char *text, int line, char *file);
 void                femErrorScan(int test, int line, char *file);
 void                femErrorGmsh(int test, int line, char *file);
 void                femWarning(char *text, int line, char *file);
 
 
 #endif
 
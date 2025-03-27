
/*
 *  fem.c
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
#include "gmshc.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define ErrorScan(a)   femErrorScan(a,__LINE__,__FILE__)
#define ErrorGmsh(a)   femErrorGmsh(a,__LINE__,__FILE__)
#define Error(a)       femError(a,__LINE__,__FILE__)
#define Warning(a)     femWarning(a,  __LINE__, __FILE__)
#define FALSE 0 
#define TRUE  1
#define MAXNAME 256


typedef enum {FEM_TRIANGLE,FEM_QUAD,FEM_EDGE} femElementType;
typedef enum {DIRICHLET_X,DIRICHLET_Y,NEUMANN_X,NEUMANN_Y} femBoundaryType;
typedef enum {PLANAR_STRESS,PLANAR_STRAIN,AXISYM} femElasticCase;
typedef enum {FEM_FULL,FEM_BAND,FEM_ITER} femSolverType;
typedef enum {FEM_NO,FEM_XNUM,FEM_YNUM} femRenumType;




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
} femSolver;


typedef struct {
    femDomain* domain;
    femBoundaryType type; 
    double value;
} femBoundaryCondition;


typedef struct {
    double E,nu,rho,g;
    double A,B,C;
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
    // femFullSystem *system;

    femSolver *solver;
    int size;
    int sizeLoc;

} femProblem;

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

femProblem*         femElasticityCreate(femGeo* theGeometry, 
    double E, double nu, double rho, double g, femElasticCase iCase , femSolverType solverType , femRenumType renumType);
void                femElasticityFree(femProblem *theProblem);
void                femElasticityPrint(femProblem *theProblem);
void                femElasticityAddBoundaryCondition(femProblem *theProblem, char *nameDomain, femBoundaryType type, double value);
void                femElasticityAssembleElements(femProblem *theProblem);
void                femElasticityAssembleNeumann(femProblem *theProblem);
double*             femElasticitySolve(femProblem *theProblem);
double*             femElasticityForces(femProblem *theProblem);
double              femElasticityIntegrate(femProblem *theProblem, double (*f)(double x, double y));


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

femFullSystem*      femFullSystemCreate(int size);
void                femFullSystemFree(femFullSystem* mySystem);
void                femFullSystemPrint(femFullSystem* mySystem);
void                femFullSystemInit(femFullSystem* mySystem);
void                femFullSystemAlloc(femFullSystem* mySystem, int size);
double*             femFullSystemEliminate(femFullSystem* mySystem);
void                femFullSystemConstrain(femFullSystem* mySystem, int myNode, double value);



double              femMin(double *x, int n);
double              femMax(double *x, int n);
void                femError(char *text, int line, char *file);
void                femErrorScan(int test, int line, char *file);
void                femErrorGmsh(int test, int line, char *file);
void                femWarning(char *text, int line, char *file);


// Fonctions de visualisation 
void               trianglePlot();
void               HexagonPlot();
void femFindBoundaryNodes(femGeo *theProblem, double targetY, double epsilon , char *name) ;

// Different type de solvers 
femSolver*           femSolverFullCreate(int size, int sizeLoc);
femSolver*           femSolverBandCreate(int size, int sizeLoc, int band);
femSolver*           femSolverIterativeCreate(int size, int sizeLoc);
void                 femSolverFree(femSolver* mySolver);
void                 femSolverInit(femSolver* mySolver);
void                 femSolverPrint(femSolver* mySolver);
void                 femSolverPrintInfos(femSolver* mySolver);
double*              femSolverEliminate(femSolver* mySolver);
void                 femFullSystemAssemble(femFullSystem *mySystem, double *Aloc, double *Bloc, int *mapX, int *mapY, int nLocal) ; 
double               femSolverGet(femSolver* mySolver, int i, int j);
int                  femSolverConverged(femSolver *mySolver);


femBandSystem*       femBandSystemCreate(int size, int band);
void                 femBandSystemFree(femBandSystem* myBandSystem);
void                 femBandSystemInit(femBandSystem *myBand);
void                 femBandSystemPrint(femBandSystem *myBand);
void                 femBandSystemPrintInfos(femBandSystem *myBand);
double*              femBandSystemEliminate(femBandSystem *myBand);
void                 femBandSystemAssemble(femBandSystem* myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc);
double               femBandSystemGet(femBandSystem* myBandSystem, int i, int j);

femIterativeSolver*  femIterativeSolverCreate(int size);
void                 femIterativeSolverFree(femIterativeSolver* mySolver);
void                 femIterativeSolverInit(femIterativeSolver* mySolver);
void                 femIterativeSolverPrint(femIterativeSolver* mySolver);
void                 femIterativeSolverPrintInfos(femIterativeSolver* mySolver);
double*              femIterativeSolverEliminate(femIterativeSolver* mySolver);
void                 femIterativeSolverAssemble(femIterativeSolver* mySolver, double *Aloc, double *Bloc, double *Uloc, int *map, int nLoc);
double               femIterativeSolverGet(femIterativeSolver* mySolver, int i, int j);
int                  femIterativeSolverConverged(femIterativeSolver *mySolver);

void femSolverAssemble(femSolver* mySolver, double *Aloc, double *Bloc, double *Uloc,int *mapX, int *mapY, int nLoc) ;
void femMeshRenumber(femMesh *theMesh, femRenumType renumType) ;

#endif

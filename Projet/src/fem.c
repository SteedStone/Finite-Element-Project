/*
 *  fem.c
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2021 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

 #include "fem.h"

 /* ======================================================================
    1. VARIABLES GLOBALES ET INITIALISATION DE LA GÉOMÉTRIE
    ====================================================================== */
 femGeo theGeometry;
 
 
 double geoGmshSize(int dim, int tag, double x, double y, double z, double lc, void *data) {
     return theGeometry.geoSize(x, y);
 }
 
 void geoInitialize() {
     int ierr;
     theGeometry.geoSize = geoSizeDefault;
     gmshInitialize(0, NULL, 1, 0, &ierr);                         ErrorGmsh(ierr);
     gmshModelAdd("MyGeometry", &ierr);                             ErrorGmsh(ierr);
     gmshModelMeshSetSizeCallback(geoGmshSize, NULL, &ierr);         ErrorGmsh(ierr);
     theGeometry.theNodes = NULL;
     theGeometry.theElements = NULL;
     theGeometry.theEdges = NULL;
     theGeometry.nDomains = 0;
     theGeometry.theDomains = NULL;
 }

 femGeo *geoGetGeometry() {
    return &theGeometry;
}
double **A_copy = NULL;
double *B_copy  = NULL;
 
 
double geoSizeDefault(double x, double y) {
    return theGeometry.h;
}

void geoSetSizeCallback(double (*geoSize)(double x, double y)) {
    theGeometry.geoSize = geoSize;
}



 double femMin(double *x, int n) 
{
    double myMin = x[0];
    int i;
    for (i=1 ;i < n; i++) 
        myMin = fmin(myMin,x[i]);
    return myMin;
}

double femMax(double *x, int n) 
{
    double myMax = x[0];
    int i;
    for (i=1 ;i < n; i++) 
        myMax = fmax(myMax,x[i]);
    return myMax;
}


void geoMeshGenerate() {


    femGeo* theGeometry = geoGetGeometry();
    int numHexY = theGeometry->NumberOfHexagonsInY;
        double hexRadius = theGeometry->hexRadius;  // Rayon de l'hexagone
    if(theGeometry->hexa_triangles == 1) {
        
        double distance = hexRadius * 1.1; // Distance entre les lignes

        theGeometry->h = numHexY * distance - hexRadius ; 
        trianglePlot();
    } else if(theGeometry->hexa_triangles == 2) {
        theGeometry->h = numHexY * sqrt(3) * hexRadius+ hexRadius/3 ;
        HexagonPlot();
    }else {

        double w = theGeometry->LxPlate;
        double h = theGeometry->LyPlate;

        int ierr;
        double r = w/4;
        int idRect = gmshModelOccAddRectangle(0.0,0.0,0.0,w,h,-1,0.0,&ierr); 
        int idDisk = gmshModelOccAddDisk(w/2.0,h/2.0,0.0,r,r,-1,NULL,0,NULL,0,&ierr); 
        int idSlit = gmshModelOccAddRectangle(w/2.0,h/2.0-r,0.0,w,2.0*r,-1,0.0,&ierr); 
        int rect[] = {2,idRect};
        int disk[] = {2,idDisk};
        int slit[] = {2,idSlit};

        gmshModelOccCut(rect,2,disk,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr); 
        gmshModelOccCut(rect,2,slit,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr); 
        gmshModelOccSynchronize(&ierr); 

        if (theGeometry->elementType == FEM_QUAD) {
            gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
            gmshOptionSetNumber("Mesh.RecombineAll",1,&ierr);
            gmshOptionSetNumber("Mesh.Algorithm",11,&ierr);  
            gmshOptionSetNumber("Mesh.SmoothRatio", 21.5, &ierr);  
            gmshOptionSetNumber("Mesh.RecombinationAlgorithm",1.0,&ierr); 
            gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  
            gmshModelMeshGenerate(2,&ierr);  }
    
        if (theGeometry->elementType == FEM_TRIANGLE) {
            gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
            gmshModelMeshGenerate(2,&ierr);  }

        return;
    
    }

    geoSetSizeCallback(geoSize);

    int ierr;
    gmshModelOccSynchronize(&ierr);
    if (theGeometry->elementType == FEM_QUAD) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshOptionSetNumber("Mesh.RecombineAll",1,&ierr);
        gmshOptionSetNumber("Mesh.Algorithm",11,&ierr);  
        gmshOptionSetNumber("Mesh.SmoothRatio", 21.5, &ierr);  
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm",1.0,&ierr); 
        gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  
        gmshModelMeshGenerate(2,&ierr);  }

    if (theGeometry->elementType == FEM_TRIANGLE) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshModelMeshGenerate(2,&ierr);  }
    if (ierr != 0) {
        printf("Erreur lors de la génération du maillage: %d\n", ierr);
        exit(1);
    }
}
 
 
void geoMeshImport() 
{
    int ierr;
    
    /* Importing nodes */
    
    size_t nNode,n,m,*node;
    double *xyz,*trash;
    gmshModelMeshGetNodes(&node,&nNode,&xyz,&n,
                         &trash,&m,-1,-1,0,0,&ierr);          ErrorGmsh(ierr);                         
    femNodes *theNodes = malloc(sizeof(femNodes));
    theNodes->nNodes = nNode;
    theNodes->X = malloc(sizeof(double)*(theNodes->nNodes));
    theNodes->Y = malloc(sizeof(double)*(theNodes->nNodes));
    for (int i = 0; i < theNodes->nNodes; i++){
        theNodes->X[i] = xyz[3*node[i]-3];
        theNodes->Y[i] = xyz[3*node[i]-2]; }
    theGeometry.theNodes = theNodes;
    gmshFree(node);
    gmshFree(xyz);
    gmshFree(trash);
    printf("Geo     : Importing %d nodes \n",theGeometry.theNodes->nNodes);
       
    /* Importing elements */
    /* Pas super joli : a ameliorer pour eviter la triple copie */
        
    size_t nElem, *elem;
    gmshModelMeshGetElementsByType(1,&elem,&nElem,
                               &node,&nNode,-1,0,1,&ierr);    ErrorGmsh(ierr);
    femMesh *theEdges = malloc(sizeof(femMesh));
    theEdges->nLocalNode = 2;
    theEdges->nodes = theNodes;
    theEdges->nElem = nElem;  
    theEdges->elem = malloc(sizeof(int)*2*theEdges->nElem);
    for (int i = 0; i < theEdges->nElem; i++)
        for (int j = 0; j < theEdges->nLocalNode; j++)
            theEdges->elem[2*i+j] = node[2*i+j]-1;  
    theGeometry.theEdges = theEdges;
    int shiftEdges = elem[0];
    gmshFree(node);
    gmshFree(elem);
    printf("Geo     : Importing %d edges \n",theEdges->nElem);
  
    gmshModelMeshGetElementsByType(2,&elem,&nElem,
                               &node,&nNode,-1,0,1,&ierr);    ErrorGmsh(ierr);
    if (nElem != 0) {
      femMesh *theElements = malloc(sizeof(femMesh));
      theElements->nLocalNode = 3;
      theElements->nodes = theNodes;
      theElements->nElem = nElem;  
      theElements->elem = malloc(sizeof(int)*3*theElements->nElem);
      for (int i = 0; i < theElements->nElem; i++)
          for (int j = 0; j < theElements->nLocalNode; j++)
              theElements->elem[3*i+j] = node[3*i+j]-1;  
      theGeometry.theElements = theElements;
      gmshFree(node);
      gmshFree(elem);
      printf("Geo     : Importing %d triangles \n",theElements->nElem); }
    
    int nElemTriangles = nElem;
    gmshModelMeshGetElementsByType(3,&elem,&nElem,
                               &node,&nNode,-1,0,1,&ierr);    ErrorGmsh(ierr);
    if (nElem != 0 && nElemTriangles != 0)  
      Error("Cannot consider hybrid geometry with triangles and quads :-(");                       
                               
    if (nElem != 0) {
      femMesh *theElements = malloc(sizeof(femMesh));
      theElements->nLocalNode = 4;
      theElements->nodes = theNodes;
      theElements->nElem = nElem;  
      theElements->elem = malloc(sizeof(int)*4*theElements->nElem);
      for (int i = 0; i < theElements->nElem; i++)
          for (int j = 0; j < theElements->nLocalNode; j++)
              theElements->elem[4*i+j] = node[4*i+j]-1;  
      theGeometry.theElements = theElements;
      gmshFree(node);
      gmshFree(elem);
      printf("Geo     : Importing %d quads \n",theElements->nElem); }

    
    /* Importing 1D entities */
  
    int *dimTags;
    gmshModelGetEntities(&dimTags,&n,1,&ierr);        ErrorGmsh(ierr);
    theGeometry.nDomains = n/2;
    theGeometry.theDomains = malloc(sizeof(femDomain*)*n/2);
    printf("Geo     : Importing %d entities \n",theGeometry.nDomains);

    for (int i=0; i < n/2; i++) {
        int dim = dimTags[2*i+0];
        int tag = dimTags[2*i+1];
        femDomain *theDomain = malloc(sizeof(femDomain)); 
        theGeometry.theDomains[i] = theDomain;
        theDomain->mesh = theEdges;
        sprintf(theDomain->name, "Entity %d ",tag-1);
         
        int *elementType;
        size_t nElementType, **elementTags, *nElementTags, nnElementTags, **nodesTags, *nNodesTags, nnNodesTags; 
        gmshModelMeshGetElements(&elementType, &nElementType, &elementTags, &nElementTags, &nnElementTags, &nodesTags, &nNodesTags, &nnNodesTags, dim, tag, &ierr);
        theDomain->nElem = nElementTags[0];
        theDomain->elem = malloc(sizeof(int)*2*theDomain->nElem); 
        for (int j = 0; j < theDomain->nElem; j++) {
            theDomain->elem[j] = elementTags[0][j] - shiftEdges; }
        printf("Geo     : Entity %d : %d elements \n",i,theDomain->nElem);
        gmshFree(nElementTags);
        gmshFree(nNodesTags);
        gmshFree(elementTags);
        gmshFree(nodesTags);
        gmshFree(elementType); }
    gmshFree(dimTags);
 
    return;

}
 
 void geoMeshPrint() {
     femNodes *theNodes = theGeometry.theNodes;
     if (theNodes != NULL) {
         printf("Number of nodes %d \n", theNodes->nNodes);
         for (int i = 0; i < theNodes->nNodes; i++)
             printf("%6d : %14.7e %14.7e \n", i, theNodes->X[i], theNodes->Y[i]);
     }
     femMesh *theEdges = theGeometry.theEdges;
     if (theEdges != NULL) {
         printf("Number of edges %d \n", theEdges->nElem);
         int *elem = theEdges->elem;
         for (int i = 0; i < theEdges->nElem; i++)
             printf("%6d : %6d %6d \n", i, elem[2 * i], elem[2 * i + 1]);
     }
     femMesh *theElements = theGeometry.theElements;
     if (theElements != NULL) {
         if (theElements->nLocalNode == 3) {
             printf("Number of triangles %d \n", theElements->nElem);
             int *elem = theElements->elem;
             for (int i = 0; i < theElements->nElem; i++)
                 printf("%6d : %6d %6d %6d\n", i, elem[3 * i], elem[3 * i + 1], elem[3 * i + 2]);
         }
         if (theElements->nLocalNode == 4) {
             printf("Number of quads %d \n", theElements->nElem);
             int *elem = theElements->elem;
             for (int i = 0; i < theElements->nElem; i++)
                 printf("%6d : %6d %6d %6d %6d\n", i, elem[4 * i], elem[4 * i + 1], elem[4 * i + 2], elem[4 * i + 3]);
         }
     }
     int nDomains = theGeometry.nDomains;
     printf("Number of domains %d\n", nDomains);
     for (int iDomain = 0; iDomain < nDomains; iDomain++) {
         femDomain *theDomain = theGeometry.theDomains[iDomain];
         printf("  Domain : %6d \n", iDomain);
         printf("  Name : %s\n", theDomain->name);
         printf("  Number of elements : %6d\n", theDomain->nElem);
         for (int i = 0; i < theDomain->nElem; i++) {
             printf("%6d", theDomain->elem[i]);
             if ((i + 1) != theDomain->nElem && (i + 1) % 10 == 0) printf("\n");
         }
         printf("\n");
     }
 }
 
 void geoMeshWrite(const char *filename) {
     FILE *file = fopen(filename, "w");
     femNodes *theNodes = theGeometry.theNodes;
     fprintf(file, "Number of nodes %d \n", theNodes->nNodes);
     for (int i = 0; i < theNodes->nNodes; i++) {
         fprintf(file, "%6d : %14.7e %14.7e \n", i, theNodes->X[i], theNodes->Y[i]);
     }
     femMesh *theEdges = theGeometry.theEdges;
     fprintf(file, "Number of edges %d \n", theEdges->nElem);
     int *elem = theEdges->elem;
     for (int i = 0; i < theEdges->nElem; i++) {
         fprintf(file, "%6d : %6d %6d \n", i, elem[2 * i], elem[2 * i + 1]);
     }
     femMesh *theElements = theGeometry.theElements;
     if (theElements->nLocalNode == 3) {
         fprintf(file, "Number of triangles %d \n", theElements->nElem);
         elem = theElements->elem;
         for (int i = 0; i < theElements->nElem; i++) {
             fprintf(file, "%6d : %6d %6d %6d\n", i, elem[3 * i], elem[3 * i + 1], elem[3 * i + 2]);
         }
     }
     if (theElements->nLocalNode == 4) {
         fprintf(file, "Number of quads %d \n", theElements->nElem);
         elem = theElements->elem;
         for (int i = 0; i < theElements->nElem; i++) {
             fprintf(file, "%6d : %6d %6d %6d %6d\n", i, elem[4 * i], elem[4 * i + 1], elem[4 * i + 2], elem[4 * i + 3]);
         }
     }
     int nDomains = theGeometry.nDomains;
     fprintf(file, "Number of domains %d\n", nDomains);
     for (int iDomain = 0; iDomain < nDomains; iDomain++) {
         femDomain *theDomain = theGeometry.theDomains[iDomain];
         fprintf(file, "  Domain : %6d \n", iDomain);
         fprintf(file, "  Name : %s\n", theDomain->name);
         fprintf(file, "  Number of elements : %6d\n", theDomain->nElem);
         for (int i = 0; i < theDomain->nElem; i++) {
             fprintf(file, "%6d", theDomain->elem[i]);
             if ((i + 1) != theDomain->nElem && (i + 1) % 10 == 0) fprintf(file, "\n");
         }
         fprintf(file, "\n");
     }
     fclose(file);
 }
 
 void geoMeshRead(const char *filename) {
     FILE *file = fopen(filename, "r");
     int trash, *elem;
     femNodes *theNodes = malloc(sizeof(femNodes));
     theGeometry.theNodes = theNodes;
     ErrorScan(fscanf(file, "Number of nodes %d \n", &theNodes->nNodes));
     theNodes->X = malloc(sizeof(double) * (theNodes->nNodes));
     theNodes->Y = malloc(sizeof(double) * (theNodes->nNodes));
     for (int i = 0; i < theNodes->nNodes; i++) {
         ErrorScan(fscanf(file, "%d : %le %le \n", &trash, &theNodes->X[i], &theNodes->Y[i]));
     }
     femMesh *theEdges = malloc(sizeof(femMesh));
     theGeometry.theEdges = theEdges;
     theEdges->nLocalNode = 2;
     theEdges->nodes = theNodes;
     ErrorScan(fscanf(file, "Number of edges %d \n", &theEdges->nElem));
     theEdges->elem = malloc(sizeof(int) * theEdges->nLocalNode * theEdges->nElem);
     for (int i = 0; i < theEdges->nElem; i++) {
         elem = theEdges->elem;
         ErrorScan(fscanf(file, "%6d : %6d %6d \n", &trash, &elem[2 * i], &elem[2 * i + 1]));
     }
     femMesh *theElements = malloc(sizeof(femMesh));
     theGeometry.theElements = theElements;
     theElements->nLocalNode = 0;
     theElements->nodes = theNodes;
     char elementType[MAXNAME];
     ErrorScan(fscanf(file, "Number of %s %d \n", elementType, &theElements->nElem));
     if (strncasecmp(elementType, "triangles", MAXNAME) == 0) {
         theElements->nLocalNode = 3;
         theElements->elem = malloc(sizeof(int) * theElements->nLocalNode * theElements->nElem);
         for (int i = 0; i < theElements->nElem; i++) {
             elem = theElements->elem;
             ErrorScan(fscanf(file, "%6d : %6d %6d %6d \n", &trash, &elem[3 * i], &elem[3 * i + 1], &elem[3 * i + 2]));
         }
     }
     if (strncasecmp(elementType, "quads", MAXNAME) == 0) {
         theElements->nLocalNode = 4;
         theElements->elem = malloc(sizeof(int) * theElements->nLocalNode * theElements->nElem);
         for (int i = 0; i < theElements->nElem; i++) {
             elem = theElements->elem;
             ErrorScan(fscanf(file, "%6d : %6d %6d %6d %6d \n", &trash, &elem[4 * i], &elem[4 * i + 1], &elem[4 * i + 2], &elem[4 * i + 3]));
         }
     }
     ErrorScan(fscanf(file, "Number of domains %d\n", &theGeometry.nDomains));
     int nDomains = theGeometry.nDomains;
     theGeometry.theDomains = malloc(sizeof(femDomain*) * nDomains);
     for (int iDomain = 0; iDomain < nDomains; iDomain++) {
         femDomain *theDomain = malloc(sizeof(femDomain));
         theGeometry.theDomains[iDomain] = theDomain;
         theDomain->mesh = theEdges;
         ErrorScan(fscanf(file, "  Domain : %6d \n", &trash));
         ErrorScan(fscanf(file, "  Name : %[^\n]s \n", (char*)&theDomain->name));
         ErrorScan(fscanf(file, "  Number of elements : %6d\n", &theDomain->nElem));
         theDomain->elem = malloc(sizeof(int) * 2 * theDomain->nElem);
         for (int i = 0; i < theDomain->nElem; i++) {
             ErrorScan(fscanf(file, "%6d", &theDomain->elem[i]));
             if ((i + 1) != theDomain->nElem && (i + 1) % 10 == 0)
                 ErrorScan(fscanf(file, "\n"));
         }
     }
     theNodes->number = malloc(sizeof(int) * theNodes->nNodes);
     printf("Geo     : Importing %d nodes \n", theNodes->nNodes);
    
     fclose(file);
 }
 
 void geoSetDomainName(int iDomain, char *name) {
     if (iDomain >= theGeometry.nDomains)  Error("Illegal domain number");
     if (geoGetDomain(name) != -1)         Error("Cannot use the same name for two domains");
     sprintf(theGeometry.theDomains[iDomain]->name, "%s", name);
 }
 
 int geoGetDomain(char *name) {
     int theIndex = -1;
     int nDomains = theGeometry.nDomains;
     for (int iDomain = 0; iDomain < nDomains; iDomain++) {
         femDomain *theDomain = theGeometry.theDomains[iDomain];
         if (strncasecmp(name, theDomain->name, MAXNAME) == 0)
             theIndex = iDomain;
     }
     return theIndex;
 }

 void geoFinalize() {
    int ierr;
    if (theGeometry.theNodes) {
        free(theGeometry.theNodes->X);
        free(theGeometry.theNodes->Y);
        free(theGeometry.theNodes);
    }
    if (theGeometry.theElements) {
        free(theGeometry.theElements->elem);
        free(theGeometry.theElements);
    }
    if (theGeometry.theEdges) {
        free(theGeometry.theEdges->elem);
        free(theGeometry.theEdges);
    }
    free(theGeometry.theDomains);
    gmshFinalize(&ierr); ErrorGmsh(ierr);
}


void _q1c0_x(double *xsi, double *eta) 
{
    xsi[0] =  1.0;  eta[0] =  1.0;
    xsi[1] = -1.0;  eta[1] =  1.0;
    xsi[2] = -1.0;  eta[2] = -1.0;
    xsi[3] =  1.0;  eta[3] = -1.0;
}

void _q1c0_phi(double xsi, double eta, double *phi)
{
    phi[0] = (1.0 + xsi) * (1.0 + eta) / 4.0;  
    phi[1] = (1.0 - xsi) * (1.0 + eta) / 4.0;
    phi[2] = (1.0 - xsi) * (1.0 - eta) / 4.0;
    phi[3] = (1.0 + xsi) * (1.0 - eta) / 4.0;
}

void _q1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
{
    dphidxsi[0] =   (1.0 + eta) / 4.0;  
    dphidxsi[1] = - (1.0 + eta) / 4.0;
    dphidxsi[2] = - (1.0 - eta) / 4.0;
    dphidxsi[3] =   (1.0 - eta) / 4.0;
    dphideta[0] =   (1.0 + xsi) / 4.0;  
    dphideta[1] =   (1.0 - xsi) / 4.0;
    dphideta[2] = - (1.0 - xsi) / 4.0;
    dphideta[3] = - (1.0 + xsi) / 4.0;

}

void _p1c0_x(double *xsi, double *eta) 
{
    xsi[0] =  0.0;  eta[0] =  0.0;
    xsi[1] =  1.0;  eta[1] =  0.0;
    xsi[2] =  0.0;  eta[2] =  1.0;
}

void _p1c0_phi(double xsi, double eta, double *phi)
{
    phi[0] = 1 - xsi - eta;  
    phi[1] = xsi;
    phi[2] = eta;
}

void _p1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
{
    dphidxsi[0] = -1.0;  
    dphidxsi[1] =  1.0;
    dphidxsi[2] =  0.0;
    dphideta[0] = -1.0;  
    dphideta[1] =  0.0;
    dphideta[2] =  1.0;
}

void _e1c0_x(double *xsi) 
{
    xsi[0] = -1.0;  
    xsi[1] =  1.0;  
}

void _e1c0_phi(double xsi,  double *phi)
{
    phi[0] = (1 - xsi) / 2.0;  
    phi[1] = (1 + xsi) / 2.0;
}

void _e1c0_dphidx(double xsi, double *dphidxsi)
{
    dphidxsi[0] = -0.5;  
    dphidxsi[1] =  0.5;
}


 /* ======================================================================
    2. PROBLÈME D'ÉLASTICITÉ
    ====================================================================== */
 
    
    
    femProblem *femElasticityCreate(femGeo* theGeometry, 
                      double E, double nu, double rho, double g, femElasticCase iCase  , femSolverType solverType , femRenumType renumType)
    {
        femProblem *theProblem = malloc(sizeof(femProblem));
        theProblem->E   = E;
        theProblem->nu  = nu;
        theProblem->g   = g;
        theProblem->rho = rho;
        
        if (iCase == PLANAR_STRESS) {
            theProblem->A = E/(1-nu*nu);
            theProblem->B = E*nu/(1-nu*nu);
            theProblem->C = E/(2*(1+nu)); }
        else if (iCase == PLANAR_STRAIN || iCase == AXISYM) {
            theProblem->A = E*(1-nu)/((1+nu)*(1-2*nu));
            theProblem->B = E*nu/((1+nu)*(1-2*nu));
            theProblem->C = E/(2*(1+nu)); }
    
        theProblem->planarStrainStress = iCase;
        theProblem->nBoundaryConditions = 0;
        theProblem->conditions = NULL;
        
        int size = 2*theGeometry->theNodes->nNodes;
        theProblem->constrainedNodes = malloc(size*sizeof(int));
        theProblem->soluce = malloc(size*sizeof(double));
        theProblem->residuals = malloc(size*sizeof(double));
        for (int i=0; i < size; i++) {
            theProblem->constrainedNodes[i] = -1;
            theProblem->soluce[i] = 0.0;
            theProblem->residuals[i] = 0.0;}
    
    
        
        theProblem->geometry = theGeometry;  
        femMesh *theMesh = theProblem->geometry->theElements;        
    
        if (theGeometry->theElements->nLocalNode == 3) {
            theProblem->space    = femDiscreteCreate(3,FEM_TRIANGLE);
            theProblem->rule     = femIntegrationCreate(3,FEM_TRIANGLE); }
        if (theGeometry->theElements->nLocalNode == 4) {
            theProblem->space    = femDiscreteCreate(4,FEM_QUAD);
            theProblem->rule     = femIntegrationCreate(4,FEM_QUAD); }
        theProblem->spaceEdge    = femDiscreteCreate(2,FEM_EDGE);
        theProblem->ruleEdge     = femIntegrationCreate(2,FEM_EDGE); 
        // theProblem->system       = femFullSystemCreate(size);
        theProblem->size = 2*theMesh->nodes->nNodes;
        theProblem->sizeLoc = 2*theMesh->nLocalNode;
        femMeshRenumber(theMesh,renumType);
        int band ;
        
        switch (solverType) {
            case FEM_FULL : 
                    theProblem->solver = femSolverFullCreate(theProblem->size,
                                                             theProblem->sizeLoc); break;
            case FEM_BAND : 
                        
                    band = femMeshComputeBand(theMesh);
                    printf("band %d" , band) ;
                    theProblem->solver = femSolverBandCreate(theProblem->size,
                                                             theProblem->sizeLoc,band); break;
            case FEM_ITER : 
                   theProblem->solver = femSolverIterativeCreate(theProblem->size,
                                                                 theProblem->sizeLoc); break;
            default : Error("Unexpected solver option"); }
            
    
        
      
        return theProblem;
    }
    
    void femElasticityFree(femProblem *theProblem)
    {
        femFullSystemFree(theProblem->solver->solver);
        femFullSystemFree(theProblem->solver->local);
    
        femIntegrationFree(theProblem->rule);
        femDiscreteFree(theProblem->space);
        femIntegrationFree(theProblem->ruleEdge);
        femDiscreteFree(theProblem->spaceEdge);
        free(theProblem->conditions);
        free(theProblem->constrainedNodes);
        free(theProblem->soluce);
        free(theProblem->residuals);
        free(theProblem);
    }
        
    
    void femElasticityPrint(femProblem *theProblem)  
    {    
        printf("\n\n ======================================================================================= \n\n");
        printf(" Linear elasticity problem \n");
        printf("   Young modulus   E   = %14.7e [N/m2]\n",theProblem->E);
        printf("   Poisson's ratio nu  = %14.7e [-]\n",theProblem->nu);
        printf("   Density         rho = %14.7e [kg/m3]\n",theProblem->rho);
        printf("   Gravity         g   = %14.7e [m/s2]\n",theProblem->g);
        
        if (theProblem->planarStrainStress == PLANAR_STRAIN)  printf("   Planar strains formulation \n");
        if (theProblem->planarStrainStress == PLANAR_STRESS)  printf("   Planar stresses formulation \n");
        if (theProblem->planarStrainStress == AXISYM)         printf("   Axisymmetric formulation \n");
    
        printf("   Boundary conditions : \n");
        for(int i=0; i < theProblem->nBoundaryConditions; i++) {
              femBoundaryCondition *theCondition = theProblem->conditions[i];
              double value = theCondition->value;
              printf("  %20s :",theCondition->domain->name);
              if (theCondition->type==DIRICHLET_X)  printf(" imposing %9.2e as the horizontal displacement  \n",value);
              if (theCondition->type==DIRICHLET_Y)  printf(" imposing %9.2e as the vertical displacement  \n",value); 
              if (theCondition->type==NEUMANN_X)    printf(" imposing %9.2e as the horizontal force desnity \n",value); 
              if (theCondition->type==NEUMANN_Y)    printf(" imposing %9.2e as the vertical force density \n",value);}
        printf(" ======================================================================================= \n\n");
    }
    
    void femElasticityAddBoundaryCondition(femProblem *theProblem, char *nameDomain, femBoundaryType type, double value)
    {
        int iDomain = geoGetDomain(nameDomain);
        if (iDomain == -1)  Error("Undefined domain :-(");
    
        femBoundaryCondition* theBoundary = malloc(sizeof(femBoundaryCondition));
        theBoundary->domain = theProblem->geometry->theDomains[iDomain];
        theBoundary->value = value;
        theBoundary->type = type;
        theProblem->nBoundaryConditions++;
        int size = theProblem->nBoundaryConditions;
        
        if (theProblem->conditions == NULL)
            theProblem->conditions = malloc(size*sizeof(femBoundaryCondition*));
        else 
            theProblem->conditions = realloc(theProblem->conditions, size*sizeof(femBoundaryCondition*));
        theProblem->conditions[size-1] = theBoundary;
        
        int shift=-1;
        if (type == DIRICHLET_X)  shift = 0;      
        if (type == DIRICHLET_Y)  shift = 1;  
        if (shift == -1) return; 
        int *elem = theBoundary->domain->elem;
        int nElem = theBoundary->domain->nElem;
        for (int e=0; e<nElem; e++) {
            for (int i=0; i<2; i++) {
                int originalNode = theBoundary->domain->mesh->elem[2*elem[e] + i];
                int renumNode = theProblem->geometry->theNodes->number[originalNode];
                theProblem->constrainedNodes[2*renumNode + shift] = size - 1;}}    
    }
    
    
    void femElasticityAssembleElements(femProblem *theProblem) {
         femIntegration *theRule  = theProblem->rule;
         femDiscrete *theSpace    = theProblem->space;
         femGeo *theGeometry      = theProblem->geometry;
         femNodes *theNodes       = theGeometry->theNodes;
         femMesh *theMesh         = theGeometry->theElements;
         femSolver *theSolver     = theProblem->solver;
         double x[4], y[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];
         int iElem, iInteg, i, j, map[4], mapX[4], mapY[4];
         int nLocal = theMesh->nLocalNode;
         int localSize = nLocal * 2;
         double a = theProblem->A, b = theProblem->B, c = theProblem->C;
         double rho = theProblem->rho, g = theProblem->g;
         double *Aloc = theSolver->local->A[0];
         double *Bloc = theSolver->local->B;
         for (iElem = 0; iElem < theMesh->nElem; iElem++) {
             for (j = 0; j < nLocal; j++) {
                 map[j] = theMesh->elem[iElem * nLocal + j];
                 x[j] = theNodes->X[map[j]];
                 y[j] = theNodes->Y[map[j]];
                 map[j] = theMesh->nodes->number[map[j]];
                 mapX[j] = 2 * map[j];
                 mapY[j] = 2 * map[j] + 1;
             }
             memset(Aloc, 0, localSize * localSize * sizeof(double));
             memset(Bloc, 0, localSize * sizeof(double));
             for (iInteg = 0; iInteg < theRule->n; iInteg++) {
                 double xsi = theRule->xsi[iInteg];
                 double eta = theRule->eta[iInteg];
                 double weight = theRule->weight[iInteg];
                 femDiscretePhi2(theSpace, xsi, eta, phi);
                 femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);
                 double dxdxsi = 0.0, dxdeta = 0.0, dydxsi = 0.0, dydeta = 0.0;
                 for (i = 0; i < theSpace->n; i++) {
                     dxdxsi += x[i] * dphidxsi[i];
                     dxdeta += x[i] * dphideta[i];
                     dydxsi += y[i] * dphidxsi[i];
                     dydeta += y[i] * dphideta[i];
                 }
                 double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
                 for (i = 0; i < theSpace->n; i++) {
                     dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
                     dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
                 }
                 double weightedJac = jac * weight;
                 for (i = 0; i < theSpace->n; i++) {
                     for (j = 0; j < theSpace->n; j++) {
                         int iX = 2 * i, iY = 2 * i + 1;
                         int jX = 2 * j, jY = 2 * j + 1;
                         Aloc[iX * localSize + jX] += (dphidx[i] * a * dphidx[j] + dphidy[i] * c * dphidy[j]) * weightedJac;
                         Aloc[iX * localSize + jY] += (dphidx[i] * b * dphidy[j] + dphidy[i] * c * dphidx[j]) * weightedJac;
                         Aloc[iY * localSize + jX] += (dphidy[i] * b * dphidx[j] + dphidx[i] * c * dphidy[j]) * weightedJac;
                         Aloc[iY * localSize + jY] += (dphidy[i] * a * dphidy[j] + dphidx[i] * c * dphidx[j]) * weightedJac;
                     }
                     int iY = 2 * i + 1;
                     Bloc[iY] -= phi[i] * g * rho * weightedJac;
                 }
             }
             femSolverAssemble(theSolver, Aloc, Bloc, NULL, mapX, mapY, nLocal);
         }
     }
     
     void femElasticityAssembleNeumann(femProblem *theProblem) {
         femFullSystem *theSystem = theProblem->solver->solver;
         femIntegration *theRule = theProblem->ruleEdge;
         femDiscrete *theSpace = theProblem->spaceEdge;
         femGeo *theGeometry = theProblem->geometry;
         femNodes *theNodes = theGeometry->theNodes;
         femMesh *theEdges = theGeometry->theEdges;
         double x[2], y[2], phi[2];
         int iBnd, iElem, iInteg, iEdge, i, j, shift, map[2], mapU[2];
         int nLocal = theSpace->n;
         int *number = theNodes->number;
         double *B = theSystem->B;
         for (iBnd = 0; iBnd < theProblem->nBoundaryConditions; iBnd++) {
             femBoundaryCondition *theCondition = theProblem->conditions[iBnd];
             femBoundaryType type = theCondition->type;
             femDomain *domain = theCondition->domain;
             double value = theCondition->value;
             femMesh *theMesh = theProblem->geometry->theElements;
             if (type == DIRICHLET_X || type == DIRICHLET_Y) continue;
             shift = (type == NEUMANN_X) ? 0 : 1;
             for (iEdge = 0; iEdge < domain->nElem; iEdge++) {
                 iElem = domain->elem[iEdge];
                 for (j = 0; j < nLocal; j++) {
                     map[j] = theEdges->elem[iElem * nLocal + j];
                     x[j] = theNodes->X[map[j]];
                     y[j] = theNodes->Y[map[j]];
                     map[j] = theMesh->nodes->number[map[j]];
                     mapU[j] = nLocal * map[j] + shift;
                 }
                 double dx = x[1] - x[0], dy = y[1] - y[0];
                 double length = sqrt(dx * dx + dy * dy);
                 double jac = length / 2;
                 for (iInteg = 0; iInteg < theRule->n; iInteg++) {
                     double xsi = theRule->xsi[iInteg];
                     double weight = theRule->weight[iInteg];
                     femDiscretePhi(theSpace, xsi, phi);
                     for (i = 0; i < theSpace->n; i++)
                         B[mapU[i]] += phi[i] * value * jac * weight;
                 }
             }
         }
     }
     
     
    double femElasticityIntegrate(femProblem *theProblem, double (*f)(double x, double y)){
        femIntegration *theRule = theProblem->rule;
        femGeo         *theGeometry = theProblem->geometry;
        femNodes       *theNodes = theGeometry->theNodes;
        femMesh        *theMesh = theGeometry->theElements;
        femDiscrete    *theSpace = theProblem->space;
    
        double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
        int iElem,iInteg,i,map[4];
        int nLocal = theMesh->nLocalNode;
        double value = 0.0;
        for (iElem = 0; iElem < theMesh->nElem; iElem++) {
            for (i=0; i < nLocal; i++) {
                map[i]  = theMesh->elem[iElem*nLocal+i];
                x[i]    = theNodes->X[map[i]];
                y[i]    = theNodes->Y[map[i]];} 
            for (iInteg=0; iInteg < theRule->n; iInteg++) {    
                double xsi    = theRule->xsi[iInteg];
                double eta    = theRule->eta[iInteg];
                double weight = theRule->weight[iInteg];  
                femDiscretePhi2(theProblem->space,xsi,eta,phi);
                femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
                double dxdxsi = 0.0;
                double dxdeta = 0.0;
                double dydxsi = 0.0; 
                double dydeta = 0.0;
                for (i = 0; i < theSpace->n; i++) {  
                    dxdxsi += x[i]*dphidxsi[i];       
                    dxdeta += x[i]*dphideta[i];   
                    dydxsi += y[i]*dphidxsi[i];   
                    dydeta += y[i]*dphideta[i]; }
                double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
                for (i = 0; i < theProblem->space->n; i++) {    
                    value += phi[i] * f(x[i],y[i]) * jac * weight; }}}
        return value;
    
    }

     double *femElasticitySolve(femProblem *theProblem) {
         femFullSystem *theSystem = theProblem->solver->solver;
         femSolver *theSolver = theProblem->solver;
         femElasticityAssembleElements(theProblem);
         femElasticityAssembleNeumann(theProblem);
         int size = theSystem->size;
         // Copies des matrices pour calculer les résidus
         if (A_copy == NULL) {
             A_copy = (double **)malloc(sizeof(double *) * size);
             for (int i = 0; i < size; i++) {
                 A_copy[i] = (double *)malloc(sizeof(double) * size);
             }
         }
         if (B_copy == NULL) {
             B_copy = (double *)malloc(sizeof(double) * size);
         }
         for (int i = 0; i < size; i++) {
             for (int j = 0; j < size; j++)
                 A_copy[i][j] = theSystem->A[i][j];
             B_copy[i] = theSystem->B[i];
         }
         int *theConstrainedNodes = theProblem->constrainedNodes;
         femGeo *theGeometry      = theProblem->geometry;
         femNodes *theNodes       = theGeometry->theNodes;
         int *number = theNodes->number;

         for (int i = 0; i < size; i++) {
             if (theConstrainedNodes[i] != -1) {
                 double value = theProblem->conditions[theConstrainedNodes[i]]->value;
                 femSolverSystemConstrain(theSolver, i, value);
             }
         }
         femSolverEliminate(theProblem->solver);
         for (int i = 0; i < theProblem->size / 2; i++) {
        
             theProblem->soluce[2 * i] += theSystem->B[2 * i + 0];
             theProblem->soluce[2 * i + 1] += theSystem->B[2 * i + 1];
         }
         return theProblem->soluce;
     }
     
     double *femElasticityForces(femProblem *theProblem) {
         double *residuals = theProblem->residuals;
         double *soluce = theProblem->soluce;
         femSolver *theSolver = theProblem->solver;
         int size = theProblem->size;
         if (residuals == NULL) {
             residuals = (double *)malloc(sizeof(double) * size);
         }
         for (int i = 0; i < size; i++) {
             residuals[i] = 0.0;
         }
         femSolverGetResidual(theSolver, residuals, soluce);

         for (int i = 0; i < size; i++) { free(A_copy[i]); A_copy[i] = NULL; }
         free(A_copy); free(B_copy);
         A_copy = NULL; B_copy = NULL;
         return residuals;
     }

     void femSolverGetResidual(femSolver *mySolver, double *residuals, double *theSoluce)
{
    switch (mySolver->type)
    {
        case FEM_FULL : femFullSystemGetResidual((femFullSystem *) mySolver->solver, mySolver->size, residuals, theSoluce); break;
        case FEM_BAND : femBandSystemGetResidual((femBandSystem *) mySolver->solver, mySolver->size, residuals, theSoluce); break;
        default :       Error("Unexpected solver type");
    }
}

void femFullSystemGetResidual(femFullSystem *mySystem, int size, double *residuals, double *theSoluce)
{
    int i, j;

    

    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++) { residuals[i] += A_copy[i][j] * theSoluce[j]; }
        residuals[i] -= B_copy[i];
    }
}

void femBandSystemGetResidual(femBandSystem *mySystem, int size, double *residuals, double *theSoluce)
{
    double  A_ij;
    int band, start, end, i, j;

   
    band = mySystem->band;

    for (i = 0; i < size; i++)
    {
        start = (i - band > 0) ? i - band : 0;
        end = (i + band < size) ? i + band : size;
        for (j = start; j < end; j++)
        {
            A_ij = (j >= i) ? femBandSystemGetA_Entry(mySystem, i, j) : femBandSystemGetA_Entry(mySystem, j, i);
            residuals[i] += A_ij * theSoluce[j];
        }
        residuals[i] -= B_copy[i];
    }
}
     

 
 /* ======================================================================
    3. DISCRÉTISATION ET INTÉGRATION
    ====================================================================== */
 /* Définitions statiques des points et poids de Gauss pour quadrangles, triangles et arêtes */
 static const double _gaussQuad4Xsi[4]    = {-0.577350269189626, -0.577350269189626,  0.577350269189626,  0.577350269189626};
 static const double _gaussQuad4Eta[4]    = { 0.577350269189626, -0.577350269189626, -0.577350269189626,  0.577350269189626};
 static const double _gaussQuad4Weight[4] = { 1.0, 1.0, 1.0, 1.0};
 static const double _gaussTri3Xsi[3]     = { 0.166666666666667,  0.666666666666667,  0.166666666666667};
 static const double _gaussTri3Eta[3]     = { 0.166666666666667,  0.166666666666667,  0.666666666666667};
 static const double _gaussTri3Weight[3]  = { 0.166666666666667,  0.166666666666667,  0.166666666666667};
 static const double _gaussEdge2Xsi[2]    = { 0.577350269189626, -0.577350269189626};
 static const double _gaussEdge2Weight[2] = { 1.0, 1.0};
 
 femIntegration *femIntegrationCreate(int n, femElementType type) {
     femIntegration *theRule = malloc(sizeof(femIntegration));
     if (type == FEM_QUAD && n == 4) {
         theRule->n      = 4;
         theRule->xsi    = _gaussQuad4Xsi;
         theRule->eta    = _gaussQuad4Eta;
         theRule->weight = _gaussQuad4Weight;
     } else if (type == FEM_TRIANGLE && n == 3) {
         theRule->n      = 3;
         theRule->xsi    = _gaussTri3Xsi;
         theRule->eta    = _gaussTri3Eta;
         theRule->weight = _gaussTri3Weight;
     } else if (type == FEM_EDGE && n == 2) {
         theRule->n      = 2;
         theRule->xsi    = _gaussEdge2Xsi;
         theRule->eta    = NULL;
         theRule->weight = _gaussEdge2Weight;
     } else
         Error("Cannot create such an integration rule !");
     return theRule;
 }
 
 void femIntegrationFree(femIntegration *theRule) {
     free(theRule);
 }
 
 /* Fonctions pour la création, l'utilisation et l'affichage des espaces discrets */
 femDiscrete *femDiscreteCreate(int n, femElementType type) {
     femDiscrete *theSpace = malloc(sizeof(femDiscrete));
     theSpace->type = type;
     theSpace->n = 0;
     theSpace->x = NULL;
     theSpace->phi = NULL;
     theSpace->dphidx = NULL;
     theSpace->x2 = NULL;
     theSpace->phi2 = NULL;
     theSpace->dphi2dx = NULL;
     if (type == FEM_QUAD && n == 4) {
         theSpace->n = 4;
         theSpace->x2 = _q1c0_x;
         theSpace->phi2 = _q1c0_phi;
         theSpace->dphi2dx = _q1c0_dphidx;
     } else if (type == FEM_TRIANGLE && n == 3) {
         theSpace->n = 3;
         theSpace->x2 = _p1c0_x;
         theSpace->phi2 = _p1c0_phi;
         theSpace->dphi2dx = _p1c0_dphidx;
     } else if (type == FEM_EDGE && n == 2) {
         theSpace->n = 2;
         theSpace->x = _e1c0_x;
         theSpace->phi = _e1c0_phi;
         theSpace->dphidx = _e1c0_dphidx;
     } else
         Error("Cannot create such a discrete space !");
     return theSpace;
 }
 
 void femDiscreteFree(femDiscrete *theSpace) {
     free(theSpace);
 }

 void femDiscretePrint(femDiscrete *mySpace) {
    int i, j;
    int n = mySpace->n;
    double xsi[4], eta[4], phi[4], dphidxsi[4], dphideta[4];
    if (mySpace->type == FEM_EDGE) {
        femDiscreteXsi(mySpace, xsi);
        for (i = 0; i < n; i++) {
            femDiscretePhi(mySpace, xsi[i], phi);
            femDiscreteDphi(mySpace, xsi[i], dphidxsi);
            for (j = 0; j < n; j++) {
                printf("(xsi=%+.1f) : ", xsi[i]);
                printf(" phi(%d)=%+.1f", j, phi[j]);
                printf("   dphidxsi(%d)=%+.1f \n", j, dphidxsi[j]);
            }
            printf(" \n");
        }
    }
    if (mySpace->type == FEM_QUAD || mySpace->type == FEM_TRIANGLE) {
        femDiscreteXsi2(mySpace, xsi, eta);
        for (i = 0; i < n; i++) {
            femDiscretePhi2(mySpace, xsi[i], eta[i], phi);
            femDiscreteDphi2(mySpace, xsi[i], eta[i], dphidxsi, dphideta);
            for (j = 0; j < n; j++) {
                printf("(xsi=%+.1f,eta=%+.1f) : ", xsi[i], eta[i]);
                printf(" phi(%d)=%+.1f", j, phi[j]);
                printf("   dphidxsi(%d)=%+.1f", j, dphidxsi[j]);
                printf("   dphideta(%d)=%+.1f \n", j, dphideta[j]);
            }
            printf(" \n");
        }
    }
}
 
 void femDiscreteXsi2(femDiscrete *mySpace, double *xsi, double *eta) {
     mySpace->x2(xsi, eta);
 }
 
 void femDiscretePhi2(femDiscrete *mySpace, double xsi, double eta, double *phi) {
     mySpace->phi2(xsi, eta, phi);
 }
 
 void femDiscreteDphi2(femDiscrete *mySpace, double xsi, double eta, double *dphidxsi, double *dphideta) {
     mySpace->dphi2dx(xsi, eta, dphidxsi, dphideta);
 }
 
 void femDiscreteXsi(femDiscrete *mySpace, double *xsi) {
     mySpace->x(xsi);
 }
 
 void femDiscretePhi(femDiscrete *mySpace, double xsi, double *phi) {
     mySpace->phi(xsi, phi);
 }
 
 void femDiscreteDphi(femDiscrete *mySpace, double xsi, double *dphidxsi) {
     mySpace->dphidx(xsi, dphidxsi);
 }
 

 
 /* ======================================================================
    4. SYSTÈMES LINÉAIRES : CRÉATION, ASSEMBLAGE, SOLUTION ET CONTRAINTE
    ====================================================================== */
 /* 4.1 Système complet (full system) */

 femFullSystem *femFullSystemCreate(int size) {
     femFullSystem *theSystem = malloc(sizeof(femFullSystem));
     femFullSystemAlloc(theSystem, size);
     femFullSystemInit(theSystem);
     return theSystem;
 }
 
 void femFullSystemFree(femFullSystem *theSystem) {
     free(theSystem->A);
     free(theSystem->B);
     free(theSystem);
 }

 void femFullSystemPrint(femFullSystem *mySystem) {
    double **A = mySystem->A;
    double *B = mySystem->B;
    int i, j, size = mySystem->size;
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++)
            printf(" %+.1e", A[i][j]);
        printf(" :  %+.1e \n", B[i]);
    }
}

 
void femFullSystemInit(femFullSystem *mySystem) {
    int i, size = mySystem->size;
    for (i = 0; i < size * (size + 1); i++)
        mySystem->B[i] = 0;
}
 
 void femFullSystemAlloc(femFullSystem *mySystem, int size) {
     int i;
     double *elem = malloc(sizeof(double) * size * (size + 1));
     mySystem->A = malloc(sizeof(double*) * size);
     mySystem->B = elem;
     mySystem->A[0] = elem + size;
     mySystem->size = size;
     for (i = 1; i < size; i++)
         mySystem->A[i] = mySystem->A[i - 1] + size;
 }


 
void luDecomposition(double **A, int *P, int size) {
    for (int k = 0; k < size; k++) {
        int maxIndex = k;
        for (int i = k + 1; i < size; i++) {
            if (fabs(A[i][k]) > fabs(A[maxIndex][k])) {
                maxIndex = i;
            }
        }
        
        if (fabs(A[maxIndex][k]) < 1e-16) {
            printf("Singular matrix detected at column %d\n", k);
            exit(EXIT_FAILURE);
        }
        
        if (maxIndex != k) {
            double *tempRow = A[k];
            A[k] = A[maxIndex];
            A[maxIndex] = tempRow;
            int tempP = P[k];
            P[k] = P[maxIndex];
            P[maxIndex] = tempP;
        }
        
        for (int i = k + 1; i < size; i++) {
            A[i][k] /= A[k][k];
            for (int j = k + 1; j < size; j++) {
                A[i][j] -= A[i][k] * A[k][j];
            }
        }
    }
}

void luSolve(double **A, int *P, double *B, int size) {
    double *Y = (double *)malloc(size * sizeof(double));
    
    for (int i = 0; i < size; i++) {
        Y[i] = B[P[i]];
        for (int j = 0; j < i; j++) {
            Y[i] -= A[i][j] * Y[j];
        }
    }
    
    for (int i = size - 1; i >= 0; i--) {
        for (int j = i + 1; j < size; j++) {
            Y[i] -= A[i][j] * B[j];
        }
        B[i] = Y[i] / A[i][i];
    }
    
    free(Y);
}
 
 double* femFullSystemEliminate(femFullSystem *mySystem) {
     double **A = mySystem->A;
     double *B = mySystem->B;
     int size = mySystem->size;
     int *P = (int *)malloc(size * sizeof(int));
     for (int i = 0; i < size; i++) P[i] = i;
     luDecomposition(A, P, size);
     luSolve(A, P, B, size);
     free(P);
     return B;
 }
 
 void femFullSystemConstrain(femFullSystem *mySystem, int myNode, double myValue , int size) {
     double **A = mySystem->A;
     double *B = mySystem->B;
     int i;
     for (i = 0; i < size; i++) {
         B[i] -= myValue * A[i][myNode];
         A[i][myNode] = 0;
     }
     for (i = 0; i < size; i++)
         A[myNode][i] = 0;
     A[myNode][myNode] = 1;
     B[myNode] = myValue;
 }

 double femFullSystemGet(femFullSystem *myFullSystem, int myRow, int myCol) {
    return myFullSystem->A[myRow][myCol];
}


void femFullSystemAssemble(femFullSystem *mySystem, 
    double *Aloc, double *Bloc, 
    int *mapX, int *mapY, int nLocal)
{
    int localSize = 2 * nLocal;
    for (int i = 0; i < nLocal; i++) {
        for (int j = 0; j < nLocal; j++) {

            mySystem->A[ mapX[i] ][ mapX[j] ] += Aloc[(2*i) * localSize + (2*j)];
            // Le bloc XY
            mySystem->A[ mapX[i] ][ mapY[j] ] += Aloc[(2*i) * localSize + (2*j+1)];
            // Le bloc YX
            mySystem->A[ mapY[i] ][ mapX[j] ] += Aloc[(2*i+1) * localSize + (2*j)];
            // Le bloc YY
            mySystem->A[ mapY[i] ][ mapY[j] ] += Aloc[(2*i+1) * localSize + (2*j+1)];
            }
        // On assemble le vecteur local de force
        mySystem->B[ mapX[i] ] += Bloc[2*i];
        mySystem->B[ mapY[i] ] += Bloc[2*i+1];
    }
}


 
void femFullSystemPrintInfos(femFullSystem *mySystem) {
    int  size = mySystem->size;
    printf(" \n");
    printf("    Full Gaussian elimination \n");
    printf("    Storage informations \n");
    printf("    Matrix size      : %8d\n",size);
    printf("    Bytes required   : %8d\n",(int)sizeof(double)*size*(size+1));     
 }


/* 4.2 Système bande (band system) */


int femMeshComputeBand(femMesh *theMesh)
{
    // Strip : BEGIN
    int myBand = 0;
    int maxNum, minNum, nodeNum, elemNum;

    for (int iElem = 0; iElem < theMesh->nElem; iElem++)
    {   
        maxNum = INT_MIN;
        minNum = INT_MAX;

        for (int j = 0; j < theMesh->nLocalNode; j++)
        {
            elemNum = theMesh->elem[iElem * theMesh->nLocalNode + j];
            nodeNum = theMesh->nodes->number[elemNum];

            maxNum = (nodeNum > maxNum) ? nodeNum : maxNum;
            minNum = (nodeNum < minNum) ? nodeNum : minNum;
        }

        if (myBand < maxNum - minNum) { myBand = maxNum - minNum; }
    }
    
    return 2 * (myBand + 1);
    // Strip : END
}


#ifndef NOBANDASSEMBLE

void femBandSystemAssemble(femBandSystem *myBandSystem, 
    double *Aloc, double *Bloc, 
    int *mapX, int *mapY, int nLoc)
{
    
int i, j;
int currRow, currCol;
int localSize = 2 * nLoc; // deux degrés de liberté par noeud local

// Assemblage de la matrice de raideur (stiffness matrix)
for (i = 0; i < nLoc; i++)
{
// Bloc X-X
currRow = mapX[i];
for (j = 0; j < nLoc; j++)
{
currCol = mapX[j];
if (currRow <= currCol)
{
// La case (2*i,2*j) dans Aloc correspond au bloc X-X
myBandSystem->A[currRow][currCol] += Aloc[(2 * i) * localSize + (2 * j)];
}
}
// Bloc X-Y
currRow = mapX[i];
for (j = 0; j < nLoc; j++)
{
currCol = mapY[j];
if (currRow <= currCol)
{
// La case (2*i,2*j+1) dans Aloc correspond au bloc X-Y
myBandSystem->A[currRow][currCol] += Aloc[(2 * i) * localSize + (2 * j + 1)];
}
}
// Bloc Y-X
currRow = mapY[i];
for (j = 0; j < nLoc; j++)
{
currCol = mapX[j];
if (currRow <= currCol)
{
// La case (2*i+1,2*j) dans Aloc correspond au bloc Y-X
myBandSystem->A[currRow][currCol] += Aloc[(2 * i + 1) * localSize + (2 * j)];
}
}
// Bloc Y-Y
currRow = mapY[i];
for (j = 0; j < nLoc; j++)
{
currCol = mapY[j];
if (currRow <= currCol)
{
// La case (2*i+1,2*j+1) dans Aloc correspond au bloc Y-Y
myBandSystem->A[currRow][currCol] += Aloc[(2 * i + 1) * localSize + (2 * j + 1)];
}
}

// Assemblage du vecteur second membre
myBandSystem->B[ mapX[i] ] += Bloc[2 * i];
myBandSystem->B[ mapY[i] ] += Bloc[2 * i + 1];
}
}


#endif
#ifndef NOBANDELIMINATE

double  *femBandSystemEliminate(femBandSystem *myBand)
{
    double  **A, *B, factor;
    int     i, j, k, jend, size, band;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;
    band = myBand->band;

    // Strip : BEGIN

    for (k = 0; k < size; k++)
    {
        if (fabs(A[k][k]) <= 1e-8) { Error("Cannot eliminate with such a pivot.\n"); }
        jend = (k + band < size) ? k + band : size;
        for (i = k + 1; i < jend; i++)
        {
            factor = A[k][i] / A[k][k];
            for (j = i ; j < jend; j++) { A[i][j] -= factor * A[k][j]; }
            B[i] -= factor * B[k];
        }    
    }
    
    for (i = size - 1; i >= 0 ; i--)
    {
        factor = 0;
        jend = (i + band < size) ? i + band : size;
        for (j = i + 1 ; j < jend; j++) { factor += A[i][j] * B[j]; }
        B[i] = ( B[i] - factor) / A[i][i];
    }
    // Strip : END

    return myBand->B;
}


femBandSystem *femBandSystemCreate(int size, int band)
{
    femBandSystem *myBandSystem = malloc(sizeof(femBandSystem));
    myBandSystem->B = malloc(sizeof(double)*size*(band+1));
    myBandSystem->A = malloc(sizeof(double*)*size);        
    myBandSystem->size = size;
    myBandSystem->band = band;
    myBandSystem->A[0] = myBandSystem->B + size;
    int i;
    for (i=1 ; i < size ; i++) 
        myBandSystem->A[i] = myBandSystem->A[i-1] + band - 1;
    femBandSystemInit(myBandSystem);
    return(myBandSystem);
}
 
void femBandSystemFree(femBandSystem *myBandSystem)
{
    free(myBandSystem->B);
    free(myBandSystem->A); 
    free(myBandSystem);
}
 
void femBandSystemInit(femBandSystem *myBandSystem)
{
    int i;
    int size = myBandSystem->size;
    int band = myBandSystem->band;
    for (i=0 ; i < size*(band+1) ; i++) 
        myBandSystem->B[i] = 0;        
}
 
void femBandSystemPrint(femBandSystem *myBand)
{
    double  **A, *B;
    int     i, j, band, size;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;
    band = myBand->band;

    for (i=0; i < size; i++) {
        for (j=i; j < i+band; j++)
            if (A[i][j] == 0) printf("         ");   
            else              printf(" %+.1e",A[i][j]);
        printf(" :  %+.1e \n",B[i]); }
}
  
void femBandSystemPrintInfos(femBandSystem *myBand)
{
    int size = myBand->size;
    int band = myBand->band;
    printf(" \n");
    printf("    Banded Gaussian elimination \n");
    printf("    Storage informations \n");
    printf("    Matrix size      : %8d\n",size);
    printf("    Matrix band      : %8d\n",band);
    printf("    Bytes required   : %8d\n",(int)sizeof(double)*size*(band+1));     
}


double femBandSystemGet(femBandSystem* myBandSystem, int myRow, int myCol)
{
    double value = 0;
    if (myCol >= myRow && myCol < myRow+myBandSystem->band)  value = myBandSystem->A[myRow][myCol]; 
    return(value);
}

int isInBand(int band, int myRow, int myCol) { return myCol >= myRow && myCol < myRow + band; }

double femBandSystemGetA_Entry(femBandSystem *mySystem, int myRow, int myCol) { return (isInBand(mySystem->band, myRow, myCol)) ? A_copy[myRow][myCol] : 0.0; }

void femBandSystemConstrain(femBandSystem *mySystem, int myNode, double myValue, int size)
{
    double **A, *B, A_entry;
    int i, band;

    A = mySystem->A;
    B = mySystem->B;
    band = mySystem->band;

    for (i = 0; i < size; i++)
    {
        A_entry = (myNode >= i) ? femBandSystemGetA_Entry(mySystem, i, myNode) : femBandSystemGetA_Entry(mySystem, myNode, i);
        if (A_entry != 0.0)
        {
            B[i] -= myValue * A_entry;
            if (myNode >= i) { A[i][myNode] = 0; }
        }
    }
    for (int i = 0; i < size; i++)
    {
        if (femBandSystemGetA_Entry(mySystem, myNode, i) != 0.0) { A[myNode][i] = 0.0; }
    }
    
    A[myNode][myNode] = 1.0;
    B[myNode] = myValue;
}

 
 /* ======================================================================
    5. SOLVERS : CRÉATION ET ASSEMBLAGE
    ====================================================================== */
/* 5.1 Solvers généraux */
femSolver *femSolverCreate(int sizeLoc){
        femSolver *mySolver = malloc(sizeof(femSolver));
        mySolver->local = femFullSystemCreate(sizeLoc);
        return (mySolver);
}

 femSolver *femSolverFullCreate(int size, int sizeLoc) {
     femSolver *mySolver = femSolverCreate(sizeLoc);
     mySolver->type = FEM_FULL;
     mySolver->solver = (void *)femFullSystemCreate(size);
     mySolver->size = size;
     return mySolver;
 }
 
 femSolver *femSolverBandCreate(int size, int sizeLoc, int band) {
     femSolver *mySolver = femSolverCreate(sizeLoc);
     mySolver->type = FEM_BAND;
     mySolver->solver = (void *)femBandSystemCreate(size, band);
     mySolver->size = size;
     return mySolver;
 }
 
 femSolver *femSolverIterativeCreate(int size, int sizeLoc) {
     femSolver *mySolver = femSolverCreate(sizeLoc);
     mySolver->type = FEM_ITER;
     mySolver->solver = (void *)femIterativeSolverCreate(size);
     return mySolver;
 }
  /* 5.2 Fonctions communes aux solvers */
 
 void femSolverAssemble(femSolver* mySolver, double *Aloc, double *Bloc, double *Uloc, int *mapX, int *mapY, int nLoc) {
     switch (mySolver->type) {
         case FEM_FULL:
             femFullSystemAssemble((femFullSystem *)mySolver->solver, Aloc, Bloc, mapX, mapY, nLoc);
             break;
         case FEM_BAND:
             femBandSystemAssemble((femBandSystem *)mySolver->solver, Aloc, Bloc, mapX, mapY, nLoc);
             break;
         case FEM_ITER:
             femIterativeSolverAssemble((femIterativeSolver *)mySolver->solver, Aloc, Bloc, Uloc, mapX, mapY, nLoc);
             break;
         default:
             Error("Unexpected solver type");
     }
 }
 
 double* femSolverEliminate(femSolver *mySolver) {
     double *soluce;
     switch (mySolver->type) {
         case FEM_FULL:
             soluce = femFullSystemEliminate((femFullSystem *)mySolver->solver);
             break;
         case FEM_BAND:
             soluce = femBandSystemEliminate((femBandSystem *)mySolver->solver);
             break;
         case FEM_ITER:
             soluce = femIterativeSolverEliminate((femIterativeSolver *)mySolver->solver);
             break;
         default:
             Error("Unexpected solver type");
     }
     return soluce;
 }
 
 double femSolverGet(femSolver *mySolver, int i, int j) {
    double value = 0;
    switch (mySolver->type) {
        case FEM_FULL:
            value = femFullSystemGet((femFullSystem *)mySolver->solver, i, j);
            break;
        case FEM_BAND:
            value = femBandSystemGet((femBandSystem *)mySolver->solver, i, j);
            break;
        default:
            Error("Unexpected solver type");
    }
    return value;
}

void femSolverPrintInfos(femSolver *mySolver) {
    switch (mySolver->type) {
        case FEM_FULL:
            femFullSystemPrintInfos((femFullSystem *)mySolver->solver);
            break;
        case FEM_BAND:
            femBandSystemPrintInfos((femBandSystem *)mySolver->solver);
            break;
        default:
            Error("Unexpected solver type");
    }
}

int femSolverConverged(femSolver *mySolver)
{
    int  testConvergence;
    switch (mySolver->type) {
        case FEM_FULL : testConvergence = 1; break;
        case FEM_BAND : testConvergence = 1; break;
        // case FEM_ITER : testConvergence = femIterativeSolverConverged((femIterativeSolver *)mySolver->solver); break;
        default : Error("Unexpected solver type"); }
    return(testConvergence);
}

void femSolverSystemConstrain(femSolver *mySolver, int node, double value)
{
    switch (mySolver->type)
    {
        case FEM_FULL : femFullSystemConstrain((femFullSystem *) mySolver->solver, node, value, mySolver->size); break;
        case FEM_BAND : femBandSystemConstrain((femBandSystem *) mySolver->solver, node, value, mySolver->size); break;
        default :       Error("Unexpected solver type");
    }
}


 /* 5.3 Solvers itératifs */

 

femIterativeSolver *femIterativeSolverCreate(int size)
{
    femIterativeSolver *mySolver = malloc(sizeof(femIterativeSolver));
    mySolver->R = malloc(sizeof(double)*size*4);      
    mySolver->D = mySolver->R + size;       
    mySolver->S = mySolver->R + size*2;       
    mySolver->X = mySolver->R + size*3;       
    mySolver->size = size;
    femIterativeSolverInit(mySolver);
    return(mySolver);
}

void femIterativeSolverFree(femIterativeSolver *mySolver)
{
    free(mySolver->R);
    free(mySolver);
}

void femIterativeSolverInit(femIterativeSolver *mySolver)
{
    int i;
    mySolver->iter = 0;
    mySolver->error = 10.0e+12;
    for (i=0 ; i < mySolver->size*4 ; i++) 
        mySolver->R[i] = 0;        
}
 
void femIterativeSolverPrint(femIterativeSolver *mySolver)
{
    double  *R;
    int     i, size;
    R    = mySolver->R;
    size = mySolver->size;

    for (i=0; i < size; i++) {
        printf("%d :  %+.1e \n",i,R[i]); }
}

void femIterativeSolverPrintInfos(femIterativeSolver *mySolver)
{
    if (mySolver->iter == 1)     printf("\n    Iterative solver \n");
    printf("    Iteration %4d : %14.7e\n",mySolver->iter,mySolver->error);
}

int femIterativeSolverConverged(femIterativeSolver *mySolver)
{
    int  testConvergence = 0;
    if (mySolver->iter  > 3000)     testConvergence = -1;
    if (mySolver->error < 10.0e-6)  testConvergence = 1;
    return(testConvergence);
}

void femIterativeSolverAssemble(femIterativeSolver* mySolver, 
                                double *Aloc, 
                                double *Bloc, 
                                double *Uloc, 
                                int *mapX, 
                                int *mapY, 
                                int nLoc)
{
    int i, j;
    int localSize = nLoc * 2; // taille du système local (x et y)
    
    for (i = 0; i < nLoc; i++) {
        int myRowX = mapX[i];  // indice global pour la DOF x du nœud i
        int myRowY = mapY[i];  // indice global pour la DOF y du nœud i
        
        // Assembler les contributions du vecteur bloc (Bloc) pour les deux composantes
        mySolver->R[myRowX] -= Bloc[2 * i];     // contribution en x
        mySolver->R[myRowY] -= Bloc[2 * i + 1];   // contribution en y
        
        // Boucle sur les nœuds locaux pour assembler la matrice locale (Aloc)
        for (j = 0; j < nLoc; j++) {
            // Contribution de la ligne associée à la composante x du nœud i
            mySolver->R[myRowX] += 
                Aloc[(2 * i) * localSize + (2 * j)]     * Uloc[2 * j]     +
                Aloc[(2 * i) * localSize + (2 * j + 1)] * Uloc[2 * j + 1];
            
            // Contribution de la ligne associée à la composante y du nœud i
            mySolver->R[myRowY] += 
                Aloc[(2 * i + 1) * localSize + (2 * j)]     * Uloc[2 * j]     +
                Aloc[(2 * i + 1) * localSize + (2 * j + 1)] * Uloc[2 * j + 1];
        }
    }
}




double *femIterativeSolverEliminate(femIterativeSolver *mySolver)
{
    mySolver->iter++;
    double error = 0.0; int i;
    for (i=0; i < mySolver->size; i++) {
        error += (mySolver->R[i])*(mySolver->R[i]);
        mySolver->X[i] = -mySolver->R[i]/5.0; 
        mySolver->R[i] = 0.0; }
        
    mySolver->error = sqrt(error);
    return(mySolver->X);
}

 /* ======================================================================
    6. RENUMÉROTATION DU MAILLAGE
    ====================================================================== */
    #ifndef NORENUMBER
    double *positionMeshNodes;
    int comparPositionNode(const void *a, const void *b) {
        const int *nodePos_a = (const int *)a;
        const int *nodePos_b = (const int *)b;
        double diff = positionMeshNodes[*nodePos_a] - positionMeshNodes[*nodePos_b];
        return (diff < 0) - (diff > 0);
    }
    void femMeshRenumber(femMesh *theMesh, femRenumType renumType) {
        int nNodes = theMesh->nodes->nNodes;
        int *mapper = (int *)malloc(nNodes * sizeof(int));
        if (mapper == NULL) { Error("Memory allocation failed !"); exit(EXIT_FAILURE); return; }
        for (int i = 0; i < nNodes; i++) { mapper[i] = i; }
        switch (renumType) {
            case FEM_NO:
                break;
            case FEM_XNUM:
                positionMeshNodes = theMesh->nodes->X;
                qsort(mapper, nNodes, sizeof(int), comparPositionNode);
                break;
            case FEM_YNUM:
                positionMeshNodes = theMesh->nodes->Y;
                qsort(mapper, nNodes, sizeof(int), comparPositionNode);
                break;
            default:
                Error("Unexpected renumbering option");
        }
        for (int i = 0; i < nNodes; i++) {
            theMesh->nodes->number[mapper[i]] = i;
        }
        free(mapper);
    }
    #endif
    

 /* ======================================================================
    7. FONCTIONS DE GESTION DES ERREURS
    ====================================================================== */
 void femError(char *text, int line, char *file) {
     printf("\n-------------------------------------------------------------------------------- ");
     printf("\n  Error in %s at line %d : \n  %s\n", file, line, text);
     printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
     exit(69);
 }
 
 void femErrorGmsh(int ierr, int line, char *file) {
     if (ierr == 0) return;
     printf("\n-------------------------------------------------------------------------------- ");
     printf("\n  Error in %s at line %d : \n  error code returned by gmsh %d\n", file, line, ierr);
     printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
     gmshFinalize(NULL);
     exit(69);
 }
 
 void femErrorScan(int test, int line, char *file) {
     if (test >= 0) return;
     printf("\n-------------------------------------------------------------------------------- ");
     printf("\n  Error in fscanf or fgets in %s at line %d : \n", file, line);
     printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
     exit(69);
 }
 
 void femWarning(char *text, int line, char *file) {
     printf("\n-------------------------------------------------------------------------------- ");
     printf("\n  Warning in %s at line %d : \n  %s\n", file, line, text);
     printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
 }
 

 

 
#include <GL/gl.h> 
#include <math.h> 
#include "fem.h"

void glfemPlotFailureNodes(femNodes *nodes, double *sigmaXX, double *sigmaYY, double *sigmaXY, double sigma_max){
    glColor3f(1.0, 0.0, 0.0); 
    glLineWidth(5);
    // Taille de la croix (en unités d'écran)
    double crossSize = 0.002;

    glBegin(GL_LINES);
    for (int i = 0; i < nodes->nNodes; i++) {
        // Calcul de la contrainte équivalente (critère de von Mises)
        double sxx = sigmaXX[i];
        double syy = sigmaYY[i];
        double sxy = sigmaXY[i];
        double vonMises = sqrt(sxx * sxx - sxx * syy + syy * syy + 3.0 * sxy * sxy);

        // Si la contrainte dépasse la valeur critique, tracer une croix
        if (vonMises > sigma_max) {
            double x = nodes->X[i];
            double y = nodes->Y[i];

            // Tracer la première diagonale
            glVertex2d(x - crossSize, y - crossSize);
            glVertex2d(x + crossSize, y + crossSize);

            // Tracer la seconde diagonale
            glVertex2d(x - crossSize, y + crossSize);
            glVertex2d(x + crossSize, y - crossSize);
        }
    }
    glEnd();

}


 #endif
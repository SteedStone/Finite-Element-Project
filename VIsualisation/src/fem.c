/*
 *  fem.c
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2021 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#include "fem.h"

femGeo theGeometry;

femGeo *geoGetGeometry()                        { return &theGeometry; }

double geoSizeDefault(double x, double y)       { return 1.0; }

double geoGmshSize(int dim, int tag, double x, double y, double z, double lc, void *data)
                                                { return theGeometry.geoSize(x,y);    }
void geoInitialize() {
    int ierr;
    theGeometry.geoSize = geoSizeDefault;
    gmshInitialize(0,NULL,1,0,&ierr);                         ErrorGmsh(ierr);
    gmshModelAdd("MyGeometry",&ierr);                         ErrorGmsh(ierr);
    gmshModelMeshSetSizeCallback(geoGmshSize,NULL,&ierr);     ErrorGmsh(ierr);
    theGeometry.theNodes = NULL;
    theGeometry.theElements = NULL;
    theGeometry.theEdges = NULL;
    theGeometry.nDomains = 0;
    theGeometry.theDomains = NULL;
}

void geoFinalize() {
    int ierr;
    
    if (theGeometry.theNodes) {
        free(theGeometry.theNodes->X);
        free(theGeometry.theNodes->Y);
        free(theGeometry.theNodes); }
    if (theGeometry.theElements) {
        free(theGeometry.theElements->elem);
        free(theGeometry.theElements); }
    if (theGeometry.theEdges) {
        free(theGeometry.theEdges->elem);
        free(theGeometry.theEdges); }
    for (int i=0; i < theGeometry.nDomains; i++) {
        free(theGeometry.theDomains[i]->elem);
        free(theGeometry.theDomains[i]);  }
    free(theGeometry.theDomains);
    gmshFinalize(&ierr); ErrorGmsh(ierr);
}


void geoSetSizeCallback(double (*geoSize)(double x, double y)) {
    theGeometry.geoSize = geoSize; }


static int cmp_size_t(const void *a, const void *b){
    return (int) (*(const size_t*) a - * (const size_t*) b);
}

void geoMeshImport() {
    int ierr;
    
    /* Importing nodes */
    
    size_t nNode,n,m,*node;
    double *xyz,*trash;
 //   gmshModelMeshRenumberNodes(&ierr);                        ErrorGmsh(ierr);
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
    printf("Geo     : Importing %d elements \n",theElements->nElem);
    
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
   printf("Number of nodes %d \n", theNodes->nNodes);
   for (int i = 0; i < theNodes->nNodes; i++) {
      printf("%6d : %14.7e %14.7e \n",i,theNodes->X[i],theNodes->Y[i]); }
   femMesh *theEdges = theGeometry.theEdges;
   printf("Number of edges %d \n", theEdges->nElem);
   int *elem = theEdges->elem;
   for (int i = 0; i < theEdges->nElem; i++) {
      printf("%6d : %6d %6d \n",i,elem[2*i],elem[2*i+1]); }
   femMesh *theElements = theGeometry.theElements;
   printf("Number of triangles %d \n", theElements->nElem);
   elem = theElements->elem;
   for (int i = 0; i < theElements->nElem; i++) {
      printf("%6d : %6d %6d %6d\n",i,elem[3*i],elem[3*i+1],elem[3*i+2]); }
 
   int nDomains = theGeometry.nDomains;
   printf("Number of domains %d\n", nDomains);
   for (int iDomain = 0; iDomain < nDomains; iDomain++) {
      femDomain *theDomain = theGeometry.theDomains[iDomain];
      printf("  Domain : %6d \n", iDomain);
      printf("  Name : %s\n", theDomain->name);
      printf("  Number of elements : %6d\n", theDomain->nElem);
      for (int i=0; i < theDomain->nElem; i++){
 //         if (i != theDomain->nElem  && (i % 10) != 0)  printf(" - ");
          printf("%6d",theDomain->elem[i]);
          if ((i+1) != theDomain->nElem  && (i+1) % 10 == 0) printf("\n"); }
      printf("\n"); }
  
  
}


void geoMeshWrite(const char *filename) {
   FILE* file = fopen(filename,"w");
 
   femNodes *theNodes = theGeometry.theNodes;
   fprintf(file, "Number of nodes %d \n", theNodes->nNodes);
   for (int i = 0; i < theNodes->nNodes; i++) {
      fprintf(file,"%6d : %14.7e %14.7e \n",i,theNodes->X[i],theNodes->Y[i]); }
      
   femMesh *theEdges = theGeometry.theEdges;
   fprintf(file,"Number of edges %d \n", theEdges->nElem);
   int *elem = theEdges->elem;
   for (int i = 0; i < theEdges->nElem; i++) {
      fprintf(file,"%6d : %6d %6d \n",i,elem[2*i],elem[2*i+1]); }
      
   femMesh *theElements = theGeometry.theElements;
   fprintf(file,"Number of triangles %d \n", theElements->nElem);
   elem = theElements->elem;
   for (int i = 0; i < theElements->nElem; i++) {
      fprintf(file,"%6d : %6d %6d %6d\n",i,elem[3*i],elem[3*i+1],elem[3*i+2]); }
     
   int nDomains = theGeometry.nDomains;
   fprintf(file, "Number of domains %d\n", nDomains);
   for (int iDomain = 0; iDomain < nDomains; iDomain++) {
      femDomain *theDomain = theGeometry.theDomains[iDomain];
      fprintf(file, "  Domain : %6d \n", iDomain);
      fprintf(file, "  Name : %s\n", theDomain->name);
      fprintf(file, "  Number of elements : %6d\n", theDomain->nElem);
      for (int i=0; i < theDomain->nElem; i++){
          fprintf(file,"%6d",theDomain->elem[i]);
          if ((i+1) != theDomain->nElem  && (i+1) % 10 == 0) fprintf(file,"\n"); }
      fprintf(file,"\n"); }
    
   fclose(file);
}


void geoSetDomainName(int iDomain, char *name) {
  if (iDomain >= theGeometry.nDomains)  Error("Illegal domain number");
  sprintf(theGeometry.theDomains[iDomain]->name,"%s",name);
} 

void geoMeshRead(const char *filename) 
{
   FILE* file = fopen(filename,"r");
   
   int trash, *elem;
   
   femNodes *theNodes = malloc(sizeof(femNodes));
   theGeometry.theNodes = theNodes;
   ErrorScan(fscanf(file, "Number of nodes %d \n", &theNodes->nNodes));
   theNodes->X = malloc(sizeof(double)*(theNodes->nNodes));
   theNodes->Y = malloc(sizeof(double)*(theNodes->nNodes));
   for (int i = 0; i < theNodes->nNodes; i++) {
       ErrorScan(fscanf(file,"%d : %le %le \n",&trash,&theNodes->X[i],&theNodes->Y[i]));} 

   femMesh *theEdges = malloc(sizeof(femMesh));
   theGeometry.theEdges = theEdges;
   theEdges->nLocalNode = 2;
   theEdges->nodes = theNodes;
   ErrorScan(fscanf(file, "Number of edges %d \n", &theEdges->nElem));
   theEdges->elem = malloc(sizeof(int)*theEdges->nLocalNode*theEdges->nElem);
   for(int i=0; i < theEdges->nElem; ++i) {
        elem = theEdges->elem;
        ErrorScan(fscanf(file, "%6d : %6d %6d \n", &trash,&elem[2*i],&elem[2*i+1])); }
  
   femMesh *theElements = malloc(sizeof(femMesh));
   theGeometry.theElements = theElements;
   theElements->nLocalNode = 0;
   theElements->nodes = theNodes;
   char elementType[MAXNAME];  
   ErrorScan(fscanf(file, "Number of %s %d \n",elementType,&theElements->nElem));  
   if (strncasecmp(elementType,"triangles",MAXNAME) == 0) {
      theElements->nLocalNode = 3;
      theElements->elem = malloc(sizeof(int)*theElements->nLocalNode*theElements->nElem);
      for(int i=0; i < theElements->nElem; ++i) {
          elem = theElements->elem;
          ErrorScan(fscanf(file, "%6d : %6d %6d %6d \n", 
                    &trash,&elem[3*i],&elem[3*i+1],&elem[3*i+2])); }}
   if (strncasecmp(elementType,"quads",MAXNAME) == 0) {
      theElements->nLocalNode = 4;
      theElements->elem = malloc(sizeof(int)*theElements->nLocalNode*theElements->nElem);
      for(int i=0; i < theElements->nElem; ++i) {
          elem = theElements->elem;
          ErrorScan(fscanf(file, "%6d : %6d %6d %6d %6d \n", 
                    &trash,&elem[4*i],&elem[4*i+1],&elem[4*i+2],&elem[4*i+3])); }}
           
   ErrorScan(fscanf(file, "Number of domains %d\n", &theGeometry.nDomains));
   int nDomains = theGeometry.nDomains;
   theGeometry.theDomains = malloc(sizeof(femDomain*)*nDomains);
   for (int iDomain = 0; iDomain < nDomains; iDomain++) {
      femDomain *theDomain = malloc(sizeof(femDomain)); 
      theGeometry.theDomains[iDomain] = theDomain;
      theDomain->mesh = theEdges; 
      ErrorScan(fscanf(file,"  Domain : %6d \n", &trash));
      ErrorScan(fscanf(file,"  Name : %[^\n]s \n", (char*)&theDomain->name));
      ErrorScan(fscanf(file,"  Number of elements : %6d\n", &theDomain->nElem));
      theDomain->elem = malloc(sizeof(int)*2*theDomain->nElem); 
      for (int i=0; i < theDomain->nElem; i++){
          ErrorScan(fscanf(file,"%6d",&theDomain->elem[i]));
          if ((i+1) != theDomain->nElem  && (i+1) % 10 == 0) ErrorScan(fscanf(file,"\n")); }}
    
   fclose(file);
}

int geoGetDomain(char *name)
{
    int theIndex = -1;
    int nDomains = theGeometry.nDomains;
    for (int iDomain = 0; iDomain < nDomains; iDomain++) {
        femDomain *theDomain = theGeometry.theDomains[iDomain];
        if (strncasecmp(name,theDomain->name,MAXNAME) == 0)
            theIndex = iDomain;  }
    return theIndex;
            
}

static const double _gaussQuad4Xsi[4]    = {-0.577350269189626,-0.577350269189626, 0.577350269189626, 0.577350269189626};
static const double _gaussQuad4Eta[4]    = { 0.577350269189626,-0.577350269189626,-0.577350269189626, 0.577350269189626};
static const double _gaussQuad4Weight[4] = { 1.000000000000000, 1.000000000000000, 1.000000000000000, 1.000000000000000};
static const double _gaussTri3Xsi[3]     = { 0.166666666666667, 0.666666666666667, 0.166666666666667};
static const double _gaussTri3Eta[3]     = { 0.166666666666667, 0.166666666666667, 0.666666666666667};
static const double _gaussTri3Weight[3]  = { 0.166666666666667, 0.166666666666667, 0.166666666666667};
static const double _gaussEdge2Xsi[2]    = { 0.577350269189626,-0.577350269189626};
static const double _gaussEdge2Weight[2] = { 1.000000000000000, 1.000000000000000};



femIntegration *femIntegrationCreate(int n, femElementType type)
{
    femIntegration *theRule = malloc(sizeof(femIntegration));
    if (type == FEM_QUAD && n == 4) {
        theRule->n      = 4;
        theRule->xsi    = _gaussQuad4Xsi;
        theRule->eta    = _gaussQuad4Eta;
        theRule->weight = _gaussQuad4Weight; }
    else if (type == FEM_TRIANGLE && n == 3) {
        theRule->n      = 3;
        theRule->xsi    = _gaussTri3Xsi;
        theRule->eta    = _gaussTri3Eta;
        theRule->weight = _gaussTri3Weight; }
    else if (type == FEM_EDGE && n == 2) {
        theRule->n      = 2;
        theRule->xsi    = _gaussEdge2Xsi;
        theRule->eta    = NULL;
        theRule->weight = _gaussEdge2Weight; }
    else Error("Cannot create such an integration rule !");
    return theRule; 
}

void femIntegrationFree(femIntegration *theRule)
{
    free(theRule);
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



femDiscrete *femDiscreteCreate(int n, femElementType type)
{
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
        theSpace->n       = 4;
        theSpace->x2      = _q1c0_x;
        theSpace->phi2    = _q1c0_phi;
        theSpace->dphi2dx = _q1c0_dphidx; }
    else if (type == FEM_TRIANGLE && n == 3) {
        theSpace->n       = 3;
        theSpace->x2      = _p1c0_x;
        theSpace->phi2    = _p1c0_phi;
        theSpace->dphi2dx = _p1c0_dphidx; }
    else if (type == FEM_EDGE && n == 2) {
        theSpace->n       = 2;
        theSpace->x       = _e1c0_x;
        theSpace->phi     = _e1c0_phi;
        theSpace->dphidx  = _e1c0_dphidx; }
    else Error("Cannot create such a discrete space !");
    return theSpace; 
}

void femDiscreteFree(femDiscrete *theSpace)
{
    free(theSpace);
}

void femDiscreteXsi2(femDiscrete* mySpace, double *xsi, double *eta)
{
    mySpace->x2(xsi,eta);
}

void femDiscretePhi2(femDiscrete* mySpace, double xsi, double eta, double *phi)
{
    mySpace->phi2(xsi,eta,phi);
}

void femDiscreteDphi2(femDiscrete* mySpace, double xsi, double eta, double *dphidxsi, double *dphideta)
{
    mySpace->dphi2dx(xsi,eta,dphidxsi,dphideta);
}

void femDiscreteXsi(femDiscrete* mySpace, double *xsi)
{
    mySpace->x(xsi);
}

void femDiscretePhi(femDiscrete* mySpace, double xsi, double *phi)
{
    mySpace->phi(xsi,phi);
}

void femDiscreteDphi(femDiscrete* mySpace, double xsi, double *dphidxsi)
{
    mySpace->dphidx(xsi,dphidxsi);
}

void femDiscretePrint(femDiscrete *mySpace)
{
    int i,j;
    int n = mySpace->n;
    double xsi[4], eta[4], phi[4], dphidxsi[4], dphideta[4];

    if (mySpace->type == FEM_EDGE) {
        femDiscreteXsi(mySpace,xsi);
        for (i=0; i < n; i++) {           
            femDiscretePhi(mySpace,xsi[i],phi);
            femDiscreteDphi(mySpace,xsi[i],dphidxsi);
            for (j=0; j < n; j++)  {
                printf("(xsi=%+.1f) : ",xsi[i]);
                printf(" phi(%d)=%+.1f",j,phi[j]);  
                printf("   dphidxsi(%d)=%+.1f \n",j,dphidxsi[j]); }
            printf(" \n"); }}
    
    if (mySpace->type == FEM_QUAD || mySpace->type == FEM_TRIANGLE) {
        femDiscreteXsi2(mySpace, xsi, eta);
        for (i = 0; i < n; i++)  {    
            femDiscretePhi2(mySpace, xsi[i], eta[i], phi);
            femDiscreteDphi2(mySpace, xsi[i], eta[i], dphidxsi, dphideta);
            for (j = 0; j < n; j++) {  
                printf("(xsi=%+.1f,eta=%+.1f) : ", xsi[i], eta[i]);  
                printf(" phi(%d)=%+.1f", j, phi[j]);
                printf("   dphidxsi(%d)=%+.1f", j, dphidxsi[j]);
                printf("   dphideta(%d)=%+.1f \n", j, dphideta[j]); }
            printf(" \n"); }}   
}

femFullSystem *femFullSystemCreate(int size)
{
    femFullSystem *theSystem = malloc(sizeof(femFullSystem));
    femFullSystemAlloc(theSystem, size);
    femFullSystemInit(theSystem);

    return theSystem; 
}

void femFullSystemFree(femFullSystem *theSystem)
{
    free(theSystem->A);
    free(theSystem->B);
    free(theSystem);
}

void femFullSystemAlloc(femFullSystem *mySystem, int size)
{
    int i;  
    double *elem = malloc(sizeof(double) * size * (size+1)); 
    mySystem->A = malloc(sizeof(double*) * size); 
    mySystem->B = elem;
    mySystem->A[0] = elem + size;  
    mySystem->size = size;
    for (i=1 ; i < size ; i++) 
        mySystem->A[i] = mySystem->A[i-1] + size;
}

void femFullSystemInit(femFullSystem *mySystem)
{
    int i,size = mySystem->size;
    for (i=0 ; i < size*(size+1) ; i++) 
        mySystem->B[i] = 0;}


void femFullSystemPrint(femFullSystem *mySystem)
{
    double  **A, *B;
    int     i, j, size;
    
    A    = mySystem->A;
    B    = mySystem->B;
    size = mySystem->size;
    
    for (i=0; i < size; i++) {
        for (j=0; j < size; j++)
            if (A[i][j] == 0)  printf("         ");   
            else               printf(" %+.1e",A[i][j]);
        printf(" :  %+.1e \n",B[i]); }
}

double* femFullSystemEliminate(femFullSystem *mySystem)
{
    double  **A, *B, factor;
    int     i, j, k, size;
    
    A    = mySystem->A;
    B    = mySystem->B;
    size = mySystem->size;
    
    /* Gauss elimination */
    
    for (k=0; k < size; k++) {
        if ( fabs(A[k][k]) <= 1e-16 ) {
            printf("Pivot index %d  ",k);
            printf("Pivot value %e  ",A[k][k]);
            Error("Cannot eliminate with such a pivot"); }
        for (i = k+1 ; i <  size; i++) {
            factor = A[i][k] / A[k][k];
            for (j = k+1 ; j < size; j++) 
                A[i][j] = A[i][j] - A[k][j] * factor;
            B[i] = B[i] - B[k] * factor; }}
    
    /* Back-substitution */
    
    for (i = size-1; i >= 0 ; i--) {
        factor = 0;
        for (j = i+1 ; j < size; j++)
            factor += A[i][j] * B[j];
        B[i] = ( B[i] - factor)/A[i][i]; }
    
    return(mySystem->B);    
}

void  femFullSystemConstrain(femFullSystem *mySystem, 
                             int myNode, double myValue) 
{
    double  **A, *B;
    int     i, size;
    
    A    = mySystem->A;
    B    = mySystem->B;
    size = mySystem->size;
    
    for (i=0; i < size; i++) {
        B[i] -= myValue * A[i][myNode];
        A[i][myNode] = 0; }
    
    for (i=0; i < size; i++) 
        A[myNode][i] = 0; 
    
    A[myNode][myNode] = 1;
    B[myNode] = myValue;
}


femProblem *femElasticityCreate(femGeo* theGeometry, 
                  double E, double nu, double rho, double g, femElasticCase iCase)
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
    if (theGeometry->theElements->nLocalNode == 3) {
        theProblem->space    = femDiscreteCreate(3,FEM_TRIANGLE);
        theProblem->rule     = femIntegrationCreate(3,FEM_TRIANGLE); }
    if (theGeometry->theElements->nLocalNode == 4) {
        theProblem->space    = femDiscreteCreate(4,FEM_QUAD);
        theProblem->rule     = femIntegrationCreate(4,FEM_QUAD); }
    theProblem->spaceEdge    = femDiscreteCreate(2,FEM_EDGE);
    theProblem->ruleEdge     = femIntegrationCreate(2,FEM_EDGE); 
    theProblem->system       = femFullSystemCreate(size); 

    femDiscretePrint(theProblem->space);   
    femDiscretePrint(theProblem->spaceEdge);  
  
    return theProblem;
}

void femElasticityFree(femProblem *theProblem)
{
    femFullSystemFree(theProblem->system);
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
            int node = theBoundary->domain->mesh->elem[2*elem[e]+i];
            theProblem->constrainedNodes[2*node+shift] = size-1; }}    
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

void femError(char *text, int line, char *file)                                  
{ 
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in %s at line %d : \n  %s\n", file, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
    exit(69);                                                 
}

void femErrorGmsh(int ierr, int line, char *file)                                  
{ 
    if (ierr == 0)  return;
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in %s at line %d : \n  error code returned by gmsh %d\n", file, line, ierr);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
    gmshFinalize(NULL);                                        
    exit(69);                                                 
}

void femErrorScan(int test, int line, char *file)                                  
{ 
    if (test >= 0)  return;
    
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in fscanf or fgets in %s at line %d : \n", file, line);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");   
    exit(69);                                       
}

void femWarning(char *text, int line, char *file)                                  
{ 
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Warning in %s at line %d : \n  %s\n", file, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");                                              
}


void trianglePlot() {

    int ierr;

    // gmshModelAdd("TriangleModel",  &ierr); // Créer un modèle pour les triangles
    // ErrorGmsh(ierr);


    femGeo* theGeometry = geoGetGeometry();
    double hexRadius = theGeometry->hexRadius;
    int numHexX = theGeometry->NumberOfTrianglesInX;
    int numHexY = theGeometry->NumberOfTrianglesInY;

    // Définition des tailles de maillage variables
    double meshSizeMin = 0.07;  // Taille de maille en haut (fine)
    double meshSizeMax = 0.3;   // Taille de maille en bas (grossière)
    double distance = hexRadius * 1.1; // Distance entre les lignes

    int innerWireTags[1000];
    int innerWireCount = 0;

    for (int i = 0; i < numHexX; i++) {
        for (int j = 0; j < numHexY; j++) {
            double x = i * 1.3 * hexRadius;
            double y = j * distance - hexRadius;

            // Calcul dynamique du meshSize en fonction de la hauteur
            double meshSize = meshSizeMax - (meshSizeMax - meshSizeMin) * ((double)j / numHexY);

            int innerPoints[3];

            if (i % 2 == 0) {
                innerPoints[0] = gmshModelOccAddPoint(x, y, 0, meshSize, -1, &ierr);
                innerPoints[1] = gmshModelOccAddPoint(x + hexRadius, y + sqrt(3) * hexRadius / 2, 0, meshSize, -1, &ierr);
                innerPoints[2] = gmshModelOccAddPoint(x - hexRadius, y + sqrt(3) * hexRadius / 2, 0, meshSize, -1, &ierr);
            } else {
                innerPoints[0] = gmshModelOccAddPoint(x, y + sqrt(3) * hexRadius / 2, 0, meshSize, -1, &ierr);
                innerPoints[1] = gmshModelOccAddPoint(x + hexRadius, y, 0, meshSize, -1, &ierr);
                innerPoints[2] = gmshModelOccAddPoint(x - hexRadius, y, 0, meshSize, -1, &ierr);
            }

            ErrorGmsh(ierr);

            int innerLines[3];
            for (int k = 0; k < 3; k++) {
                innerLines[k] = gmshModelOccAddLine(innerPoints[k], innerPoints[(k + 1) % 3], -1, &ierr);
                ErrorGmsh(ierr);
            }

            int innerWire = gmshModelOccAddWire(innerLines, 3, -1, 1, &ierr);
            ErrorGmsh(ierr);

            innerWireTags[innerWireCount++] = innerWire;
        }
    }

    int rectPoints[4];
    rectPoints[0] = gmshModelOccAddPoint(0, 0, 0, meshSizeMax, -1, &ierr);
    rectPoints[1] = gmshModelOccAddPoint(numHexX * 1.3 * hexRadius - (1.5) * hexRadius, 0, 0, meshSizeMax, -1, &ierr);
    rectPoints[2] = gmshModelOccAddPoint(numHexX * 1.3 * hexRadius - (1.5) * hexRadius, numHexY * distance - hexRadius, 0, meshSizeMin, -1, &ierr);
    rectPoints[3] = gmshModelOccAddPoint(0, numHexY * distance - hexRadius, 0, meshSizeMin, -1, &ierr);

    int rectLines[4];
    for (int k = 0; k < 4; k++) {
        rectLines[k] = gmshModelOccAddLine(rectPoints[k], rectPoints[(k + 1) % 4], -1, &ierr);
        ErrorGmsh(ierr);
    }

    int rectWire = gmshModelOccAddWire(rectLines, 4, -1, 1, &ierr);
    ErrorGmsh(ierr);

    int plateSurface = gmshModelOccAddPlaneSurface(&rectWire, 1, -1, &ierr);
    ErrorGmsh(ierr);

    int triSurfaces[100];
    for (int i = 0; i < innerWireCount; i++) {
        triSurfaces[i] = gmshModelOccAddPlaneSurface(&innerWireTags[i], 1, -1, &ierr);
        ErrorGmsh(ierr);
    }

    int notch[innerWireCount][2];
    for (int i = 0; i < innerWireCount; i++) {
        notch[i][0] = 2;
        notch[i][1] = triSurfaces[i];
    }

    int plate[2] = {2, plateSurface};
    gmshModelOccCut(plate, 2, (int *)notch, 2 * innerWireCount, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    ErrorGmsh(ierr);

    gmshModelOccSynchronize(&ierr);
    ErrorGmsh(ierr);


    

    gmshModelMeshGenerate(2, &ierr);
    ErrorGmsh(ierr);

    // gmshFltkRun(&ierr);
    // ErrorGmsh(ierr);
}



void HexagonPlot(){

//
//  -1- Construction de la g�om�trie avec OpenCascade
//      On cr�e le rectangle
//      On cr�e les deux cercles
//      On soustrait les cercles du rectangle :-)
//

    int ierr;

    // gmshModelAdd("HexagonModel",  &ierr);
    // ErrorGmsh(ierr);
    
    femGeo* theGeometry = geoGetGeometry();

    double hexRadius = theGeometry->hexRadius;  // Rayon de l'hexagone
    int numHexX = theGeometry->NumberOfHexagonsInX;  // Nombre d'hexagones en largeur
    int numHexY = theGeometry->NumberOfHexagonsInY;  // Nombre d'hexagones en hauteur
    double meshSizeMin = 0.07;  // Taille de maille en haut (fine)
    double meshSizeMax = 0.3;   // Taille de maille en bas (grossière)
    double f = 0.81;

    int mainWireTags[1000]; // Contient tous les hexagones extérieurs
    int innerWireTags[1000]; // Contient les petits hexagones à soustraire
    int trianglecoteTags[1000]; // Contient les triangles à soustraire
    int quadrilatereTags[1000]; // Contient les quadrilateres à soustraire
    int losangeTags[1000]; // Contient les losanges à soustraire
    

    int wireCount = 0, innerWireCount = 0, triangleWireCount = 0, quadrilatereWireCount = 0, losangeWireCount = 0;

    for (int i = 0; i < numHexX; i++) {
        for (int j = 0; j < numHexY; j++) {
            double x = i * 1.5 * hexRadius;
            double y = j * sqrt(3) * hexRadius + (i % 2) * sqrt(3) * hexRadius / 2;
            int trianglecotePoints[3], trianglecoteLines[3], previoustrianglepoints[numHexY];
            double meshSize = meshSizeMax - (meshSizeMax - meshSizeMin) * ((double)j / numHexY);
            if (j == 0 && i == 0) {
                for (int k = 0; k < 2; k++)
                {
                    double angle = M_PI / 3 * (k+3);
                    double px = x + hexRadius * cos(angle);
                    double py = y + hexRadius * sin(angle);
                    trianglecotePoints[k] = gmshModelOccAddPoint(px, py, 0, meshSize, -1, &ierr);
                    ErrorGmsh(ierr);
                }
                trianglecotePoints[2] = gmshModelOccAddPoint(-hexRadius, -hexRadius*sqrt(3)/2, 0, meshSize, -1, &ierr);
                previoustrianglepoints[0] = gmshModelOccAddPoint(hexRadius * cos(M_PI), hexRadius*sin(M_PI), 0, meshSize, -1, &ierr);
                ErrorGmsh(ierr);
                for (int k = 0; k < 3; k++)
                {
                    trianglecoteLines[k] = gmshModelOccAddLine(trianglecotePoints[k], trianglecotePoints[(k+1)%3], -1, &ierr);
                    ErrorGmsh(ierr);
                }
                
            }
            if(i == 0 && j > 0 ) {
                for (int k = 0; k < 2; k++)
                {
                    double angle = M_PI / 3 * (k+3);
                    double px = x + hexRadius * cos(angle);
                    double py = y + hexRadius * sin(angle);
                    trianglecotePoints[k] = gmshModelOccAddPoint(px, py, 0, meshSize, -1, &ierr);
                    ErrorGmsh(ierr);
                }
                trianglecotePoints[2] = previoustrianglepoints[j-1];
                previoustrianglepoints[j] = gmshModelOccAddPoint(hexRadius * cos(M_PI) + x, hexRadius*sin(M_PI)+ y, 0, meshSize, -1, &ierr);
                ErrorGmsh(ierr);
                for (int k = 0; k < 3; k++)
                {
                    trianglecoteLines[k] = gmshModelOccAddLine(trianglecotePoints[k], trianglecotePoints[(k+1)%3], -1, &ierr);
                    ErrorGmsh(ierr);
                }
            }
            if(i == numHexX-1 && j == 0) {
                for (int k = 0; k < 2; k++)
                {
                    double angle = -M_PI / 3 * (k);
                    double px = x + hexRadius * cos(angle);
                    double py = y + hexRadius * sin(angle);
                    trianglecotePoints[k] = gmshModelOccAddPoint(px, py, 0, meshSize, -1, &ierr);
                    ErrorGmsh(ierr);
                }
                trianglecotePoints[2] = gmshModelOccAddPoint((numHexX -1)* 1.5* hexRadius + hexRadius, -hexRadius*sqrt(3)/2, 0, meshSize, -1, &ierr);
                previoustrianglepoints[0] = gmshModelOccAddPoint(hexRadius * cos(0)+ x, hexRadius*sin(0), 0, meshSize, -1, &ierr);
                ErrorGmsh(ierr);
                for (int k = 0; k < 3; k++)
                {
                    trianglecoteLines[k] = gmshModelOccAddLine(trianglecotePoints[k], trianglecotePoints[(k+1)%3], -1, &ierr);
                    ErrorGmsh(ierr);
                }
            }
            if(i == numHexX-1 && j > 0 ) {
                for (int k = 0; k < 2; k++)
                {
                    double angle = -M_PI / 3 * (k);
                    double px = x + hexRadius * cos(angle);
                    double py = y + hexRadius * sin(angle);
                    trianglecotePoints[k] = gmshModelOccAddPoint(px, py, 0, meshSize, -1, &ierr);
                    ErrorGmsh(ierr);
                }
                trianglecotePoints[2] = previoustrianglepoints[j-1];
                previoustrianglepoints[j] = gmshModelOccAddPoint(hexRadius * cos(0) + x, hexRadius*sin(0)+ y, 0, meshSize, -1, &ierr);
                ErrorGmsh(ierr);
                for (int k = 0; k < 3; k++)
                {
                    trianglecoteLines[k] = gmshModelOccAddLine(trianglecotePoints[k], trianglecotePoints[(k+1)%3], -1, &ierr);
                    ErrorGmsh(ierr);
                }
            }
            int quadrilaterePoints[5], quadrilatereLines[5];
            if( i == 0 && j == numHexY - 1) {
                for (int k = 0; k < 3; k++)
                {
                    double angle = M_PI / 3 * (k+1);
                    double px = x + hexRadius * cos(angle);
                    double py = y + hexRadius * sin(angle);
                    quadrilaterePoints[k] = gmshModelOccAddPoint(px, py, 0, meshSize, -1, &ierr);
                    ErrorGmsh(ierr);
                }
                quadrilaterePoints[3] = gmshModelOccAddPoint(-hexRadius, numHexY * sqrt(3) * hexRadius, 0, meshSize, -1, &ierr);
                quadrilaterePoints[4] = gmshModelOccAddPoint(1 * 1.5 * hexRadius + hexRadius * cos(M_PI/3 * 2), j * sqrt(3) * hexRadius + (1 % 2) * sqrt(3) * hexRadius / 2 + hexRadius * sin(M_PI/3 * 2), 0, meshSize, -1, &ierr);
                ErrorGmsh(ierr);
                for (int k = 0; k < 5; k++)
                {
                    quadrilatereLines[k] = gmshModelOccAddLine(quadrilaterePoints[k], quadrilaterePoints[(k+1)%5], -1, &ierr);
                    ErrorGmsh(ierr);
                }
            }
            if( i == numHexX - 1 && j == numHexY - 1 ){
                for (int k = 0; k < 3; k++)
                {
                    double angle = M_PI / 3 * (k);
                    double px = x + hexRadius * cos(angle);
                    double py = y + hexRadius * sin(angle);
                    quadrilaterePoints[k] = gmshModelOccAddPoint(px, py, 0, meshSize, -1, &ierr);
                    ErrorGmsh(ierr);
                }
                quadrilaterePoints[4] = gmshModelOccAddPoint((numHexX -1 ) * 1.5 * hexRadius + hexRadius, numHexY * sqrt(3) * hexRadius, 0, meshSize, -1, &ierr);
                quadrilaterePoints[3] = gmshModelOccAddPoint((numHexX - 2) * 1.5 * hexRadius + hexRadius * cos(M_PI/3 ), j * sqrt(3) * hexRadius + ((numHexX - 2) % 2) * sqrt(3) * hexRadius / 2 + hexRadius * sin(M_PI/3 ), 0, meshSize, -1, &ierr);
                ErrorGmsh(ierr);
                for (int k = 0; k < 5; k++)
                {
                    quadrilatereLines[k] = gmshModelOccAddLine(quadrilaterePoints[k], quadrilaterePoints[(k+1)%5], -1, &ierr);
                    ErrorGmsh(ierr);
                }
            }
            int losangePoints[4], losangeLines[4];
            if(i % 2 != 0 && j == 0) {
                for (int k = 1; k < 3; k++) {
                    double angle = -M_PI / 3 * (k);
                    double px = x + hexRadius * cos(angle);
                    double py = y + hexRadius * sin(angle);
                    losangePoints[k] = gmshModelOccAddPoint(px, py, 0, meshSize, -1, &ierr);
                    ErrorGmsh(ierr);
                }
                losangePoints[3] = gmshModelOccAddPoint((i - 1) * 1.5 * hexRadius + hexRadius * cos(-M_PI/3 ), j * sqrt(3) * hexRadius + ((i- 1) % 2) * sqrt(3) * hexRadius / 2 + hexRadius * sin(-M_PI/3 ), 0, meshSize, -1, &ierr);
                losangePoints[0] = gmshModelOccAddPoint((i + 1) * 1.5 * hexRadius + hexRadius * cos(-M_PI/3 *2), j * sqrt(3) * hexRadius + ((i+ 1) % 2) * sqrt(3) * hexRadius / 2 + hexRadius * sin(-M_PI/3 *2), 0, meshSize, -1, &ierr);
                for (int k = 0; k < 4; k++) {
                    losangeLines[k] = gmshModelOccAddLine(losangePoints[k], losangePoints[(k + 1) % 4], -1, &ierr);
                    ErrorGmsh(ierr);
                }
                int losangeWire = gmshModelOccAddWire(losangeLines, 4, -1, 1, &ierr);
                ErrorGmsh(ierr);
                losangeTags[losangeWireCount++] = losangeWire;
            }
            if( i % 2 == 0 && j == numHexY - 1 && i != 0 && i != numHexX - 1) {
                for (int k = 1; k < 3; k++) {
                    double angle = M_PI / 3 * (k);
                    double px = x + hexRadius * cos(angle);
                    double py = y + hexRadius * sin(angle);
                    losangePoints[k] = gmshModelOccAddPoint(px, py, 0, meshSize, -1, &ierr);
                    ErrorGmsh(ierr);
                }
                losangePoints[3] = gmshModelOccAddPoint((i - 1) * 1.5 * hexRadius + hexRadius * cos(M_PI/3 ), j * sqrt(3) * hexRadius + ((i- 1) % 2) * sqrt(3) * hexRadius / 2 + hexRadius * sin(M_PI/3 ), 0, meshSize, -1, &ierr);
                losangePoints[0] = gmshModelOccAddPoint((i + 1) * 1.5 * hexRadius + hexRadius * cos(M_PI/3 *2), j * sqrt(3) * hexRadius + ((i+ 1) % 2) * sqrt(3) * hexRadius / 2 + hexRadius * sin(M_PI/3 *2), 0, meshSize, -1, &ierr);
                for (int k = 0; k < 4; k++) {
                    losangeLines[k] = gmshModelOccAddLine(losangePoints[k], losangePoints[(k + 1) % 4], -1, &ierr);
                    ErrorGmsh(ierr);
                }
                int losangeWire = gmshModelOccAddWire(losangeLines, 4, -1, 1, &ierr);
                ErrorGmsh(ierr);
                losangeTags[losangeWireCount++] = losangeWire;
            }
        
           
            // triangle interne 
            int points[6], innerPoints[6];
            for (int k = 0; k < 6; k++) {
                double angle = M_PI / 3 * k;
                double px = x + hexRadius * cos(angle);
                double py = y + hexRadius * sin(angle);
                points[k] = gmshModelOccAddPoint(px, py, 0, meshSize, -1, &ierr);
                ErrorGmsh(ierr);

                double innerPx = x + (hexRadius * f) * cos(angle);
                double innerPy = y + (hexRadius * f) * sin(angle);
                innerPoints[k] = gmshModelOccAddPoint(innerPx, innerPy, 0, meshSize, -1, &ierr);
                ErrorGmsh(ierr);
            }

            int lines[6], innerLines[6];
            for (int k = 0; k < 6; k++) {
                lines[k] = gmshModelOccAddLine(points[k], points[(k + 1) % 6], -1, &ierr);
                ErrorGmsh(ierr);
                innerLines[k] = gmshModelOccAddLine(innerPoints[k], innerPoints[(k + 1) % 6], -1, &ierr);
                ErrorGmsh(ierr);
            }

            int wire = gmshModelOccAddWire(lines, 6, -1, 1, &ierr);
            ErrorGmsh(ierr);
            int innerWire = gmshModelOccAddWire(innerLines, 6, -1, 1, &ierr);
            ErrorGmsh(ierr);

            mainWireTags[wireCount++] = wire;
            innerWireTags[innerWireCount++] = innerWire;
            if (j == 0 && i==0) {
                int trianglecoteWire = gmshModelOccAddWire(trianglecoteLines, 3, -1, 1, &ierr);
                ErrorGmsh(ierr);
                trianglecoteTags[triangleWireCount++] = trianglecoteWire;
            }
            if(i == 0 && j > 0 ) {
                int trianglecoteWire = gmshModelOccAddWire(trianglecoteLines, 3, -1, 1, &ierr);
                ErrorGmsh(ierr);
                trianglecoteTags[triangleWireCount++] = trianglecoteWire;
            }
            if(i == numHexX-1 && j == 0) {
                int trianglecoteWire = gmshModelOccAddWire(trianglecoteLines, 3, -1, 1, &ierr);
                ErrorGmsh(ierr);
                trianglecoteTags[triangleWireCount++] = trianglecoteWire;
            }
            if(i == numHexX-1 && j > 0) {
                int trianglecoteWire = gmshModelOccAddWire(trianglecoteLines, 3, -1, 1, &ierr);
                ErrorGmsh(ierr);
                trianglecoteTags[triangleWireCount++] = trianglecoteWire;
            }
            if( i == 0 && j == numHexY - 1) {
                int quadrilatereWire = gmshModelOccAddWire(quadrilatereLines, 5, -1, 1, &ierr);
                ErrorGmsh(ierr);
                quadrilatereTags[0] = quadrilatereWire;
            }
            if( i == numHexX - 1 && j == numHexY - 1) {
                int quadrilatereWire = gmshModelOccAddWire(quadrilatereLines, 5, -1, 1, &ierr);
                ErrorGmsh(ierr);
                quadrilatereTags[1] = quadrilatereWire;
            }
            
            
        }
    }

    // Création d'un grand rectangle englobant la structure
    int rectPoints[4];
    rectPoints[0] = gmshModelOccAddPoint(-hexRadius, -hexRadius*sqrt(3)/2 - hexRadius/3, 0, meshSizeMax, -1, &ierr);
    rectPoints[1] = gmshModelOccAddPoint((numHexX -1)* 1.5* hexRadius + hexRadius, -hexRadius*sqrt(3)/2 - hexRadius/3, 0, meshSizeMax, -1, &ierr);
    rectPoints[2] = gmshModelOccAddPoint((numHexX -1)* 1.5* hexRadius + hexRadius, numHexY * sqrt(3) * hexRadius+ hexRadius/3, 0, meshSizeMin, -1, &ierr);
    rectPoints[3] = gmshModelOccAddPoint(-hexRadius, numHexY * sqrt(3) * hexRadius+ hexRadius/3, 0, meshSizeMin, -1, &ierr);

    // int newRectPoints[4];
    // newRectPoints[0] = gmshModelOccAddPoint(-hexRadius, -hexRadius*sqrt(3)/2 , 0, meshSizeMax, -1, &ierr);
    // newRectPoints[1] = gmshModelOccAddPoint((numHexX -1)* 1.5* hexRadius + hexRadius, -hexRadius*sqrt(3)/2 , 0, meshSizeMax, -1, &ierr);
    // newRectPoints[2] = gmshModelOccAddPoint((numHexX -1 ) * 1.5 * hexRadius + hexRadius, numHexY * sqrt(3) * hexRadius, 0, meshSizeMin, -1, &ierr);
    // newRectPoints[3] = gmshModelOccAddPoint(-hexRadius, numHexY * sqrt(3) * hexRadius, 0, meshSizeMin, -1, &ierr);
    

    int rectLines[4];
    for (int k = 0; k < 4; k++) {
        rectLines[k] = gmshModelOccAddLine(rectPoints[k], rectPoints[(k + 1) % 4], -1, &ierr);
        ErrorGmsh(ierr);
    }

    // int newRectLines[4];
    // for (int k = 0; k < 4; k++) {
    //     newRectLines[k] = gmshModelOccAddLine(newRectPoints[k], newRectPoints[(k + 1) % 4], -1, &ierr);
    //     ErrorGmsh(ierr);
    // }


    int rectWire = gmshModelOccAddWire(rectLines, 4, -1, 1, &ierr);
    ErrorGmsh(ierr);

    

    // Création de la surface principale
    int plateSurface = gmshModelOccAddPlaneSurface(&rectWire, 1, -1, &ierr);
    ErrorGmsh(ierr);

    // Création des surfaces hexagonales
    int hexSurfaces[100];
    for (int i = 0; i < innerWireCount; i++) {
        hexSurfaces[i] = gmshModelOccAddPlaneSurface(&innerWireTags[i], 1, -1, &ierr);
        ErrorGmsh(ierr);
    }

    // Création des surfaces triangulaires
    int triangleSurfaces[100];
    for (int i = 0; i < triangleWireCount; i++) {
        triangleSurfaces[i] = gmshModelOccAddPlaneSurface(&trianglecoteTags[i], 1, -1, &ierr);
        ErrorGmsh(ierr);


    }
    // Creation des surfaces quadrilateres
    int quadrilatereSurfaces[2];
    for (int i = 0; i < 2; i++) {
        quadrilatereSurfaces[i] = gmshModelOccAddPlaneSurface(&quadrilatereTags[i], 1, -1, &ierr);
        ErrorGmsh(ierr);
    }
    int plate[2] = {2, plateSurface};
    // // Découper les quadrilateres
    int notchquadrilatere[2][2];
    for (int i = 0; i < 2; i++) {
        notchquadrilatere[i][0] = 2 ;
        notchquadrilatere[i][1] = quadrilatereSurfaces[i] ;
    }
    gmshModelOccCut(plate, 2, (int *)notchquadrilatere, 2*2  ,NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccSynchronize(&ierr);

    // Creation des surfaces losanges
    int losangeSurfaces[100];
    for (int i = 0; i < losangeWireCount; i++) {
        losangeSurfaces[i] = gmshModelOccAddPlaneSurface(&losangeTags[i], 1, -1, &ierr);
        ErrorGmsh(ierr);
    }
    // Découper les losanges
    int notchlosange[losangeWireCount][2];
    for (int i = 0; i < losangeWireCount; i++) {
        notchlosange[i][0] = 2 ;
        notchlosange[i][1] = losangeSurfaces[i] ;
    }
    gmshModelOccCut(plate, 2, (int *)notchlosange, 2*losangeWireCount  ,NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccSynchronize(&ierr);

    // Découper les triangles 
    int notchtriangle[triangleWireCount][2];
    for (int i = 0; i < triangleWireCount; i++) {
        notchtriangle[i][0] = 2 ;
        notchtriangle[i][1] = triangleSurfaces[i] ;
    }
    gmshModelOccCut(plate, 2, (int *)notchtriangle, 2*triangleWireCount  ,NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccSynchronize(&ierr);

    // Découper les hexagones de la plaque
    int notch[innerWireCount][2];
    for (int i = 0; i < innerWireCount; i++) {
        notch[i][0] = 2;
        notch[i][1] = hexSurfaces[i];
    }

    gmshModelOccCut(plate, 2, (int *)notch, 2*innerWireCount, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    ErrorGmsh(ierr);

    // Synchronisation apres ajout des surfaces
    gmshModelOccSynchronize(&ierr);
    ErrorGmsh(ierr);

    // Génération du maillage
    gmshModelMeshGenerate(2, &ierr);
    ErrorGmsh(ierr);

    // Affichage dans Gmsh
    // gmshFltkRun(&ierr);
    // ErrorGmsh(ierr);



}

void femFindBoundaryNodes(femGeo *theProblem, double targetY, double epsilon)
{
    
    femMesh *theEdges = theProblem->theEdges;

    int size = theProblem->theElements->nElem;
    printf("Size: %d\n", size);
    int *map = malloc(size * sizeof(int));
    if (map == NULL) { Error("Memory allocation failed !"); exit(EXIT_FAILURE); return; }
    for (int i = 0; i < size; i++) { map[i] = 0; }

    int nDomains = theProblem->nDomains;
    for (int iDomain = 0; iDomain < nDomains; iDomain++)
    {
        femDomain *theDomain = theProblem->theDomains[iDomain];
        for (int i = 0; i < theDomain->nElem; i++)
        {
           int iEdge = theDomain->elem[i];
            for (int j = 0; j < 2; j++)
            {
               int iNode = theEdges->elem[iEdge * 2 + j];
               map[iNode] = 1;
            }
        }
    }

    // Compter les nœuds de la frontière qui sont à targetY
    int nFilteredBoundary = 0;
    for (int i = 0; i < size; i++)
    {
        printf("Node %d: %f\n", i, theProblem->theNodes->Y[i]);
        if (map[i] == 1 && fabs(theProblem->theNodes->Y[i] - targetY) < epsilon) { 
            nFilteredBoundary++; 
        }
    }

    // Création du domaine pour les nœuds filtrés
    femDomain *theFilteredBoundary = malloc(sizeof(femDomain));
    if (theFilteredBoundary == NULL) { Error("Memory allocation failed !"); exit(EXIT_FAILURE); return; }
    theProblem->nDomains++;
    theProblem->theDomains = realloc(theProblem->theDomains, theProblem->nDomains * sizeof(femDomain *));
    if (theProblem->theDomains == NULL) { Error("Memory allocation failed !"); exit(EXIT_FAILURE); return; }
    theProblem->theDomains[theProblem->nDomains - 1] = theFilteredBoundary;
    theFilteredBoundary->nElem = nFilteredBoundary;
    theFilteredBoundary->elem = malloc(nFilteredBoundary * sizeof(int));
    if (theFilteredBoundary->elem == NULL) { Error("Memory allocation failed !"); exit(EXIT_FAILURE); return; }

    theFilteredBoundary->mesh = NULL;
    sprintf(theFilteredBoundary->name, "FilteredBoundary");

    // Stocker les nœuds filtrés
    int index = 0;
    for (int i = 0; i < size; i++)
    {
        if (map[i] == 1 && fabs(theProblem->theNodes->Y[i] - targetY) < epsilon) {
            theFilteredBoundary->elem[index++] = i; 
            printf("Filtered node: %d\n", i);
        }
    }
    
    free(map);
}

void geoMeshGenerate() {


    femGeo* theGeometry = geoGetGeometry();
    if(theGeometry->hexa_triangles == 1) {
        trianglePlot();
    } else {
        HexagonPlot();
    }
   

    geoSetSizeCallback(geoSize);

    int ierr;
    gmshModelOccSynchronize(&ierr);
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);
    if (ierr != 0) {
        printf("Erreur lors de la génération du maillage: %d\n", ierr);
        exit(1);
    }
}
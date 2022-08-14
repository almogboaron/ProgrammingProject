#define PY_SSIZE_T_CLEAN
#include <Python/Python.h>

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <float.h>


/*Declorations*/
int K;
char* filename;
int max_iter = 300;
int n;
int d;
double EPSILON = 0;
int iter_num = 0;
int convergenceCount;
double delta;
int* count;

/*Centroids Decloration*/
double* Centroids;
double** CenetroidAloc;

/*Clusters Can be represented Only By The Sum of The Data Points*/
double* ClusterSum;
double** ClusterSumAloc;

/*Data Points decloration*/
double* data;
double** dataAloc;

/*Weighted Adjacency Matrix decloration*/
double* wam;
double** wamAloc;

/*The Diagonal Degree Matrix*/
double* ddg;
double** ddgAloc;

/*The Normalized Graph Laplacian*/
double* lnorm;
double** lnormAloc;

/*Python Objects*/
PyObject* data_arr;
PyObject* centroid_arr;

/*Vectors matrix*/
double* vectors;
double** vectorsAloc;

/*define assert with message*/
#define assert__(x) for (;!(x);assert(x))


/* -------------------------------- Handling Data ------------------- */

/* Function Cols Dims. */
int numOfCols(){
    int countcols=0;
    char ch;
    FILE *fp = fopen(filename,"r");
    assert__(fp!=NULL){
        printf("An Error Has Occurred");
    }

    while(!feof(fp)){
        ch = fgetc(fp);
        if(ch == ','){
            countcols++;
        }
        if(ch=='\n'){break;}
    }

    countcols++;
    fclose(fp);
    return countcols;

}

/*Counting Rows Dims.*/
int numbOfRows(){
    char ch ;
    int countrows=0;
    FILE *fp = fopen(filename,"r");
    assert__(fp!=NULL){
        printf("An Error Has Occurred");
    }

    while(!feof(fp)){
        ch = fgetc(fp);
        if(ch == '\n'){
            countrows++;
        }
    }
    fclose(fp);
    return countrows;
}

/*reads file and initialize data matrix*/
void init_DataMat(){
    FILE* fp;
    int i,j;

    /*Initialization*/
    data = calloc(n*d,sizeof(double));
    assert__(data!=NULL){
        printf("An Error Has Occurred");
    }
    dataAloc = calloc(n,sizeof(double));
    assert__(dataAloc!=NULL){
        printf("An Error Has Occurred");
    }

    /*Initializing Pointer to Data Mat (k*d) (Or Vectors X1,...Xn == R1....Rn)*/
    for( i=0; i<n ; i++){
        dataAloc[i] = data + i*d ;
    }

    /*Initialize Data Matrix( Data Array n*d);*/
    fp = fopen(filename,"r");
    assert__(fp!=NULL){
        printf("An Error Has Occurred");
    }

    for(i=0;i<n;i++){
        for(j=0;j<d;j++){
            fscanf(fp,"%lf",&dataAloc[i][j]);
            fgetc(fp);
        }
    }
}

void read_file(){
    n = numbOfRows();
    assert__(K<=n){
        printf("Invalid Input!");
    }
    d = numOfCols();
    init_DataMat();
}

/* ---------------------------------- k-means functions ----------------------------------*/

/*Function Initalizeation Centroids -> Matrix Formation as first k data points */
void Init_Centroids(){
    int i;
    int j;

    /*Initializing Centroids Mat*/
    for(i=0;i<K;i++){
        for(j=0;j<d;j++){
            Centroids[i*d+j] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(centroid_arr,i),j));
        }
    }
    /*Initializing Pointer to Centroids Mat (k*d)*/
    for(i=0; i < K ; i++){
        CenetroidAloc[i] = Centroids + i*d ;
    }
}

/*Initializeation of ClusterSum*/
void init_Clusters(){
    int i;

    /*Initializing Centroids Mat*/
    for(i=0 ; i<K*d ;i++){
        ClusterSum[i] = 0;
    }
    /*Initializing Pointer to Centroids Mat (k*d)*/
    for(i=0; i < K ; i++){
        ClusterSumAloc[i] = ClusterSum + i*d ;
    }
}

/*Sum Data Point with appropriate Cluster*/
void SumWithCluster(double* clusterSum , double* dataPoint){
    int i;
    for(i=0; i<d; i++){
        clusterSum[i] = clusterSum[i] + dataPoint[i];
    }
}

/* Sum DataPoint with Appropriate Cluster Sum return idx of cluster*/
int assign(double *datapoint){
    double minDist = DBL_MAX;
    double norm=0 ;
    int i,j;
    int idx=0;

    /*Go through all Centroids and choose the closest to dataPoint*/
    for(i=0; i<K; i++){
        norm = 0;

        /*Calculating Delta */
        for(j=0; j<d; j++){
            norm += pow(datapoint[j] - CenetroidAloc[i][j],2);

        }

        /*Checking MinDistance*/
        if (norm < minDist){
            minDist = norm;
            idx = i;
        }
    }
    /*Summing the DataPoint With the closest Cluster.*/
    SumWithCluster(ClusterSumAloc[idx], datapoint);
    return idx;
}

/*Devide Cluster with number of points*/
void NormelizeClusterSums(int* counter){
    int i,j;
    for(i=0; i<K; i++){
        for(j=0;j<d;j++){
            ClusterSumAloc[i][j] = ClusterSumAloc[i][j]/counter[i];
        }
    }
}

/*Update Centroinds, Resets Cluster  And returns the EuclidianDistance.*/
double update_centroid(double* Centroid , double* ClusterSum){
    double norm=0;
    int i;
    for(i=0; i<d; i++){
        norm +=pow(Centroid[i] - ClusterSum[i],2);
        Centroid[i] = ClusterSum[i];
        ClusterSum[i] = 0;
    }
    return sqrt(norm);
}

/*-------------------Wam-------------------*/
/*Norm Implementation*/
double norm2(int i, int j){
    int k;
    double res = 0;

    if(i==j){
        return 0;
    }

    for(k=0;k<n;k++){
        res+=pow(dataAloc[i][k]-dataAloc[j][k],2);
    }

    return sqrt(res);
}

/*Wam Initialization and Implementation*/
void wamInit(){
    int i,j;
    /*Initialize*/
    wam = calloc(n*n,sizeof(double));
    assert__(wam!=NULL){
        printf("An Error Has Occurred");
    }
    wamAloc = calloc(n,sizeof(double*));
    assert__(wamAloc!=NULL){
        printf("An Error Has Occurred");
    }
    
    for(i=0;i<n;i++){
        wamAloc[i] = wam + i*n;
    }
    /*Implementation*/

    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            wamAloc[i][j] = exp(-(norm2(i,j)/2));
        }
    }
}

/*------------------Diagonal Degree Matrix--------------*/
void ddgInit(){
    int i,z;
    /*Initialize*/
    ddg = calloc(n*n,sizeof(double));
    assert__(ddg!=NULL){
        printf("An Error Has Occurred");
    }
    ddgAloc = calloc(n,sizeof(double*));
    assert__(ddgAloc!=NULL){
        printf("An Error Has Occurred");
    }

    for(i=0;i<n;i++){
        ddgAloc[i] = ddg + i*n;
    }

    /*Implementation*/
    for(i=0;i<n;i++){
        for(z=0;z<n;z++){
            ddgAloc[i][i] += wamAloc[i][z];
        }
    }
}

/*---------------The Normalized Graph Laplacian---------*/
/*Matrix Multiplacation*/
void MatMulti(double** mat1,double** mat2 , double** matRes){
    int i,j,k;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            for(k=0;k<n;k++){
                matRes[i][j]+=mat1[i][k]*mat2[k][j];
            }
        }
    }
}

/*Matrix Substraction*/
void MatSubEye( double** matRes){
    int i,j;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(i==j){ matRes[i][j] = 1 - matRes[i][j];}
            else{matRes[i][j]= 0 - matRes[i][j];}
        }
    }
}

/*Lnorm Implementation*/
void lnormInit(){
    int i;
    double* res;
    double** resAloc;
    /*Initialize*/
    lnorm = calloc(n*n,sizeof(double));
    assert__(lnorm!=NULL){
        printf("An Error Has Occurred");
    }
    lnormAloc = calloc(n,sizeof(double*));
    assert__(lnormAloc!=NULL){
        printf("An Error Has Occurred");
    }
    res = calloc(n*n,sizeof(double));
    assert__(lnorm!=NULL){
        printf("An Error Has Occurred");
    }
    resAloc = calloc(n,sizeof(double*));
    assert__(lnormAloc!=NULL){
        printf("An Error Has Occurred");
    }
    for(i=0;i<n;i++){
        lnormAloc[i] = lnorm + i*n;
        resAloc[i] = res +i*n;
    }

    /*Implementation*/
    for(i=0;i<n;i++){
        ddgAloc[i][i] = 1/sqrt(ddgAloc[i][i]);
    }
    MatMulti(wamAloc,ddgAloc,resAloc);
    MatMulti(ddgAloc,resAloc,lnormAloc);
    MatSubEye(lnormAloc);
    free(res);free(resAloc);
}

/*---------------Jacobi Algorithem---------*/

/*---------------------------PYTHON C API-----------------------------*/
/*Create PyCentroids matrix*/
PyObject* cMatToPy(double** mat, int rows){
    int i;
    int j;
    PyObject* pyMat = PyList_New(0);
    for(i=0;i<rows;i++){
        PyObject* PyPoint = PyList_New(0);
        for(j=0;j<d;j++){
            if(PyList_Append(PyPoint,PyFloat_FromDouble(mat[i][j]))==-1){
                return NULL;
            }
        }
        if(PyList_Append(pyMat,PyPoint)==-1){
            return NULL;
        }
    }
    return pyMat;
}
/*Main Function of Kmeans Algorithem .*/
static PyObject* kmeans() {
    int i;
    int idx;

    /*Initialization*/
    ClusterSum = calloc(K*d,sizeof(double));
    assert__(ClusterSum!=NULL){
        printf("An Error Has Occurred");
    }
    ClusterSumAloc = calloc(K,sizeof(double));
    assert__(ClusterSumAloc!=NULL){
        printf("An Error Has Occurred");
    }
    Centroids = calloc(K*d,sizeof(double));
    assert__(Centroids!=NULL){
        printf("An Error Has Occurred");
    }
    CenetroidAloc = calloc(K,sizeof(double));
    assert__(CenetroidAloc!=NULL){
        printf("An Error Has Occurred");
    }

    Init_Centroids();
    init_Clusters();

    count = calloc(K,sizeof(int));
    assert__(count!=NULL){
        printf("An Error Has Occurred");
    }
    while (iter_num < max_iter){

        /*Reset:countPerCluster , ConvergenceCount;*/
        convergenceCount = 0;
        memset(count, 0, K*sizeof(int));

        /*Assigning Data to clustersSum.*/
        for(i=0; i < n; i++){
            idx = assign(dataAloc[i]);
            count[idx] += 1;
        }
        /*Normalize Cluster*/
        NormelizeClusterSums(count);

        /*Updates Centroid, Reset ClusterSum, Returns Delta.*/
        for(i=0; i<K; i++){
            delta = update_centroid(CenetroidAloc[i] ,ClusterSumAloc[i]);
            if (delta < EPSILON){
                convergenceCount+=1;
            }
        }

        /*Convergence or Continue Iteration;*/
        if (convergenceCount == K){break;}
        iter_num++;
    }
/*Update pyCentroids And return that object*/
PyObject* pyCentroids = cMatToPy(CenetroidAloc,K);

/*Free Memory*/
free(ClusterSum);
free(ClusterSumAloc);
free(Centroids);
free(CenetroidAloc);
free(data);
free(dataAloc);
free(count);
return pyCentroids;
}

static PyObject* spk() {
    /*getting the data from the file*/
    read_file();
    /*Asserts*/
    assert__(K<=n){
        printf("Invalid Input!");
    }
    wamInit(); /*step 1*/
    lnormInit(); /*step 2*/
    /*step 3*/
    if (K == 0){
        /*K = Eigngap_Heuristic(); need to be implemented*/
    }
    /*step 4*/
    /*jacobi(); need to be implemented*/

    /*step 5 - renormalize U's rows - NEED TO ADD*/

    /*step 6 - passing the data back to python*/
    PyObject* U = cMatToPy(vectorsAloc,n);
    return U;
}

void print_output(double** mat){
    int i,j;
    for(i=0; i<n ;i++){
        for(j=0; j<n; j++){
            /*Not Last Column*/
            if(j<n-1){printf("%.4f,",mat[i][j]);}
            else    {printf("%.4f",mat[i][j]);}
        }
        printf("%s","\n");
    }
}

int main(int argc, char *argv[]) {
    assert__(argc == 3){
        printf("Invalid Input!");
    }
    filename = argv[2];
    read_file();

    if(strcmp(argv[1],"wam") == 0){
        wamInit();
        print_output(wamAloc);
    }
    else if(strcmp(argv[1],"ddg") == 0){
        ddgInit();
        print_output(ddgAloc);
    }
    else if(strcmp(argv[1],"lnorm") == 0){
        lnormInit();
        print_output(lnormAloc);
    }

    free(data);
    free(dataAloc);
    free(wam);
    free(wamAloc);
    free(ddg);
    free(ddgAloc);
    free(lnorm);
    free(lnormAloc);
    return 0;
}




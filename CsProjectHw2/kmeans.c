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
char* filename_in;
char* filename_out;
int max_iter;
int n;
int d;
double EPSILON = 0.001;
int iter_num = 0;
int convergenceCount;
int countcols=0;
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

/*Python Objects*/
PyObject* data_arr; 
PyObject* centroid_arr; 

/*define assert with message*/
#define assert__(x) for (;!(x);assert(x))

/*Function Initalizeation Centroids -> Matrix Formation as first k data points */
void Init_Centroids(){
    int i;
    int j;

    /*Initializing Centroids Mat*/
    for(i=0;i<K;i++){
        for(j=0;j<d;j++){
            Centroids[i*d+j] = PyFloat_AsDouble(PYList_getItem(PyList_GetItem(centroid_arr,i),j))
        }
    }
    /*Initializing Pointer to Centroids Mat (k*d)*/
    for(i=0; i < K ; i++){
        CenetroidAloc[i] = Centroids + i*d ;
    }
}

/*Function Read PythonData From List of Lists (Matrix) and initialize C Matrix (n*d) of Information*/ 
void init_DataMat(){
    int i;
    int j;

    for(i=0;i<n;i++){
        for(j=0;j<d;j++){
            data[i*d+j] = PyFloat_AsDouble(PYList_getItem(PyList_GetItem(data_arr,i),j))
        }
    }
    /*Initializing Pointer to Data Mat (k*d) (Or Vectors X1,...Xn == R1....Rn)*/
    for( i=0; i<n ; i++){
        dataAloc[i] = data + i*d ; 
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

/*Sum Data Point with apropriate Cluster*/
void SumWithCluster(double* clusterSum , double* dataPoint){
    int i;
    for(i=0; i<d; i++){
        clusterSum[i] = clusterSum[i] + dataPoint[i];
    }
}

/* Sum DataPoint with Apropritate Cluster Sum return idx of cluster*/
int assign(double *datapoint){
    double minDist = DBL_MAX;
    double norm=0 ;
    int i,j;
    int idx;

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

/*Create PyCentroids matrix*/
Pyobject* cCentroidsToPy(){
    int i;
    int j;
    PyObject* pyCentroids = PyList_New(0);
    for(i=0;i<K;i++){
        PyObject* PyCentroid = PyList_New(0);
        for(j=0;j<d;j++){
            if(PyList_Append(PyCentroid,PyFloat_FromDouble(CenetroidAloc[i][j]))==-1){
                return NULL;
            }
        }
        if(PyList_Append(pyCentroids,pyCentroid)==-1){
            return NULL;
        }
    }
    return pyCentroids;
}
/*Main Function of Kmeans Algorithem .*/
static PyObject* kmeans(int argc, char *argv[]) {
    int i;
    int idx;
    /*Asserts*/
    assert__(K<=n){
        printf("Invalid Input!");
    }
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
    data = calloc(n*d,sizeof(double));
    assert__(data!=NULL){
        printf("An Error Has Occurred");
    }
    dataAloc = calloc(n,sizeof(double));
    assert__(dataAloc!=NULL){
        printf("An Error Has Occurred");
    }

    init_DataMat();
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
PyObject* pyCentroids = cCentroidsToPy();
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

/*the func that operates when we call the mykmeanssp module from python,
 * the fuction will get the cmd arguments that we proccesed in python and call the main function
 * we wrote in hw1 with some diffrences.
 * I renamed the main func to kmeans.
 */

static PyObject* kmeans_capi(PyObject *self, PyObject *args){
    /* passing arguments from py- K,num of rows, num of cols, max_iter, epsilon,
     * data and initial centroinds*/

    if (!PyArg_ParseTuple(args,"iiiidO!O!", &K, &n, &d, &max_iter, &EPSILON, &data_arr, &centroid_arr)){
        return NULL;
    }
    /* don't forget to change the epsilon init*/

    return kmeans();
    /* returning double or array*/
}

/* declaring the kmeans function */
static PyMethodDef capiMethods[] = {
        {"fit", (PyCFunction) kmeans_capi, METH_VARARGS, PyDoc_STR("kmeans algorithem")},
                 {NULL,NULL,0,NULL}
};

/*defining the module*/
static struct PyModuleDef capimodule = {
        PyModuleDef_HEAD_INIT, "mykmeanssp", NULL, -1, capiMethods
};

/* creating the module */
PyMODINIT_FUNC
PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if(!m) {
        return NULL;
    }
    return m;
}

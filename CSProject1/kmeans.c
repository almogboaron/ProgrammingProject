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

/*Function Initalizeation Centroids -> Matrix Formation as first k data points */
void Init_Centroids(){
    int i;

    /*Initializing Centroids Mat*/
    for(i=0 ; i<K*d ;i++){
        Centroids[i] = data[i];
    }
    /*Initializing Pointer to Centroids Mat (k*d)*/
    for(i=0; i < K ; i++){
        CenetroidAloc[i] = Centroids + i*d ;
    }
}

/*Function Read Files and initialize Matrix (n*d) of Information*/ 
void init_DataMat(){
    FILE* fp;
    int i,j;

    /*Initializing Pointer to Data Mat (k*d) (Or Vectors X1,...Xn == R1....Rn)*/
    for( i=0; i<n ; i++){
        dataAloc[i] = data + i*d ; 
    }

    /*Initialize Data Matrix( Data Array n*d);*/
    fp = fopen(filename_in,"r");
    assert(fp!=NULL);

    for(i=0;i<n;i++){
        for(j=0;j<d;j++){
            fscanf(fp,"%lf",&dataAloc[i][j]);
            fgetc(fp);
        }
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

/* Function Cols Dims. */
int numOfCols(){
    int countcols=0;
    char ch;
    FILE *fp = fopen(filename_in,"r");
    assert(fp!=NULL);
    
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
    FILE *fp = fopen(filename_in,"r");
    assert(fp!=NULL);
    
    while(!feof(fp)){
        ch = fgetc(fp);
        if(ch == '\n'){
            countrows++;
        }
    }
    fclose(fp);
    return countrows;
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
    double norm;
    int i;
    for(i=0; i<d; i++){
        norm +=pow(Centroid[i] - ClusterSum[i],2);
        Centroid[i] = ClusterSum[i];
        ClusterSum[i] = 0;
    }
    return sqrt(norm);
}

void WriteBackCentroids(){
    int i,j;
    FILE *fp = fopen(filename_out,"w");
    assert(fp!=NULL);

    for(i=0; i<K ;i++){
        for(j=0; j<d; j++){
                /*Not Last Cloumn*/
                if(j<d-1){fprintf(fp,"%.4f,",CenetroidAloc[i][j]);}
                else    {fprintf(fp,"%.4f",CenetroidAloc[i][j]);}
        }
        fprintf(fp,"%s","\n");
    }
    fclose(fp);
}

/*Main Function of Kmeans Algorithem .*/
int main(int argc, char *argv[]) {
    int i;
    int idx;
    /*Make Assertions!*/

    /*Without MaxIter*/
    if (argc == 4) { 
        K = atoi(argv[1]);
        max_iter = 200;
        filename_in = argv[2];
        filename_out = argv[3];
    }
    /*With MaxIter*/
    else {
        K = atoi(argv[1]);
        max_iter = atoi(argv[2]);
        filename_in = argv[3];
        filename_out = argv[4]; 
    }

    /*Initialization*/
    n = numbOfRows();
    d = numOfCols();

    ClusterSum = calloc(K*d,sizeof(double));
    ClusterSumAloc = calloc(K,sizeof(double));
    Centroids = calloc(K*d,sizeof(double));
    CenetroidAloc = calloc(K,sizeof(double));
    data = calloc(n*d,sizeof(double));
    dataAloc = calloc(n,sizeof(double));

    init_DataMat();
    Init_Centroids();
    init_Clusters();
    
    count = calloc(K,sizeof(int));
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

WriteBackCentroids();
free(ClusterSum);
free(ClusterSumAloc);
free(Centroids);
free(CenetroidAloc);
free(data);
free(dataAloc);
free(count);
return 1;
}




#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

//Declorations
int K;
char* filename_in;
char* filename_out;
int max_iter;
int n;
int d;
double EPSILON = 0.001;

//Centroids Decloration
double* Centroids;
double** CenetroidAloc;

//Clusters Can be represented Only By The Sum of The Data Points
double* ClusterSum;
double** ClusterSumAloc;

//Data Points decloration
double* data;
double** dataAloc;

//Function Initalizeation Centroids -> Matrix Formation as first k data points 
void Init_Centroids(){
    Centroids = calloc(K*d,sizeof(double));
    CenetroidAloc = calloc(K,sizeof(int));
    //Initializing Centroids Mat
    int i;
    for(i=0 ; i<K*d ;i++){
        Centroids[i] = data[i];
    }
    //Initializing Pointer to Centroids Mat (k*d)
    for(i=0; i < K ; i++){
        CenetroidAloc[i] = Centroids + i*d ;
    }
}

//Function Read Files and initialize Matrix (n*d) of Information 
void init_DataMat(){
    double a;
    int i=0;
    int j=0;
    char line , number ;
    
    //Call location
    data = calloc(n*d,sizeof(double));
    dataAloc = calloc(n,sizeof(int));
    
    //Initializing Pointer to Data Mat (k*d) (Or Vectors X1,...Xn == R1....Rn)
    for( i=0; i<n ; i++){
        dataAloc[i] = data + i*d ; 
    }

    //Initialize Data Matrix( Data Array n*d); 
    FILE *fp = fopen(filename_in,"r");
    assert(fp!=NULL);
    while(!feof(fp)){
        fscanf(fp,"%s",line);
        number = *strtok(&line,",");
        while(number!=NULL){
            data[i++] = (double)atof(&number);
            number = *strtok(NULL,",");
        }
    }
}

//Initializeation of ClusterSum
void init_Clusters(){
    ClusterSum = calloc(K*d,sizeof(double));
    ClusterSumAloc = calloc(K,sizeof(int));
    //Initializing Centroids Mat
    int i;
    for(i=0 ; i<K*d ;i++){
        ClusterSum[i] = 0;
    }
    //Initializing Pointer to Centroids Mat (k*d)
    for(i=0; i < K ; i++){
        ClusterSumAloc[i] = ClusterSum + i*d ;
    }
}

// Function Cols Dims. 
int numOfCols(){
    int countcols=0;
    char *token,*line;
    FILE *fp = fopen(filename_in,"r");
    assert(fp!=NULL);
    fscanf(fp,"%s",line);
    fclose(fp);
    token = strtok(line,",");
    while (token != NULL){
        countcols++;
        token = strtok(NULL,",");
    }
    return countcols;

}

//Counting Rows Dims.
int numbOfRows(){
    char *line;
    int countrows=0;
    FILE *fp = fopen(filename_in,"r");
    assert(fp!=NULL);
    while(feof(fp)){
        fscanf(fp,"%s",line);
        ++countrows;
    }
    fclose(fp);
    return countrows;
}

// Sum DataPoint with Apropritate Cluster Sum return idx of cluster
int assign(double *datapoint){
    double minDist = 1.7976931348623158E+308;
    double norm ;
    int idx;
    int j;
    int i;

    //Go through all Centroids and choose the closest to dataPoint
    for(j=0; j<K; j++){
        
        //Calculating Delta 
        for(i=0; i<d; i++){
            norm += datapoint[i] - CenetroidAloc[j][i];
        }

        //Calculation Norm
        norm = sqrt(norm);
        
        //Checking MinDistance
        if (norm < minDist){
            minDist = norm;
            idx = j;
        }
    }

    //Summing the DataPoint With the closest Cluster.
   SumWithCluster(ClusterSumAloc[idx], datapoint);
   return idx;
}

//Sum Data Point with apropriate Cluster
void SumWithCluster(double* clusterSum , double* dataPoint){
    int i;
    for(i=0; i<d; i++){
        clusterSum[i] = clusterSum[i] + dataPoint[i];
    }
}

//Devide Cluster with number of points
void NormelizeClusterSums(int* counter){
    int i,j;
    for(i=0; i<K; i++){
        for(j=0;j<d;j++){
            ClusterSumAloc[i][j] = ClusterSumAloc[i][j]/counter[i];
        }
    }
}

//Update Centroinds, Resets Cluster  And returns the EuclidianDistance.
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
    FILE *fp = fopen(filename_out,"w");
    assert(fp!=NULL);
    int i,j;
    for(i=0; i<K ;i++){
        for(j=0; j<d; j++){
                if(j==d-1){fprintf(fp,"%d,",CenetroidAloc[i][j]);}
                else    {fprintf(fp,"%d",CenetroidAloc[i][j]);}
        }
        fprintf(fp,"%s","\n");
    }
    fclose(fp);
}

//Main Function of Kmeans Algorithem .
int main(int argc, char *argv[]) {
    
    //Without MaxIter
    if (argc == 2) { 
        K = atoi(argv[0]);
        max_iter = 200;
        filename_in = argv[1];
        filename_out = argv[2];
    }
    
    //With MaxIter
    else { 
        K = atoi(argv[0]);
        max_iter = atoi(argv[1]);
        filename_in = argv[2];
        filename_out = argv[3]; 
    }

    //Initialization
    n = numbOfRows();
    d = numOfCols();

    init_DataMat();
    Init_Centroids();
    init_Clusters();
    
    int i;
    int idx;
    int iter_num = 0;
    int count[K];
    int convergenceCount;
    while (iter_num < max_iter){
        //Reset:countPerCluster , ConvergenceCount;
        convergenceCount = 0;
        memset(count, 0, sizeof count);

        //Assigning Data to clustersSum.
        for(i=0; i < n; i++){
            idx = assign(dataAloc[i]);
            count[idx] += 1;
        }
        
        //Normalize Cluster
        NormelizeClusterSums(count);

        //Updates Centroid, Reset ClusterSum, Returns Delta.
        double delta;
        for(i=0; i<K; i++){
            delta = update_centroid(CenetroidAloc[i] ,ClusterSumAloc[i]);
            if (delta < EPSILON){
                convergenceCount+=1;
            }
        }
        //Convergence or Continue Iteration;
        if (convergenceCount == K){break;}
        iter_num++;
    }
WriteBackCentroids();  
}




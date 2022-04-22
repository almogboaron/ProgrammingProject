#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

//DataPoint type(Pointer or index in matrixData)
typedef int DATA; 

//Struct for ClusterPed list = Node = Cluster.
struct linked_list{
    DATA d;
    struct linked_list *next;
};

//Calling a linked_list A cluster 
typedef struct linked_list Cluster;

//Pointer of Cluster is a ClusterP.
typedef Cluster* ClusterP;


//Don't think we need ( just adding to a new clusters each time while once cheking distance between centroid to datapoint)
//Iterative List Creation (Lecture 4 slide 6)
//Did not finish right (did not finish slide)
ClusterP centroid_to_list(int DataIndex){
    ClusterP head=NULL, tail=NULL;
    head = (Cluster*)malloc(sizeof(Cluster));
    assert(head != NULL);
    head-> d=centroid;
    return head;

}
// Adding A datapoint to A list 
void add_to_list(ClusterP list, int datapoint){
    ClusterP tail = list->next;
    tail-> next = (ELEMENT*) malloc(sizeof(ELEMENT));
    tail = tail->next;
    assert(tail!=NULL);
    tail-> d=datapoint;
    tail-> next = NULL;
}

//Declorations
int K;
char* filename;
int max_iter;
int n;
int d;
int *Centroids;
int **CenetroidAloc;
int *data;
int **dataAloc;
ClusterP *clusters;
double EPSILON = 0.001;


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
    FILE *fp = fopen(filename,"r");
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

// Function Cols Dims. 
int numOfCols(){
    int countcols=0;
    char *token,*line;
    FILE *fp = fopen(filename,"r");
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
    FILE *fp = fopen(filename,"r");
    assert(fp!=NULL);
    while(feof(fp)){
        fscanf(fp,"%s",line);
        ++countrows;
    }
    fclose(fp);
    return countrows;
}

//Euclidian norm calc *No sqrt needed*
double norm_calc(double *delta){
    double s = 0;
    int i;
    for(i=0; i<d; i++){
        s += pow(delta[i],2);
    }
    //Changed to no Sqrt
    return s;
}

//Assign DataPoint pointer to Cluster
void assign(int *datapoint){
    double val = 1.7976931348623158E+308;
    int idx;
    int j;
    //Go through all Centroids and choose the closest to dataPoint
    for(j=0; j<K; j++){
        double delta[d];
        int i;
        //Calculating Delta 
        for(i=0; i<d; i++){
            delta[i] = datapoint[i] - CenetroidAloc[j][i];
        }
        //Calculation Norm
        double norm = norm_calc(delta);
        if (norm < val){
            val = norm;
            idx = j;
        }
    }
    //Adding The DataPoint to A cluster.
    add_to_list(clusters[idx],*datapoint);
}

//Update Centroinds and returns the Change.
double update_centroid(){
    int i;
    int j;
    // Allocating Array for calculation of data in clusters.
    double dataSum[d]
    
    for (i=0;i<k;i++){
         calc_cluster(i,dataSum)
        for(j=0; j<d; j++){
            CentroidsAloc[i][j] = datasum[j]
        }
    }
    //Freeing Datasum Array
    free(dataSum)
}

//Calculating sum of all datapoints and then deviding in number of datapoints in cluster.
void calc_cluster(int i, double* dataSum){

}

//Initialize Clusters as K first DataPoints.
void init_clusters(){
    clusters = calloc(K,sizeof(ClusterP));
    int i;
    for(i=0; i<K; i++){
        ClusterP list_i = centroid_to_list(*CenetroidAloc[i]);
        clusters[i] = list_i;
    }
}

//Main Function of Kmeans Algorithem .
int main(int argc, char *argv[]) {
    if (argc == 2) { max_iter = 200; }
    else { max_iter = atoi(argv[2]); }
    K = atoi(argv[0]);
    filename = argv[1];
    n = numbOfRows();
    d = numOfCols();
    init_DataMat();
    Init_Centroids();
    init_clusters();
    
    int i;
    int iter_num = 0;
    while (iter_num < max_iter){
        //Assigning Data to clusters.
        for(i=0; i < n; i++){
            assign(data[i]);
        }
        
        //For every Cetroid update from cluster and calculate delta.
        int cnt = 0;
        double* delta;
        for (i=0; i < K; i++){
            delta = update_centroid(k); //updates centroids and calculate the change
            if (norm_calc(delta) < EPSILON) {
                cnt++;
            }
        }
        if (cnt == K){break;}
        iter_num++;
    }
    
    WriteOutToFile();
}




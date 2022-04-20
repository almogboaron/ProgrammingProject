#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>


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
int *clusters;
int **clustersAloc;
double EPSILON = 0.001;




//Function Initalizeation Centroids -> Matrix Formation as first k data points 
void Init_Centroids(int k , int d){
    Centroids = calloc(k*d,sizeof(float));
    CenetroidAloc = calloc(k,sizeof(int));
    //Initializing Centroids Mat
    int i;
    for(i=0 ; i<k*d ;i++){
        Centroids[i] = data[i];
    }
    //Initializing Pointer to Centroids Mat (k*d)
    for(i=0; i < k ; i++){
        CenetroidAloc[i] = Centroids + i*d ;
    }
}

//Function Read Files and initialize Matrix (n*d) of Information 
void init_DataMat( char *filename ,int n , int d){
    float a;
    int i=0;
    int j=0;
    char line , number ;
    
    //Call location
    data = calloc(n*d,sizeof(float));
    dataAloc = calloc(n,sizeof(int));
    
    //Initializing Pointer to Data Mat (k*d) (Or Vectors X1,...Xn == R1....Rn)
    for( i=0; i<n ; i++){
        dataAloc[i] = data + i*d ; 
    }

    //Initialize Data Matrix( Data Array n*d); 
    File *fp = fopen(filename,"r");
    Assert(fp!=NULL);
    while(!feof(fp)){
        fscanf(fp,"%s",line);
        number = strtok(line,",");
        while(number!=NULL){
            data[i++] = (float)atof(number);
            number = strtok(NULL,",");
        }
    }
}

// Function Cols Dims. 
int numOfCols(char *filename){
    int countcols=0;
    char *token,*line;
    File *fp = fopen(filename,"r");
    Assert(fp!=NULL);
    fscanf(fp,"%s",line);
    fclose(fp);
    token = strtok(line,",");
    while (token != NULL){
        countcols++;
        token = strtok(NULL,",");
    }
    return countcols;

}

//Counting Rows Dim
int numbOfRows(char *filename){
    char *line;
    int countrows=0;
    File *fp = fopen(filename,"r");
    Assert(fp!=NULL);
    while(feof(fp)){
        fscanf(fp,"%s",line);
        ++countrows;
    }
    fclose(fp);
    return countrows;
}

double norm_calc(double *delta){
    double s = 0;
    int i;
    for(i=0; i<d; i++){
        s += pow(delta[i],2);
    }
    return sqrt(s);
}

void assign(int *datapoint){
    double val = 1.7976931348623157E+308;
    int idx;
    int j;
    for(j=0; j<K; j++){
        double delta[d];
        int i;
        for(i=0; i<d; i++){
            delta[i] = datapoint[i] - CenetroidAloc[j][i];
        }
        double norm = norm_calc(delta);
        if (norm < val){
            val = norm;
            idx = j;
        }
    }
    int i;
    for (i=0; i<n; i++) {
        if (clustersAloc[idx][i] == 0) {
            clustersAloc[idx][i] = *datapoint;
        }
    }
}

int update_centroids(int k){
    return 0;
}


int main(int argc, char *argv[]) {
    if (argc == 2) { max_iter = 200; }
    else { max_iter = atoi(argv[2]); }
    K = atoi(argv[0]);
    filename = argv[1];
    n = numbOfRows(filename);
    d = numOfCols(filename);
    init_DataMat(filename,n,d);
    Init_Centroids(K,d);
    clusters = calloc(n*K,sizeof(double));
    clustersAloc = calloc(n,sizeof (double *));
    int i;
    for (i=0; i<n; i++){
        clustersAloc[i] = clusters + i*d;
        int j;
        for (j=0; j<n; j++){
            clustersAloc[i][j] = 0;
        }
    }
    int iter_num = 0;
    while (iter_num < max_iter){
        for(i=0; i < n; i++){
            assign(&data[i]);
        }
        int cnt = 0;
        int k;
        for (k=0; k < K; k++){
            double *delta;
            *delta = update_centroids(k); //updates centroids and calculate the change
            if (norm_calc(delta) < EPSILON) {
                cnt++;
            }
        }
        if (cnt == K){break;}
        iter_num++;
    }
}




#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>


//Declorations
int *Centroids;
int **CenetroidAloc;
int *data;
int **dataAloc;


//Function Initalizeation Centroids -> Matrix Formation as first k data points 
void Init_Centroids(int k , int d){
    centroids = calloc(k*d,sizeof(float));
    centroidAloc = calloc(k,sizof(int));
    //Initializing Centroids Mat
    for( i=0 ; i<k*d ;i++){
        centroids[i] = data[i];
    }
    //Initializing Pointer to Centroids Mat (k*d)
    for( i=0 , i<k , i++){
        centroidAloc[i] = centroids + i*d ; 
    }
}

//Function Read Files and initialize Matrix (n*d) of Information 
void init_DataMat( char *filename ,int n , int d){
    float a;
    int i=0;
    int j=0;
    *char line , number ;
    
    //Call location
    data = calloc(n*d,sizeof(float));
    dataAloc = calloc(n,sizeof(int));
    
    //Initializing Pointer to Data Mat (k*d) (Or Vectors X1,...Xn == R1....Rn)
    for( i=0 , i<n , i++){
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
        countcols++
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
        fscanf(fp,"%s",line)
        ++countrows;
    }
    fclose(fp);
    return countrow;
}

//Function Argmin

//Update Centroids

// Euclidean delta Centruids -> True or False 


















int main(int k , c[] filename, Optional Iter )
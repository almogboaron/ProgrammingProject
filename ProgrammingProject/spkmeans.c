#define PY_SSIZE_T_CLEAN

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <float.h>

/*Declorations*/
int K=0;
char* filename;
int max_iter = 300;
int n;
int d;
double EPSILON = 0;
const double eps =0.00001;
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

/*All eigenVectors*/
double* V;
double** VAloc;

/*k eigenVectors*/
double* U;
double** UAloc;

/*Normalized by Rows U -> T*/
double* T;
double** TAloc;

/*Eigen Values and there indexes*/
double* EigenValues;
int* indexes;

/*define assert with message*/
#define assert__(x) for (;!(x);assert(x))

/* -------------------------------- Handling Data ------------------- */
double MAX(double x,double y){ 
  return (x >= y) ? x : y;
  }

/*Prints Functions for Arrays and Matrix*/
void print_output(double** mat,int r,int c){
    int i,j;
    for(i=0; i<r ;i++){
        for(j=0; j<c; j++){
            /*Not Last Column*/
            if(j<c-1){printf("%.4f,",mat[i][j]);}
            else    {printf("%.4f",mat[i][j]);}
        }
        printf("%s","\n");
    }
}

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
    fclose(fp);
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
    for(k=0;k<d;k++){
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
            if(i==j){continue;}
            wamAloc[i][j] = exp(-(norm2(i,j)/2));
        }
    }
    free(data);free(dataAloc);
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
    free(wamAloc);free(wam);
}

/*---------------The Normalized Graph Laplacian---------*/
/*Matrix Multiplacation*/
void MatMulti(double** mat1,double** mat2 , double** matRes){
    int i,j,k;
    double sum=0;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            for(k=0;k<n;k++){
                sum+=mat1[i][k]*mat2[k][j];
            }
            matRes[i][j] = sum;
            sum=0;
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
    free(ddg);free(ddgAloc);
}

/*---------------Jacobi Algorithem---------*/
/*Jacobi Helper Functions*/
/*Initalize PAloc*/
void initPAloc(double s,double c, int piv_i , int piv_j, double** PAloc){
    int i,j;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(i==j){
                PAloc[i][j]=1;
            }else{
                PAloc[i][j]=0;
            }
        } 
    }
    PAloc[piv_i][piv_i] = c;
    PAloc[piv_j][piv_j] = c;
    PAloc[piv_i][piv_j] = -s;
    PAloc[piv_j][piv_i] = s;
}

/*Sign Function*/
int sign(double teta){
    if(teta>=0){
        return 1;
    }else{
        return -1;
    }
}

/*Off function*/
double offMat(double** mat){
    int i,j;
    double sum=0;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(i!=j){
                sum+=pow(mat[i][j],2);
            }
        }
    }
    return sum;
}

/*Compare function - Reversed*/
int cmpfuncindex (const void * a, const void * b) {
   return ( EigenValues[*(int*)b] - EigenValues[*(int*)a] );
}

/*Compare function - Reversed*/
int cmpfunc (const void * a, const void * b) {
   return ( *(double*)b - *(double*)a );
}

/*Eigengap Heuristic*/
void EigengapHeuristic(){
    double sigma;
    int i;
    double max = -DBL_MAX;
    
/*Finidng k*/
    for(i=0; i<floor(n/2);i++){
        sigma = abs(EigenValues[i]-EigenValues[i+1]);
        if(max < sigma){
            max = sigma;
            K = i+1;  /*Could Be A BUG************************/
        }
    }
}

/*Jacobi Main Function - Creates V*/
void jacobiInit(){
    int k,r,i,piv_i,j,piv_j;
    double teta, t, s, c, temp, maxVal, off_prev,off_after, *P, **PAloc, *res, **resAloc, *tempmat, **tempMat;

    /* Creating mat*/
    P = calloc(n*n,sizeof(double));
    assert__(P!=NULL){
        printf("An Error Has Occurred");
    }
    PAloc = calloc(n,sizeof(double*));
    assert__(PAloc!=NULL){
        printf("An Error Has Occurred");
    }
    V = calloc(n*n,sizeof(double));
    assert__(V!=NULL){
        printf("An Error Has Occurred");
    }
    VAloc = calloc(n,sizeof(double*));
    assert__(VAloc!=NULL){
        printf("An Error Has Occurred");
    }
    res = calloc(n*n,sizeof(double));
    assert__(res!=NULL){
        printf("An Error Has Occurred");
    }
    resAloc = calloc(n,sizeof(double*));
    assert__(resAloc!=NULL){
        printf("An Error Has Occurred");
    }
    for(i=0;i<n;i++){
        PAloc[i] = P + i*n;
        VAloc[i] = V +i*n;
        VAloc[i][i] = 1;
        resAloc[i] = res +i*n;
        resAloc[i][i] = 1;
    }
    
    off_prev = offMat(lnormAloc);
    
    /*Diagonalization*/
    for(k=0;k<100;k++){
        /*Pivot*/
        maxVal = -DBL_MAX;
        for(i=0;i<n;i++){
            for(j=0;j<i;j++){
                if(MAX(maxVal,fabs(lnormAloc[i][j]))==fabs(lnormAloc[i][j])){
                    piv_i = i; 
                    piv_j = j;
                    maxVal = fabs(lnormAloc[i][j]);
                }
            }
        }
        
        printf("maxVal %f",maxVal);
        printf("piv_i = %d, piv_j = %d\n",piv_i,piv_j);
        printf("lnormIter%d before change\n\n",k);
        print_output(lnormAloc,n,n);

        /*c And s*/
        teta = (lnormAloc[piv_j][piv_j]-lnormAloc[piv_i][piv_i])/(2*lnormAloc[piv_i][piv_j]);
        t = sign(teta)/(fabs(teta) + sqrt(pow(teta,2)+1)); 
        c = 1/sqrt(pow(t,2)+1);
        s = c*t;
        printf("teta %f,t %f,c %f,s %f\n",teta,t,c,s);
      
        /*Transform lnorm | Temp used here to avoid overiding*/
        for(r=0;r<n;r++){
            if((r!=piv_i)&&(r!=piv_j)){
                temp = lnormAloc[r][piv_i];
                lnormAloc[r][piv_i] = c*temp - s*lnormAloc[r][piv_j];
                lnormAloc[r][piv_j] = c*lnormAloc[r][piv_j] + s*temp;
                lnormAloc[piv_i][r] = lnormAloc[r][piv_i];
                lnormAloc[piv_j][r] = lnormAloc[r][piv_j];
            }
        }
        temp = lnormAloc[piv_i][piv_i];
        printf("%f,%f,%f\n",lnormAloc[piv_i][piv_i],lnormAloc[piv_j][piv_j],lnormAloc[piv_i][piv_j]);
        
        lnormAloc[piv_i][piv_i] = pow(c,2)*temp + pow(s,2)*lnormAloc[piv_j][piv_j] - 2*s*c*lnormAloc[piv_i][piv_j];
        lnormAloc[piv_j][piv_j] = pow(s,2)*temp + pow(c,2)*lnormAloc[piv_j][piv_j] + 2*s*c*lnormAloc[piv_i][piv_j];
        lnormAloc[piv_i][piv_j] = 0;
        lnormAloc[piv_j][piv_i] = 0;
        printf("%f,%f,%f\n",lnormAloc[piv_i][piv_i],lnormAloc[piv_j][piv_j],lnormAloc[piv_i][piv_j]);
        
        /*Initialize P*/
        initPAloc(s, c, piv_i, piv_j, PAloc);
        
        /*Update V*/
        MatMulti(VAloc,PAloc,resAloc);
        tempMat = VAloc;
        tempmat = V;
        VAloc = resAloc;
        V = res;
        resAloc = tempMat;
        res = tempmat;

        /*Checking Convergence*/         
        off_after = offMat(lnormAloc);
        if( (off_prev - off_after) <= eps){break;}
        off_prev = off_after;
    }

    /*Sorting EigenValues and indexes for matrix U*/
    indexes = calloc(n,sizeof(int));
    assert__(indexes!=NULL){
        printf("An Error Has Occurred");
    }
    EigenValues = calloc(n,sizeof(double));
    assert__(EigenValues!=NULL){
        printf("An Error Has Occurred");
    }
    
    /*Initializing*/
    for(i=0;i<n;i++){
        EigenValues[i] = lnormAloc[i][i];
        indexes[i] = i;
    }

    /*Sorting indexes and values*/
    qsort(indexes, n, sizeof(int), cmpfuncindex);
    qsort(EigenValues,n,sizeof(double),cmpfunc);

    /*EigengapHeuristic*/
    if(K == 0){
        EigengapHeuristic();
    }
    free(P);free(PAloc);
    free(res);free(resAloc);
    free(lnorm);free(lnormAloc);
}

/*---------------------------PYTHON C API-----------------------------*/


/*Main Function of Kmeans Algorithem .*/
void kmeans() {
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

/*Free Memory*/
free(ClusterSum);
free(ClusterSumAloc);
free(data);
free(dataAloc);
free(count);
}

/*spk Implemntation for Python*/
void spk() {
    int i,j;
    double norma;

    /*getting the data from the file*/
    read_file();
    
    /*Asserts*/
    assert__(K<=n){
        printf("Invalid Input!");
    }
    
    /*step 1*/
    wamInit();

    /*step 2*/
    ddgInit();
    lnormInit();

    /*step 3*/
    jacobiInit();

    /*step 4 - renormalize U's rows - NEED TO ADD*/
    U = calloc(n*K,sizeof(double));
    assert__(U!=NULL){
        printf("An Error Has Occurred");
    }
    UAloc = calloc(n,sizeof(double*));
    assert__(UAloc!=NULL){
        printf("An Error Has Occurred");
    }
    /*Initialize UAloc*/
    for(i=0;i<n;i++){
        UAloc[i] = U +i*K;
    }
    
    /*Initialize choesen EigenVectors by indexes*/
    for(i=0;i<n;i++){
        for(j=0;j<K;j++){
            UAloc[i][j] = VAloc[i][indexes[j]]; 
        }
    }
    free(indexes);
    free(V);free(VAloc);
    

    /*Normaliazing*/
    T = calloc(n*K,sizeof(double));
    assert__(T!=NULL){
        printf("An Error Has Occurred");
    }
    TAloc = calloc(n,sizeof(double*));
    assert__(TAloc!=NULL){
        printf("An Error Has Occurred");
    }
    for(i=0;i<n;i++){
        TAloc[i] = T +i*K;
    }
    
    /*Norma Function*/
    for(i=0;i<n;i++){
        norma =0;
        for(j=0;j<K;j++){
            norma+=pow(UAloc[i][j],2);
        }
        norma = sqrt(norma);
        for(j=0;j<K;j++){
            TAloc[i][j] = UAloc[i][j]/norma;
        }
    }

    /*step 5 - passing the data back to python*/
    free(U);free(UAloc);
}


void print_outputV(double** mat,int r,int c){
    int i,j;
    for(i=0; i<r ;i++){
        for(j=0; j<c; j++){
            /*Take the Index of the currect EigenValue*/
            if(j<c-1){printf("%.4f,",mat[i][indexes[j]]);}
            else    {printf("%.4f",mat[i][indexes[j]]);}
        }
        printf("%s","\n");
    }
}

void print_outputArr(double* arr,int r){
    int i;
    for(i=0; i<r; i++){
        if(i<r-1){printf("%.4f,",arr[i]);}
        else    {printf("%.4f\n",arr[i]);}
    }
}

/*Main Function for C*/
int main(int argc, char *argv[]) {
    assert__(argc == 3){
        printf("Invalid Input!");
    }
    filename = argv[2];
    read_file();

    if(strcmp(argv[1],"wam") == 0){
        wamInit();
        print_output(wamAloc,n,n);
        free(wam);free(wamAloc);
    }
    else if(strcmp(argv[1],"ddg") == 0){
        wamInit();
        ddgInit();
        print_output(ddgAloc,n,n);
        free(ddg);free(ddgAloc);
    }
    else if(strcmp(argv[1],"lnorm") == 0){
        wamInit();
        ddgInit();
        lnormInit();
        print_output(lnormAloc,n,n);
        free(lnorm);free(lnormAloc);
    }
    else if(strcmp(argv[1],"jacobi")==0){
        wamInit();
        ddgInit();
        lnormInit();
        jacobiInit();
        print_outputArr(EigenValues,n);
        print_outputV(VAloc,n,n);
        free(EigenValues);
        free(indexes);
        free(V);free(VAloc);
    }
    else{printf("Invalid Input!");}
    return 0;
}
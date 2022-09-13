/*vars*/
int K;
char* filename;
double* data;
double** dataAloc;
double* wam;
double** wamAloc;
double* ddg;
double** ddgAloc;
double* lnorm;
double** lnormAloc;
double* V;
double** VAloc;
double* EigenValues;
int* indexes;
int n;
double* Centroids;
double** CenetroidAloc;
double* T;
double** TAloc;
int d;
double* data;
double** dataAloc;


/*functions*/

void read_file();
void wamInit();
void print_output();
void print_outputV();
void print_outputArr();
void ddgInit();
void lnormInit();
void jacobiInit();
void kmeans();
void spk();

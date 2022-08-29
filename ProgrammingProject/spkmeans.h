/*vars*/
PyObject* data_arr;
PyObject* centroid_arr;
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
int* indexs;
int n;

/*functions*/
void read_file();
void wamInit();
void print_output();
void print_outputV();
void print_outputArr();
void ddgInit();
void lnormInit();
void jacobiInit();
PyObject* kmeans();
PyObject* spk();

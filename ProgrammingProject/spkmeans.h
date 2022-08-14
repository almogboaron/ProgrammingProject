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

/*functions*/
static PyObject* kmeans();
static PyObject* spk();
void read_file();
void wamInit();
void print_output();
void ddgInit();
void lnormInit();


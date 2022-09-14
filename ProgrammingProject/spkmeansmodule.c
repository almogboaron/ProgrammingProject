#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <float.h>


/*Python Objects*/
PyObject* centroid_arr;
PyObject* pyCentroids;
PyObject* pyT;

#define assert__(x) for (;!(x);assert(x))
/*Function Initialization Centroids -> Matrix Formation as first k data points */
void Init_Centroids(){
    int i;
    int j;
    Centroids = calloc(K*K,sizeof(double));
    assert__(Centroids!=NULL){
        printf("An Error Has Occurred");
    }
    CenetroidAloc = calloc(K,sizeof(double*));
    assert__(CenetroidAloc!=NULL){
        printf("An Error Has Occurred");
    }
    /*Initializing Centroids Mat*/
    for(i=0;i<K;i++){
        for(j=0;j<K;j++){
            Centroids[i*d+j] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(centroid_arr,i),j));
        }
    }
    /*Initializing Pointer to Centroids Mat (k*d)*/
    for(i=0; i < K ; i++){
        CenetroidAloc[i] = Centroids + i*K ;
    }
}

/*Create PyCentroids matrix*/
static PyObject* cMatToPy(double** mat, int rows, int cols){
    int i;
    int j;
    PyObject* pyMat = PyList_New(0);
    for(i=0;i<rows;i++){
        PyObject* PyPoint = PyList_New(0);
        for(j=0;j<cols;j++){
            if(PyList_Append(PyPoint,PyFloat_FromDouble(mat[i][j]))==-1){
                return NULL;
            }
        }
        if(PyList_Append(pyMat,PyPoint)==-1){
            return NULL;
        }
    }
    return pyMat;
}

static PyObject* fit(PyObject* Py_UNUSED(self), PyObject* args){
    /* data and initial centroids*/
    if (!PyArg_ParseTuple(args,"O!",&PyList_Type, &centroid_arr)){
        return NULL;
    }
    d = K;
    Init_Centroids();
    kmeans();
    pyCentroids = cMatToPy(CenetroidAloc,K,K);
    free(Centroids);
    free(CenetroidAloc);
    return pyCentroids;
}

static PyObject* spkC(PyObject* self, PyObject* args){
    if(!PyArg_ParseTuple(args,"si",&filename,&K)){
        return NULL;
    }
    spk();
    pyT = cMatToPy(TAloc,n,K);
    return pyT;
}

static PyObject* wamC(PyObject* self, PyObject* args){
    if(!PyArg_ParseTuple(args,"s",&filename)){
        return NULL;
    }
    read_file();
    wamInit();
    print_output(wamAloc,n,n);
    free(wam);
    free(wamAloc);
    Py_RETURN_NONE;
}

static PyObject* ddgC(PyObject* self, PyObject* args){
    if(!PyArg_ParseTuple(args,"s",&filename)){
        return NULL;
    }
    read_file();
    wamInit();
    ddgInit();
    print_output(ddgAloc,n,n);
    free(wam);free(wamAloc);
    free(ddg);free(ddgAloc);
    Py_RETURN_NONE;
}

static PyObject* lnormC(PyObject* self, PyObject* args){
    if(!PyArg_ParseTuple(args,"s",&filename)){
        return NULL;
    }
    read_file();
    wamInit();
    ddgInit();
    read_file();
    lnormInit();
    print_output(lnormAloc,n,n);
    free(lnorm);
    free(lnormAloc);
    Py_RETURN_NONE;
}

static PyObject* jacobiC(PyObject* self,PyObject* args){
    if(!PyArg_ParseTuple(args,"s",&filename)){
        return NULL;
    }
    read_file();
    lnormAloc = dataAloc;
    jacobiInit();
    print_outputArr(EigenValues,n);
    print_output(VAloc,n,n);
    free(EigenValues);
    free(V);free(VAloc);
    Py_RETURN_NONE;
}

/* declaring the kmeans function */
static PyMethodDef capiMethods[] = {
        {"fit", (PyCFunction) fit, METH_VARARGS, PyDoc_STR("Kmeans Algorithem")},
        {"spkC", (PyCFunction) spkC, METH_VARARGS, PyDoc_STR("Spk Algorithem")},
        {"wamC", (PyCFunction) wamC, METH_VARARGS, PyDoc_STR("Weighted Adjacency Matrix")},
        {"ddgC", (PyCFunction) ddgC, METH_VARARGS, PyDoc_STR("Diagonal Degree Matrix")},
        {"lnormC",(PyCFunction) lnormC, METH_VARARGS, PyDoc_STR("Normalized Graph Laplacian")},
        {"jacobiC",(PyCFunction) jacobiC, METH_VARARGS, PyDoc_STR("Jacobi Algorithm")},
        {NULL,NULL,0,NULL}
};

/*defining the module*/
static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT, "spkmeansmodule", NULL, -1, capiMethods
};

/* creating the module */
PyMODINIT_FUNC
PyInit_spkmeansmodule(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if(!m) {
        return NULL;
    }
    return m;
}
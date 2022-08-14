#include <Python/Python.h>
#include "spkmeans.h"


static PyObject* fit(PyObject* Py_UNUSED(self), PyObject* args){
    /* data and initial centroinds*/
    if (!PyArg_ParseTuple(args,"O!O!", &PyList_Type ,&data_arr,&PyList_Type, &centroid_arr)){
        return NULL;
    }
    return kmeans();
}

static PyObject* spkC(PyObject* self, PyObject* args){
    if(!PyArg_ParseTuple(args,"Oi",&filename,&K)){
        return NULL;
    }
    return spk();
}

static PyObject* wamC(PyObject* self, PyObject* args){
    if(!PyArg_ParseTuple(args,"O",&filename)){
        return NULL;
    }
    read_file();
    wamInit();
    print_output(wamAloc);
    free(data);
    free(dataAloc);
    free(wam);
    free(wamAloc);
    return 0;
}

static PyObject* ddgC(PyObject* self, PyObject* args){
    if(!PyArg_ParseTuple(args,"O",&filename)){
        return NULL;
    }
    read_file();
    ddgInit();
    print_output(ddgAloc);
    free(data);
    free(dataAloc);
    free(ddg);
    free(ddgAloc);
    return 0;
}

static PyObject* lnormC(PyObject* self, PyObject* args){
    if(!PyArg_ParseTuple(args,"O",&filename)){
        return NULL;
    }
    read_file();
    lnormInit();
    print_output(lnormAloc);
    free(data);
    free(dataAloc);
    free(lnorm);
    free(lnormAloc);
    return 0;
}

/* declaring the kmeans function */
static PyMethodDef capiMethods[] = {
        {"fit", (PyCFunction) fit, METH_VARARGS, PyDoc_STR("kmeans algorithem")},
        {"spkC", (PyCFunction) spkC, METH_VARARGS, PyDoc_STR("spk algorithem")},
        {"wamC", (PyCFunction) wamC, METH_VARARGS, PyDoc_STR("Weighted Adjacency Matrix")},
        {"ddgC", (PyCFunction) ddgC, METH_VARARGS, PyDoc_STR("Diagonal Degree Matrix")},
        {"lnormC",(PyCFunction) lnormC, METH_VARARGS, PyDoc_STR("Normalized Graph Laplacian")},
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
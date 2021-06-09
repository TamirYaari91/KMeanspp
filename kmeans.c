#define PY_SSIZE_T_CLEAN

//#include <python3.7/Python.h>
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>


double distance(double[], double *, int, int);

void set_cluster(int, int, double *, double *, int);

double *cluster_mean(int, int *, double *, int, int);

int update_centroids(int, int, double *, double *, int);

int equal(double *, double *, int);

static PyObject *kmeans(int, int, int, int, PyObject *, PyObject *, int, int);

static PyObject *fit(PyObject *self, PyObject *args);

#define FUNC(_flag, _name, _docstring) { #_name, (PyCFunction)_name, _flag, PyDoc_STR(_docstring) }


static PyMethodDef capiMethods[] = {
        FUNC(METH_VARARGS, fit, "calculate k-means"),
        {NULL, NULL, 0, NULL}   /* sentinel */
};


static struct PyModuleDef moduleDef = {
        PyModuleDef_HEAD_INIT, "mykmeanssp", NULL, -1, capiMethods};

PyMODINIT_FUNC
PyInit_mykmeanssp(void) {
    return PyModule_Create(&moduleDef);
}

static PyObject *fit(PyObject *self, PyObject *args) {
    int k;
    int max_iter;
    int num_of_lines;
    int dim;
    PyObject *centroids_py;
    PyObject *points_to_cluster_py;
    int centroids_length;
    int points_to_cluster_length;

    if (!PyArg_ParseTuple(args, "llllOOll:fit", &k, &max_iter, &num_of_lines, &dim, &centroids_py,
                          &points_to_cluster_py,
                          &centroids_length, &points_to_cluster_length)) {
        return NULL;
    }

    return Py_BuildValue("O", kmeans(k, max_iter, num_of_lines, dim, centroids_py, points_to_cluster_py,
                                     centroids_length, points_to_cluster_length));
}

static PyObject *kmeans(int k, int max_iter, int num_of_lines, int dim, PyObject *centroids_py,
                        PyObject *points_to_cluster_py, int centroids_length, int points_to_cluster_length) {

    double *centroids;
    double *points_to_cluster;
    int i;

    if (centroids_length < 0 || points_to_cluster_length < 0) {
        return NULL;
    }

    centroids = (double *) calloc(centroids_length, sizeof(double));
    assert(centroids != NULL && "Problem in allocating centroids memory");

    points_to_cluster = (double *) calloc(points_to_cluster_length, sizeof(double));
    assert(points_to_cluster != NULL && "Problem in allocating points_to_cluster memory");

    for (i = 0; i < centroids_length; i++) {
        PyObject *item;
        item = PyList_GetItem(centroids_py, i);
        centroids[i] = PyFloat_AsDouble(item);
    }
    for (i = 0; i < points_to_cluster_length; i++) {
        PyObject *item;
        item = PyList_GetItem(points_to_cluster_py, i);
        points_to_cluster[i] = PyFloat_AsDouble(item);
    }

    int changed;
    int iters;

    iters = 0;
    while (1) {
        for (i = 0; i < num_of_lines; ++i) {
            set_cluster(i, k, points_to_cluster, centroids, dim);
        }
        changed = update_centroids(k, num_of_lines, points_to_cluster, centroids, dim);
        iters++;
        if (changed == 0 || iters == max_iter) {
            break;
        }
    }

    PyObject *list = PyList_New(centroids_length);
    for (i = 0; i < centroids_length; i++) {
        PyList_SetItem(list, i, PyFloat_FromDouble(centroids[i]));
    }
    free(centroids);
    free(points_to_cluster);
    return list;
}

double distance(double *p, double *centroids, int cluster, int dim) {
    double d = 0;
    int i;

    for (i = 0; i < dim; ++i) {
        double multi;
        multi = (p[i] - *((centroids + (cluster * dim) + i)));
        d += multi * multi;
    }
    return d;
}

void set_cluster(int p_index, int k, double *point_to_cluster, double *centroids, int dim) {
    int i;
    int min_index = 0;
    double *distances;
    distances = (double *) calloc(k, sizeof(double));
    assert(distances);

    for (i = 0; i < k; ++i) {
        distances[i] = distance((point_to_cluster + p_index * (dim + 1)), centroids, i, dim);
        if (distances[i] < distances[min_index]) {
            min_index = i;
        }
    }
    point_to_cluster[(p_index * (dim + 1) + dim)] = min_index;
}

double *cluster_mean(int cluster, int *c2p, double *p2c, int dim, int num_of_points) {
    int size;
    double val;
    double p_val;
    int i, j;
    static double *center;
    center = (double *) calloc(dim, sizeof(double));
    assert(center);
    size = 0;
    val = 0.0;

    for (i = 0; i < dim; ++i) {
        for (j = 0; j < num_of_points; ++j) {

            int p_index = *(c2p + (cluster * num_of_points) + j);
            if (p_index < 0) {
                break;
            }
            size++;
            p_val = *(p2c + p_index * (dim + 1) + i);
            val += p_val;
        }
        center[i] = val / size;
        size = 0;
        val = 0.0;
    }

    return center;
}

int update_centroids(int k, int num_of_points, double *p2c, double *centroids, int dim) {
    int *c2p;
    int i, j, changed;
    c2p = (int *) calloc(k * num_of_points, sizeof(int));
    assert(c2p);

    for (i = 0; i < k; ++i) {
        for (j = 0; j < num_of_points; ++j) {
            c2p[i * num_of_points + j] = -1;
        }
    }
    for (i = 0; i < num_of_points; ++i) {
        int cluster = p2c[i * (dim + 1) + dim];
        for (j = 0; j < num_of_points; ++j) {
            if (c2p[(cluster * num_of_points) + j] == -1) {
                c2p[(cluster * num_of_points) + j] = i;

                break;
            }
        }
    }
    changed = 0;
    for (i = 0; i < k; ++i) {
        double *new_centroid = cluster_mean(i, c2p, p2c, dim, num_of_points);
        if (equal((centroids + dim * i), new_centroid, dim) == 0) {
            for (j = 0; j < dim; ++j) {
                *(centroids + dim * i + j) = new_centroid[j];
            }
            changed = 1;
        }
        free(new_centroid);
    }
    free(c2p);
    return changed;
}

int equal(double *arr1, double *arr2, int dim) {
    int i;
    for (i = 0; i < dim; ++i) {
        if (arr1[i] != arr2[i]) {
            return 0;
        }
    }
    return 1;
}
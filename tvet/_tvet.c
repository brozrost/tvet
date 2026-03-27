// clang -O3 -shared -fPIC \
$(python3-config --includes) \
$(python3-config --embed --ldflags) \
_tvet.c \
csrc/vector_math.c csrc/intersect_AB_t.c csrc/bounding_box.c csrc/shadowing.c \
-o _tvet$(python3-config --extension-suffix)

#include <Python.h>
#include <stdbool.h>
#include <stddef.h>

#include "csrc/vector_math.h"
#include "csrc/intersect_AB_t.h"
#include "csrc/bounding_box.h"
#include "csrc/shadowing.h"

static PyObject *py_dotProduct3D(PyObject *self, PyObject *args) {
    PyObject *x = NULL;
    PyObject *y = NULL;
    if (!PyArg_ParseTuple(args, "OO", &x, &y)) return NULL;

    Py_buffer x_buf;
    Py_buffer y_buf;

    if (PyObject_GetBuffer(x, &x_buf, PyBUF_FORMAT) != 0) return NULL;
    if (PyObject_GetBuffer(y, &y_buf, PyBUF_FORMAT) != 0) {
        PyBuffer_Release(&x_buf);
        return NULL;
    }

    double result = dotProduct3D((const double *)x_buf.buf, (const double *)y_buf.buf);

    PyBuffer_Release(&x_buf);
    PyBuffer_Release(&y_buf);

    return PyFloat_FromDouble(result);
}

static PyObject *py_crossProduct3D(PyObject *self, PyObject *args) {
    PyObject *x = NULL;
    PyObject *y = NULL;
    PyObject *out = NULL;
    if (!PyArg_ParseTuple(args, "OOO", &x, &y, &out)) return NULL;

    Py_buffer x_buf;
    Py_buffer y_buf;
    Py_buffer out_buf;

    if (PyObject_GetBuffer(x, &x_buf, PyBUF_FORMAT) != 0) return NULL;
    if (PyObject_GetBuffer(y, &y_buf, PyBUF_FORMAT) != 0) {
        PyBuffer_Release(&x_buf);
        return NULL;
    }
    if (PyObject_GetBuffer(out, &out_buf, PyBUF_WRITABLE | PyBUF_FORMAT) != 0) {
        PyBuffer_Release(&x_buf);
        PyBuffer_Release(&y_buf);
        return NULL;
    }

    crossProduct3D((const double *)x_buf.buf, (const double *)y_buf.buf, (double *)out_buf.buf);

    PyBuffer_Release(&x_buf);
    PyBuffer_Release(&y_buf);
    PyBuffer_Release(&out_buf);

    Py_RETURN_NONE;
}

static PyObject *py_normalizeVector3D(PyObject *self, PyObject *args) {
    PyObject *v = NULL;
    if (!PyArg_ParseTuple(args, "O", &v)) return NULL;

    Py_buffer v_buf;
    if (PyObject_GetBuffer(v, &v_buf, PyBUF_WRITABLE | PyBUF_FORMAT) != 0) return NULL;

    int ok = normalizeVector3D((double *)v_buf.buf);

    PyBuffer_Release(&v_buf);

    if (ok) Py_RETURN_TRUE;
    Py_RETURN_FALSE;
}

static PyObject *py_intersect_AB_t(PyObject *self, PyObject *args) {
    PyObject *A = NULL;
    PyObject *B = NULL;
    PyObject *t = NULL;
    PyObject *C = NULL;
    if (!PyArg_ParseTuple(args, "OOOO", &A, &B, &t, &C)) return NULL;

    Py_buffer A_buf;
    Py_buffer B_buf;
    Py_buffer t_buf;
    Py_buffer C_buf;

    if (PyObject_GetBuffer(A, &A_buf, PyBUF_FORMAT) != 0) return NULL;
    if (PyObject_GetBuffer(B, &B_buf, PyBUF_FORMAT) != 0) {
        PyBuffer_Release(&A_buf);
        return NULL;
    }
    if (PyObject_GetBuffer(t, &t_buf, PyBUF_FORMAT) != 0) {
        PyBuffer_Release(&A_buf);
        PyBuffer_Release(&B_buf);
        return NULL;
    }
    if (PyObject_GetBuffer(C, &C_buf, PyBUF_WRITABLE | PyBUF_FORMAT) != 0) {
        PyBuffer_Release(&A_buf);
        PyBuffer_Release(&B_buf);
        PyBuffer_Release(&t_buf);
        return NULL;
    }

    bool has_solution = false;
    const double (*tri)[3] = (const double (*)[3])t_buf.buf;

    intersect_AB_t(
        (const double *)A_buf.buf,
        (const double *)B_buf.buf,
        tri,
        (double *)C_buf.buf,
        &has_solution
    );

    PyBuffer_Release(&A_buf);
    PyBuffer_Release(&B_buf);
    PyBuffer_Release(&t_buf);
    PyBuffer_Release(&C_buf);

    if (has_solution) Py_RETURN_TRUE;
    Py_RETURN_FALSE;
}

static PyObject *py_boundingBox(PyObject *self, PyObject *args) {
    PyObject *faces = NULL;
    Py_ssize_t nof_faces = 0;

    PyObject *nodes = NULL;
    Py_ssize_t nof_nodes = 0;

    PyObject *s = NULL;
    PyObject *boxes = NULL;

    if (!PyArg_ParseTuple(args, "OnOnOO", &faces, &nof_faces, &nodes, &nof_nodes, &s, &boxes)) return NULL;

    Py_buffer faces_buf;
    Py_buffer nodes_buf;
    Py_buffer s_buf;
    Py_buffer boxes_buf;

    if (PyObject_GetBuffer(faces, &faces_buf, PyBUF_FORMAT) != 0) return NULL;
    if (PyObject_GetBuffer(nodes, &nodes_buf, PyBUF_FORMAT) != 0) {
        PyBuffer_Release(&faces_buf);
        return NULL;
    }
    if (PyObject_GetBuffer(s, &s_buf, PyBUF_FORMAT) != 0) {
        PyBuffer_Release(&faces_buf);
        PyBuffer_Release(&nodes_buf);
        return NULL;
    }
    if (PyObject_GetBuffer(boxes, &boxes_buf, PyBUF_WRITABLE | PyBUF_FORMAT) != 0) {
        PyBuffer_Release(&faces_buf);
        PyBuffer_Release(&nodes_buf);
        PyBuffer_Release(&s_buf);
        return NULL;
    }

    boundingBox(
        (const int (*)[3])faces_buf.buf, (size_t)nof_faces,
        (const double (*)[3])nodes_buf.buf, (size_t)nof_nodes,
        (const double *)s_buf.buf,
        (double (*)[4])boxes_buf.buf
    );

    PyBuffer_Release(&faces_buf);
    PyBuffer_Release(&nodes_buf);
    PyBuffer_Release(&s_buf);
    PyBuffer_Release(&boxes_buf);

    Py_RETURN_NONE;
}

static PyObject *py_mu(PyObject *self, PyObject *args) {
    PyObject *normals = NULL;
    Py_ssize_t nof_faces = 0;
    PyObject *s = NULL;
    PyObject *mu_i = NULL;

    if (!PyArg_ParseTuple(args, "OnOO", &normals, &nof_faces, &s, &mu_i)) return NULL;

    Py_buffer normals_buf;
    Py_buffer s_buf;
    Py_buffer mu_i_buf;

    if (PyObject_GetBuffer(normals, &normals_buf, PyBUF_FORMAT) != 0) return NULL;
    if (PyObject_GetBuffer(s, &s_buf, PyBUF_FORMAT) != 0) {
        PyBuffer_Release(&normals_buf);
        return NULL;
    }
    if (PyObject_GetBuffer(mu_i, &mu_i_buf, PyBUF_WRITABLE | PyBUF_FORMAT) != 0) {
        PyBuffer_Release(&normals_buf);
        PyBuffer_Release(&s_buf);
        return NULL;
    }

    mu(
        (const double (*)[3])normals_buf.buf,
        (size_t)nof_faces,
        (const double *)s_buf.buf,
        (double *)mu_i_buf.buf
    );

    PyBuffer_Release(&normals_buf);
    PyBuffer_Release(&s_buf);
    PyBuffer_Release(&mu_i_buf);

    Py_RETURN_NONE;
}

static PyObject *py_non(PyObject *self, PyObject *args) {
    PyObject *mu_i = NULL;
    PyObject *mu_e = NULL;
    Py_ssize_t nof_faces = 0;
    PyObject *nu_i = NULL;
    PyObject *nu_e = NULL;

    if (!PyArg_ParseTuple(args, "OOnOO", &mu_i, &mu_e, &nof_faces, &nu_i, &nu_e)) return NULL;

    Py_buffer mu_i_buf;
    Py_buffer mu_e_buf;
    Py_buffer nu_i_buf;
    Py_buffer nu_e_buf;

    if (PyObject_GetBuffer(mu_i, &mu_i_buf, PyBUF_FORMAT) != 0) return NULL;
    if (PyObject_GetBuffer(mu_e, &mu_e_buf, PyBUF_FORMAT) != 0) {
        PyBuffer_Release(&mu_i_buf);
        return NULL;
    }
    if (PyObject_GetBuffer(nu_i, &nu_i_buf, PyBUF_WRITABLE | PyBUF_FORMAT) != 0) {
        PyBuffer_Release(&mu_i_buf);
        PyBuffer_Release(&mu_e_buf);
        return NULL;
    }
    if (PyObject_GetBuffer(nu_e, &nu_e_buf, PyBUF_WRITABLE | PyBUF_FORMAT) != 0) {
        PyBuffer_Release(&mu_i_buf);
        PyBuffer_Release(&mu_e_buf);
        PyBuffer_Release(&nu_i_buf);
        return NULL;
    }

    non(
        (const double *)mu_i_buf.buf,
        (const double *)mu_e_buf.buf,
        (size_t)nof_faces,
        (double *)nu_i_buf.buf,
        (double *)nu_e_buf.buf
    );

    PyBuffer_Release(&mu_i_buf);
    PyBuffer_Release(&mu_e_buf);
    PyBuffer_Release(&nu_i_buf);
    PyBuffer_Release(&nu_e_buf);

    Py_RETURN_NONE;
}

static PyObject *py_nu(PyObject *self, PyObject *args) {
    PyObject *faces = NULL;
    Py_ssize_t nof_faces = 0;

    PyObject *nodes = NULL;
    Py_ssize_t nof_nodes = 0;

    PyObject *normals = NULL;
    PyObject *centers = NULL;
    PyObject *s = NULL;
    PyObject *nu_i = NULL;

    if (!PyArg_ParseTuple(args, "OnOnOOOO", &faces, &nof_faces, &nodes, &nof_nodes, &normals, &centers, &s, &nu_i)) {
        return NULL;
    }

    Py_buffer faces_buf;
    Py_buffer nodes_buf;
    Py_buffer normals_buf;
    Py_buffer centers_buf;
    Py_buffer s_buf;
    Py_buffer nu_i_buf;

    if (PyObject_GetBuffer(faces, &faces_buf, PyBUF_FORMAT) != 0) return NULL;
    if (PyObject_GetBuffer(nodes, &nodes_buf, PyBUF_FORMAT) != 0) {
        PyBuffer_Release(&faces_buf);
        return NULL;
    }
    if (PyObject_GetBuffer(normals, &normals_buf, PyBUF_FORMAT) != 0) {
        PyBuffer_Release(&faces_buf);
        PyBuffer_Release(&nodes_buf);
        return NULL;
    }
    if (PyObject_GetBuffer(centers, &centers_buf, PyBUF_FORMAT) != 0) {
        PyBuffer_Release(&faces_buf);
        PyBuffer_Release(&nodes_buf);
        PyBuffer_Release(&normals_buf);
        return NULL;
    }
    if (PyObject_GetBuffer(s, &s_buf, PyBUF_FORMAT) != 0) {
        PyBuffer_Release(&faces_buf);
        PyBuffer_Release(&nodes_buf);
        PyBuffer_Release(&normals_buf);
        PyBuffer_Release(&centers_buf);
        return NULL;
    }
    if (PyObject_GetBuffer(nu_i, &nu_i_buf, PyBUF_WRITABLE | PyBUF_FORMAT) != 0) {
        PyBuffer_Release(&faces_buf);
        PyBuffer_Release(&nodes_buf);
        PyBuffer_Release(&normals_buf);
        PyBuffer_Release(&centers_buf);
        PyBuffer_Release(&s_buf);
        return NULL;
    }

    nu(
        (const int (*)[3])faces_buf.buf, (size_t)nof_faces,
        (const double (*)[3])nodes_buf.buf, (size_t)nof_nodes,
        (const double (*)[3])normals_buf.buf,
        (const double (*)[3])centers_buf.buf,
        (const double *)s_buf.buf,
        (double *)nu_i_buf.buf
    );

    PyBuffer_Release(&faces_buf);
    PyBuffer_Release(&nodes_buf);
    PyBuffer_Release(&normals_buf);
    PyBuffer_Release(&centers_buf);
    PyBuffer_Release(&s_buf);
    PyBuffer_Release(&nu_i_buf);

    Py_RETURN_NONE;
}

static PyMethodDef methods[] = {
    {"dotProduct3D",      py_dotProduct3D,      METH_VARARGS, NULL},
    {"crossProduct3D",    py_crossProduct3D,    METH_VARARGS, NULL},
    {"normalizeVector3D", py_normalizeVector3D, METH_VARARGS, NULL},

    {"intersect_AB_t",    py_intersect_AB_t,    METH_VARARGS, NULL},

    {"boundingBox",       py_boundingBox,       METH_VARARGS, NULL},

    {"mu",                py_mu,                METH_VARARGS, NULL},
    {"non",               py_non,               METH_VARARGS, NULL},
    {"nu",                py_nu,                METH_VARARGS, NULL},

    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef module = {
    PyModuleDef_HEAD_INIT,
    "_tvet",
    NULL,
    -1,
    methods
};

PyMODINIT_FUNC PyInit__tvet(void) {
    return PyModule_Create(&module);
}

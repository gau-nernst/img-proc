// https://docs.python.org/3/extending/extending.html
#include "img_proc.h"
#include <stdint.h>

#define PY_SSIZE_T_CLEAN
#include <Python.h>

static PyObject *py_image_resize(PyObject *self, PyObject *args) {
  const uint8_t *image;
  Py_ssize_t image_size; // unused
  int width, height, depth, new_width, new_height;
  Interpolation interpolation;
  Py_buffer output;

  // y# -> read-only bytes-like object
  // i  -> int
  // y* -> bytes-like object
  if (!PyArg_ParseTuple(args, "y#iiiiiiy*", &image, &image_size, &width, &height, &depth, &new_width, &new_height,
                        &interpolation, &output))
    return NULL;

  Py_BEGIN_ALLOW_THREADS;
  image_resize(image, width, height, depth, new_width, new_height, interpolation, output.buf);
  Py_END_ALLOW_THREADS;

  PyBuffer_Release(&output);

  Py_INCREF(Py_None);
  return Py_None;
}

// method table
static PyMethodDef ImgProcMethods[] = {
    {"image_resize", py_image_resize, METH_VARARGS, ""}, //
    {NULL, NULL, 0, NULL},                               // sentinel
};

// module definition
static struct PyModuleDef img_proc_module = {
    PyModuleDef_HEAD_INIT, "img_proc", "module docstring", -1, ImgProcMethods,
};

// module initialization
PyMODINIT_FUNC PyInit_img_proc(void) { return PyModule_Create(&img_proc_module); }

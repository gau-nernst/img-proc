// This file is auto-generated from generate_binding.py
#include "img_proc.h"
#define PY_SSIZE_T_CLEAN
#include <Python.h>

static PyObject *py_image_resize(PyObject *self, PyObject *args) {
  const char * image;
  Py_ssize_t image_size;
  int width;
  int height;
  int depth;
  int new_width;
  int new_height;
  int interpolation;
  Py_buffer output;

  if (!PyArg_ParseTuple(args, "y#iiiiiiy*", &image, &image_size, &width, &height, &depth, &new_width, &new_height, &interpolation, &output))
    return NULL;

  Py_BEGIN_ALLOW_THREADS;
  image_resize(image, width, height, depth, new_width, new_height, interpolation, output.buf);
  Py_END_ALLOW_THREADS;

  PyBuffer_Release(&output);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *py_get_rotation_matrix_2d(PyObject *self, PyObject *args) {
  float cx;
  float cy;
  float angle;
  float scale;
  Py_buffer output;

  if (!PyArg_ParseTuple(args, "ffffy*", &cx, &cy, &angle, &scale, &output))
    return NULL;

  Py_BEGIN_ALLOW_THREADS;
  get_rotation_matrix_2d(cx, cy, angle, scale, output.buf);
  Py_END_ALLOW_THREADS;

  PyBuffer_Release(&output);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *py_invert_affine_transform(PyObject *self, PyObject *args) {
  const char * m;
  Py_ssize_t m_size;
  Py_buffer m_inv;

  if (!PyArg_ParseTuple(args, "y#y*", &m, &m_size, &m_inv))
    return NULL;

  Py_BEGIN_ALLOW_THREADS;
  invert_affine_transform(m, m_inv.buf);
  Py_END_ALLOW_THREADS;

  PyBuffer_Release(&m_inv);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *py_image_warp_affine(PyObject *self, PyObject *args) {
  const char * src;
  Py_ssize_t src_size;
  int width;
  int height;
  int depth;
  const char * transform;
  Py_ssize_t transform_size;
  int new_width;
  int new_height;
  int interpolation;
  Py_buffer dst;

  if (!PyArg_ParseTuple(args, "y#iiiy#iiiy*", &src, &src_size, &width, &height, &depth, &transform, &transform_size, &new_width, &new_height, &interpolation, &dst))
    return NULL;

  Py_BEGIN_ALLOW_THREADS;
  image_warp_affine(src, width, height, depth, transform, new_width, new_height, interpolation, dst.buf);
  Py_END_ALLOW_THREADS;

  PyBuffer_Release(&dst);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *py_invert_matrix_3x3(PyObject *self, PyObject *args) {
  const char * m;
  Py_ssize_t m_size;
  Py_buffer m_inv;

  if (!PyArg_ParseTuple(args, "y#y*", &m, &m_size, &m_inv))
    return NULL;

  Py_BEGIN_ALLOW_THREADS;
  invert_matrix_3x3(m, m_inv.buf);
  Py_END_ALLOW_THREADS;

  PyBuffer_Release(&m_inv);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *py_image_warp_perspective(PyObject *self, PyObject *args) {
  const char * src;
  Py_ssize_t src_size;
  int width;
  int height;
  int depth;
  const char * transform;
  Py_ssize_t transform_size;
  int new_width;
  int new_height;
  int interpolation;
  Py_buffer dst;

  if (!PyArg_ParseTuple(args, "y#iiiy#iiiy*", &src, &src_size, &width, &height, &depth, &transform, &transform_size, &new_width, &new_height, &interpolation, &dst))
    return NULL;

  Py_BEGIN_ALLOW_THREADS;
  image_warp_perspective(src, width, height, depth, transform, new_width, new_height, interpolation, dst.buf);
  Py_END_ALLOW_THREADS;

  PyBuffer_Release(&dst);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *py_image_box_filter(PyObject *self, PyObject *args) {
  const char * image;
  Py_ssize_t image_size;
  int width;
  int height;
  int depth;
  int kw;
  int kh;
  Py_buffer output;
  int impl;

  if (!PyArg_ParseTuple(args, "y#iiiiiy*i", &image, &image_size, &width, &height, &depth, &kw, &kh, &output, &impl))
    return NULL;

  Py_BEGIN_ALLOW_THREADS;
  image_box_filter(image, width, height, depth, kw, kh, output.buf, impl);
  Py_END_ALLOW_THREADS;

  PyBuffer_Release(&output);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *py_image_gaussian_filter(PyObject *self, PyObject *args) {
  const char * image;
  Py_ssize_t image_size;
  int width;
  int height;
  int channels;
  int kx;
  int ky;
  float sigma_x;
  float sigma_y;
  Py_buffer output;

  if (!PyArg_ParseTuple(args, "y#iiiiiffy*", &image, &image_size, &width, &height, &channels, &kx, &ky, &sigma_x, &sigma_y, &output))
    return NULL;

  Py_BEGIN_ALLOW_THREADS;
  image_gaussian_filter(image, width, height, channels, kx, ky, sigma_x, sigma_y, output.buf);
  Py_END_ALLOW_THREADS;

  PyBuffer_Release(&output);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyMethodDef ImgProcMethods[] = {
  {"image_resize", py_image_resize, METH_VARARGS, ""},
  {"get_rotation_matrix_2d", py_get_rotation_matrix_2d, METH_VARARGS, ""},
  {"invert_affine_transform", py_invert_affine_transform, METH_VARARGS, ""},
  {"image_warp_affine", py_image_warp_affine, METH_VARARGS, ""},
  {"invert_matrix_3x3", py_invert_matrix_3x3, METH_VARARGS, ""},
  {"image_warp_perspective", py_image_warp_perspective, METH_VARARGS, ""},
  {"image_box_filter", py_image_box_filter, METH_VARARGS, ""},
  {"image_gaussian_filter", py_image_gaussian_filter, METH_VARARGS, ""},
  {NULL, NULL, 0, NULL},
};

static struct PyModuleDef img_proc_module = {
    PyModuleDef_HEAD_INIT, "img_proc", "module docstring", -1, ImgProcMethods,
};

PyMODINIT_FUNC PyInit_img_proc(void) { return PyModule_Create(&img_proc_module); }

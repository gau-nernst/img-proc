#include "img_proc.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#define CLAMP(x, low, high) MIN(MAX((x), (low)), (high))

// there are different ways to handle out-of-bound values
static uint8_t image_value(const uint8_t *image, int width, int height, int depth, int col, int row, int d) {
  // return zero
  // return 0 <= col && col < width && 0 <= row && row < height ? image[(row * width + col) * depth + d] : 0;
  // clamp to border
  return image[(CLAMP(row, 0, height - 1) * width + CLAMP(col, 0, width - 1)) * depth + d];
}

static int round_half_down(double x) { return (int)ceil(x - 0.5); }
static int round_half_up(double x) { return (int)floor(x + 0.5); }

static void image_interpolate_nearest(const uint8_t *image, int width, int height, int depth, double x, double y,
                                      uint8_t *output) {
  int col = round_half_down(x - 0.5);
  int row = round_half_down(y - 0.5);
  for (int d = 0; d < depth; d++)
    output[d] = image_value(image, width, height, depth, col, row, d);
}

static void image_interpolate_bilinear(const uint8_t *image, int width, int height, int depth, double x, double y,
                                       uint8_t *output) {
  x -= 0.5;
  y -= 0.5;
  int col1 = (int)floor(x);
  int row1 = (int)floor(y);
  int col2 = col1 + 1;
  int row2 = row1 + 1;

  double dx = x - floor(x);
  double dy = y - floor(y);

  for (int d = 0; d < depth; d++) {
    double p1 = (double)image_value(image, width, height, depth, col1, row1, d);
    double p2 = (double)image_value(image, width, height, depth, col2, row1, d);
    double p3 = (double)image_value(image, width, height, depth, col1, row2, d);
    double p4 = (double)image_value(image, width, height, depth, col2, row2, d);

    double p = p1 * (1.0 - dx) * (1.0 - dy) + p2 * dx * (1.0 - dy) + p3 * (1.0 - dx) * dy + p4 * dx * dy;
    output[d] = (uint8_t)round(p);
  }
}

static double filter_bicubic(double x0, double x1, double x2, double x3, double t) {
  return x1 + (-0.5 * x0 + 0.5 * x2) * t + (x0 - 2.5 * x1 + 2.0 * x2 - 0.5 * x3) * t * t +
         (-0.5 * x0 + 1.5 * x1 - 1.5 * x2 + 0.5 * x3) * t * t * t;
}

// bicubic spline: derivatives at the corners/boundaries are maintained
static void image_interpolate_bicubic(const uint8_t *image, int width, int height, int depth, double x, double y,
                                      uint8_t *output) {
  x -= 0.5;
  y -= 0.5;
  int cols[4], rows[4];
  for (int i = 0; i < 4; i++) {
    cols[i] = (int)floor(x) + i - 1;
    rows[i] = (int)floor(y) + i - 1;
  }

  double dx = x - floor(x);
  double dy = y - floor(y);

  double ys[4];
  for (int d = 0; d < depth; d++) {
    for (int j = 0; j < 4; j++) {
      double x0 = (double)image_value(image, width, height, depth, cols[0], rows[j], d);
      double x1 = (double)image_value(image, width, height, depth, cols[1], rows[j], d);
      double x2 = (double)image_value(image, width, height, depth, cols[2], rows[j], d);
      double x3 = (double)image_value(image, width, height, depth, cols[3], rows[j], d);
      ys[j] = filter_bicubic(x0, x1, x2, x3, dx);
    }
    double out = filter_bicubic(ys[0], ys[1], ys[2], ys[3], dy);
    output[d] = (uint8_t)CLAMP(out, 0.0, 255.0);
  }
}

static void image_interpolate(const uint8_t *image, int width, int height, int depth, double x, double y,
                              uint8_t *output, Interpolation interpolation) {
  switch (interpolation) {
  case NEAREST:
    return image_interpolate_nearest(image, width, height, depth, x, y, output);
  case BILINEAR:
    return image_interpolate_bilinear(image, width, height, depth, x, y, output);
  case BICUBIC:
    return image_interpolate_bicubic(image, width, height, depth, x, y, output);
  }
}

void image_resize(const uint8_t *image, int width, int height, int depth, int new_width, int new_height,
                  Interpolation interpolation, uint8_t *output) {
  double x_scale = (double)width / (double)new_width;
  double y_scale = (double)height / (double)new_height;

  for (int row = 0; row < new_height; row++)
    for (int col = 0; col < new_width; col++) {
      double x = ((double)col + 0.5) * x_scale;
      double y = ((double)row + 0.5) * y_scale;
      uint8_t *pixel_output = output + (row * new_width + col) * depth;
      image_interpolate(image, width, height, depth, x, y, pixel_output, interpolation);
    }
}

void get_rotation_matrix_2d(double cx, double cy, double angle, double scale, double *output) {
  angle *= M_PI / 180.0;
  double cos_a = scale * cos(angle);
  double sin_a = scale * sin(angle);

  output[0] = cos_a;
  output[1] = sin_a;
  output[2] = (1.0 - cos_a) * cx - sin_a * cy;
  output[3] = -sin_a;
  output[4] = cos_a;
  output[5] = sin_a * cx + (1.0 - cos_a) * cy;
}

void invert_affine_transform(const double *m, double *m_inv) {
  double m00 = m[0];
  double m01 = m[1];
  double m02 = m[2];
  double m10 = m[3];
  double m11 = m[4];
  double m12 = m[5];

  double inv_det = 1.0 / (m00 * m11 - m01 * m10);
  m_inv[0] = m11 * inv_det;
  m_inv[1] = -m01 * inv_det;
  m_inv[3] = -m10 * inv_det;
  m_inv[4] = m00 * inv_det;
  m_inv[2] = -(m_inv[0] * m02 + m_inv[1] * m12);
  m_inv[5] = -(m_inv[3] * m02 + m_inv[4] * m12);
}

void image_warp_affine(const uint8_t *src, int width, int height, int depth, const double *transform, int new_width,
                       int new_height, Interpolation interpolation, uint8_t *dst) {
  double inv_transform[6];
  invert_affine_transform(transform, inv_transform);

  double m00 = inv_transform[0];
  double m01 = inv_transform[1];
  double m02 = inv_transform[2];
  double m10 = inv_transform[3];
  double m11 = inv_transform[4];
  double m12 = inv_transform[5];

  for (int dst_y = 0; dst_y < new_height; dst_y++) {
    for (int dst_x = 0; dst_x < new_width; dst_x++) {
      double src_x = m00 * ((double)dst_x + 0.5) + m01 * ((double)dst_y + 0.5) + m02;
      double src_y = m10 * ((double)dst_x + 0.5) + m11 * ((double)dst_y + 0.5) + m12;

      if (src_x < 0.0 || src_x > (double)width || src_y < 0.0 || src_y > (double)height)
        continue;

      uint8_t *pixel_dst = dst + (dst_y * new_width + dst_x) * depth;
      image_interpolate(src, width, height, depth, src_x, src_y, pixel_dst, interpolation);
    }
  }
}

void invert_matrix_3x3(const double *m, double *m_inv) {
  double m00 = m[0];
  double m01 = m[1];
  double m02 = m[2];
  double m10 = m[3];
  double m11 = m[4];
  double m12 = m[5];
  double m20 = m[6];
  double m21 = m[7];
  double m22 = m[8];

  double det = m00 * (m11 * m22 - m12 * m21)   //
               - m01 * (m10 * m22 - m12 * m20) //
               + m02 * (m10 * m21 - m11 * m20);
  double inv_det = 1.0 / det;

  m_inv[0] = (m11 * m22 - m12 * m21) * inv_det;
  m_inv[1] = -(m01 * m22 - m02 * m21) * inv_det;
  m_inv[2] = (m01 * m12 - m02 * m11) * inv_det;
  m_inv[3] = -(m10 * m22 - m12 * m20) * inv_det;
  m_inv[4] = (m00 * m22 - m02 * m20) * inv_det;
  m_inv[5] = -(m00 * m12 - m02 * m10) * inv_det;
  m_inv[6] = (m10 * m21 - m11 * m20) * inv_det;
  m_inv[7] = -(m00 * m21 - m01 * m20) * inv_det;
  m_inv[8] = (m00 * m11 - m01 * m10) * inv_det;
}

void image_warp_perspective(const uint8_t *src, int width, int height, int depth, const double *transform,
                            int new_width, int new_height, Interpolation interpolation, uint8_t *dst) {
  double inv_transform[9];
  invert_matrix_3x3(transform, inv_transform);

  double m00 = inv_transform[0];
  double m01 = inv_transform[1];
  double m02 = inv_transform[2];
  double m10 = inv_transform[3];
  double m11 = inv_transform[4];
  double m12 = inv_transform[5];
  double m20 = inv_transform[6];
  double m21 = inv_transform[7];
  double m22 = inv_transform[8];

  for (int dst_y = 0; dst_y < new_height; dst_y++) {
    for (int dst_x = 0; dst_x < new_width; dst_x++) {
      double scale = 1.0 / (m20 * ((double)dst_x + 0.5) + m21 * ((double)dst_y + 0.5) + m22);
      double src_x = scale * (m00 * ((double)dst_x + 0.5) + m01 * ((double)dst_y + 0.5) + m02);
      double src_y = scale * (m10 * ((double)dst_x + 0.5) + m11 * ((double)dst_y + 0.5) + m12);

      if (src_x < 0.0 || src_x > (double)width || src_y < 0.0 || src_y > (double)height)
        continue;

      uint8_t *pixel_dst = dst + (dst_y * new_width + dst_x) * depth;
      image_interpolate(src, width, height, depth, src_x, src_y, pixel_dst, interpolation);
    }
  }
}

void image_box_filter(const uint8_t *image, int width, int height, int depth, int kw, int kh, uint8_t *output) {
  int rw = kw / 2;
  int rh = kh / 2;

  // implementations:
  // - naive                             O(width x height x kw x kh) -> no extra memory
  // - integral image                    O(width x height)           -> will overflow, extra memory (height x width x depth)
  // - separable conv                    O(width x kw + height x kh) -> extra memory (kh x width x depth)
  // - separable + online moving average O(width x height)           -> extra memory (kh x width x depth)
  for (int dst_row = 0; dst_row < height; dst_row++) {
    for (int dst_col = 0; dst_col < width; dst_col++) {
      for (int d = 0; d < depth; d++) {
        // border handling: normalize over visible area only
        int value = 0;
        int count = 0;

        for (int src_row = MAX(dst_row - rh, 0); src_row < MIN(dst_row + rh + 1, height - 1); src_row++) {
          for (int src_col = MAX(dst_col - rw, 0); src_col < MIN(dst_col + rw + 1, width - 1); src_col++) {
            count += 1;
            value += image[(src_row * width + src_col) * depth + d];
          }
        }

        double value_f64 = (double)value / (double)count;
        output[(dst_row * width + dst_col) * depth + d] = (uint8_t)value_f64;
      }
    }
  }
}

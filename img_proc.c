#include "img_proc.h"
#include <math.h>

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

void invert_affine_transform(const double *transform, double *output) {
  double a = transform[0];
  double b = transform[1];
  double c = transform[2];
  double d = transform[3];
  double e = transform[4];
  double f = transform[5];

  double inv_det = 1.0 / (a * e - b * d);
  output[0] = e * inv_det;
  output[1] = -b * inv_det;
  output[3] = -d * inv_det;
  output[4] = a * inv_det;
  output[2] = -(output[0] * c + output[1] * f);
  output[5] = -(output[3] * c + output[4] * f);
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

void image_warp_affine(const uint8_t *src, int width, int height, int depth, const double *transform, int new_width,
                       int new_height, Interpolation interpolation, uint8_t *dst) {
  double inv_transform[6];
  invert_affine_transform(transform, inv_transform);

  double a = inv_transform[0];
  double b = inv_transform[1];
  double c = inv_transform[2];
  double d = inv_transform[3];
  double e = inv_transform[4];
  double f = inv_transform[5];

  for (int dst_y = 0; dst_y < new_height; dst_y++) {
    for (int dst_x = 0; dst_x < new_width; dst_x++) {
      double src_x = a * ((double)dst_x + 0.5) + b * ((double)dst_y + 0.5) + c;
      double src_y = d * ((double)dst_x + 0.5) + e * ((double)dst_y + 0.5) + f;

      if (src_x < 0.0 || src_x > (double)width || src_y < 0.0 || src_y > (double)height)
        continue;

      uint8_t *pixel_dst = dst + (dst_y * new_width + dst_x) * depth;
      image_interpolate(src, width, height, depth, src_x, src_y, pixel_dst, interpolation);
    }
  }
}

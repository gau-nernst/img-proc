#include "img_proc.h"
#include <math.h>

#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#define CLAMP(x, low, high) MIN(MAX((x), (low)), (high))

static void interpolate_nearest(const uint8_t *image, int width, int height, int depth, double x, double y,
                                uint8_t *output) {
  int col = CLAMP((int)x, 0, width - 1);
  int row = CLAMP((int)y, 0, height - 1);
  for (int i = 0; i < depth; i++)
    output[i] = image[(row * width + col) * depth + i];
}

static void interpolate_bilinear(const uint8_t *image, int width, int height, int depth, double x, double y,
                                 uint8_t *output) {
  x -= 0.5;
  y -= 0.5;
  int col1 = CLAMP((int)floor(x), 0, width - 1);
  int row1 = CLAMP((int)floor(y), 0, height - 1);
  int col2 = CLAMP((int)ceil(x), 0, width - 1);
  int row2 = CLAMP((int)ceil(y), 0, height - 1);

  double dx = CLAMP(x, 0.0, (double)(width - 1)) - (double)col1;
  double dy = CLAMP(y, 0.0, (double)(height - 1)) - (double)row1;

  for (int i = 0; i < depth; i++) {
    double p1 = (double)image[(row1 * width + col1) * depth + i];
    double p2 = (double)image[(row1 * width + col2) * depth + i];
    double p3 = (double)image[(row2 * width + col1) * depth + i];
    double p4 = (double)image[(row2 * width + col2) * depth + i];

    double p = p1 * (1.0 - dx) * (1.0 - dy) + p2 * dx * (1.0 - dy) + p3 * (1.0 - dx) * dy + p4 * dx * dy;
    output[i] = (uint8_t)round(p);
  }
}

static void interpolate_bicubic(const uint8_t *image, int width, int height, int depth, double x, double y,
                                uint8_t *output) {
  // TODO
}

static void interpolate(const uint8_t *image, int width, int height, int depth, double x, double y, uint8_t *output,
                        Interpolation interpolation) {
  switch (interpolation) {
  case NEAREST:
    return interpolate_nearest(image, width, height, depth, x, y, output);
  case BILINEAR:
    return interpolate_bilinear(image, width, height, depth, x, y, output);
  case BICUBIC:
    return interpolate_bicubic(image, width, height, depth, x, y, output);
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
      interpolate(image, width, height, depth, x, y, pixel_output, interpolation);
    }
}

#include <stdint.h>

typedef enum Interpolation {
  NEAREST,
  BILINEAR,
  BICUBIC,
} Interpolation;

void image_resize(const uint8_t *image, int width, int height, int depth, int new_width, int new_height,
                  Interpolation interpolation, uint8_t *output);
void get_rotation_matrix_2d(double cx, double cy, double angle, double scale, double *output);
void invert_affine_transform(const double *m, double *m_inv);
void image_warp_affine(const uint8_t *src, int width, int height, int depth, const double *transform, int new_width,
                       int new_height, Interpolation interpolation, uint8_t *dst);
void invert_matrix_3x3(const double *m, double *m_inv);
void image_warp_perspective(const uint8_t *src, int width, int height, int depth, const double *transform,
                            int new_width, int new_height, Interpolation interpolation, uint8_t *dst);

#include <stdint.h>

typedef enum Interpolation {
  NEAREST,
  BILINEAR,
  BICUBIC,
} Interpolation;

void image_resize(const uint8_t *image, int width, int height, int channels, int new_width, int new_height,
                  Interpolation interpolation, uint8_t *output);
void get_rotation_matrix_2d(float cx, float cy, float angle, float scale, float *output);
void invert_affine_transform(const float *m, float *m_inv);
void image_warp_affine(const uint8_t *src, int width, int height, int channels, const float *transform, int new_width,
                       int new_height, Interpolation interpolation, uint8_t *dst);
void invert_matrix_3x3(const float *m, float *m_inv);
void image_warp_perspective(const uint8_t *src, int width, int height, int channels, const float *transform,
                            int new_width, int new_height, Interpolation interpolation, uint8_t *dst);

// only supports odd kernel size
void image_box_filter(const uint8_t *image, int width, int height, int channels, int kx, int ky, uint8_t *output,
                      int impl);
void image_gaussian_filter(const uint8_t *image, int width, int height, int channels, int kx, int ky, float sigma_x,
                           float sigma_y, uint8_t *output);

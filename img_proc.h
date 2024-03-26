#include <stdint.h>

typedef enum Interpolation {
  NEAREST,
  BILINEAR,
  BICUBIC,
} Interpolation;

void image_resize(const uint8_t *image, int width, int height, int depth, int new_width, int new_height,
                  Interpolation interpolation, uint8_t *output);

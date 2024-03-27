#include "img_proc.h"
#include <stdio.h>
#include <string.h>

#define CHECK(cond)                                                                                                    \
  if (!(cond)) {                                                                                                       \
    printf("Check fails %s:%d\n", __FILE__, __LINE__);                                                                 \
  }

void test_interpolate_nearest() {
  const uint8_t image[] = {1, 2, 3, //
                           4, 5, 6};
  uint8_t output[6 * 4];
  uint8_t *expected;

  image_resize(image, 3, 2, 1, 3, 2, NEAREST, output);
  CHECK(memcmp(output, image, 6) == 0);

  image_resize(image, 3, 2, 1, 6, 4, NEAREST, output);
  expected = (uint8_t[]){1, 1, 2, 2, 3, 3, //
                         1, 1, 2, 2, 3, 3, //
                         4, 4, 5, 5, 6, 6, //
                         4, 4, 5, 5, 6, 6};
  CHECK(memcmp(output, expected, 6 * 4) == 0);

  image_resize(image, 3, 2, 1, 1, 1, NEAREST, output);
  CHECK(output[0] == 2);
}

void test_interpolate_bilinear() {
  const uint8_t image[] = {1, 2, 3, //
                           4, 5, 6};
  uint8_t output[6 * 4];
  uint8_t *expected;

  image_resize(image, 3, 2, 1, 3, 2, BILINEAR, output);
  CHECK(memcmp(output, image, 6) == 0);

  // match OpenCV. Pillow makes more sense though.
  image_resize(image, 3, 2, 1, 6, 4, BILINEAR, output);
  expected = (uint8_t[]){1, 1, 2, 2, 3, 3, //
                         2, 2, 3, 3, 4, 4, //
                         3, 4, 4, 5, 5, 5, //
                         4, 4, 5, 5, 6, 6};
  CHECK(memcmp(output, expected, 6 * 4) == 0);

  image_resize(image, 3, 2, 1, 1, 1, BILINEAR, output);
  CHECK(output[0] == 4); // (2 + 5) / 2
}

int main() {
  test_interpolate_nearest();
  test_interpolate_bilinear();

  return 0;
}

# Image processing

[![Build](https://github.com/gau-nernst/img-proc/actions/workflows/build.yaml/badge.svg)](https://github.com/gau-nernst/img-proc/actions/workflows/build.yaml)

TODO:

- [ ] Image resampling:
    - [x] Interpolation: nearest ✅, bilinear ✅, bicubic ✅
    - [ ] Convolution-based (like pillow)
- [ ] Geometric transform:
    - [x] Affine transform ✅, Perspective transform ✅
    - [x] Handle border correctly to not cause jagged edges ✅
    - [ ] Estimate affine transform: full, partial (scale, rotation, and translation only)
    - [ ] Estimate perspective transform: exact from 4 pairs of points (`cv.getPerspectiveTransform`)
- [ ] Filters:
    - [ ] Basic: box, Gaussian, median
    - [ ] Edge preserving: bilateral, guided
- [ ] Image decoders:
    - [ ] JPEG (migrate from https://github.com/gau-nernst/jpeg.c)
    - [ ] PNG
    - [ ] TIFF

Resources:

- https://leimao.github.io/article/Interpolation/
- https://github.com/python-pillow/Pillow/blob/10.2.0/src/libImaging/Resample.c

## Python binding

Since I don't use custom structs, and don't use `malloc()` within the core library, creating binding for Python is super simple. In fact, it is so simple that I write a small script `generate_binding.py` to auto-generate binding from the header file `img_proc.h`. Numpy array works out-of-the-box thanks to [buffer protocol](https://docs.python.org/3/c-api/buffer.html). I also use Numpy to allocate and manage memory.

Note: there are no safety checks e.g. bounds check, correct dtype.

TODO:

- [ ] Auto-generate stubs `.pyi`

```bash
python generate_binding.py
python setup.py build_ext -i
```

```python
import numpy as np
from PIL import Image
import img_proc

img = Image.open("image.jpg")
img = np.array(img)
height, width, channels = img.shape
output = np.zeros((height, width, channels), dtype=np.uint8)

mat = np.empty((2, 3), dtype=np.float64)
img_proc.get_rotation_matrix_2d(width / 2, height / 2, 45.0, 1.0, mat)
img_proc.image_warp_affine(img, width, height, channels, mat, width, height, 2, output)

Image.fromarray(output)
```

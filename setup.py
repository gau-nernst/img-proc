# https://setuptools.pypa.io/en/latest/userguide/ext_modules.html
from setuptools import Extension, setup

SOURCES = ["img_proc.c", "py_img_proc.c"]


setup(
    ext_modules=[Extension(name="img_proc", sources=SOURCES)],
)

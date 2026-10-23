# Builds the C library (../c/md.c) into the package, where the wrapper loads
# it with ctypes. Install from this directory with `pip install .`
import os

from setuptools import Extension, setup

csrc = os.path.join("..", "c")

setup(
    ext_modules=[
        Extension(
            "minimum_degree._libmd",
            sources=[os.path.join("src", "minimum_degree", "_stub.c"), os.path.join(csrc, "md.c")],
            include_dirs=[csrc],
            define_macros=[("MD_BUILD_DLL", None)],
        )
    ]
)

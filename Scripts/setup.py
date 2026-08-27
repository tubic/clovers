"""
Setup script for building the Zcurve Cython extension module.

This module provides DNA sequence encoding functionality using the Z-curve method
for bioinformatics applications, particularly useful in gene prediction.
"""
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext
from Cython.Build import cythonize
import numpy as np


class BuildExt(_build_ext):
    """Apply compiler-specific flags (MSVC vs GCC/MinGW)."""

    def build_extensions(self):
        compiler = self.compiler.compiler_type  # 'msvc' | 'mingw32' | 'unix'
        for ext in self.extensions:
            if compiler == 'msvc':
                ext.extra_compile_args = ['/O2', '/openmp', '/arch:AVX2']
                ext.extra_link_args = []
            else:
                ext.extra_compile_args = ['-O3', '-mavx', '-mfma', '-fopenmp']
                ext.extra_link_args = ['-fopenmp']
        super().build_extensions()


# Extension module definition (Cython source, OpenMP-parallelized)
zcurve_extension = Extension(
    "Zcurve",
    sources=["Zcurve.pyx"],
    include_dirs=[".", np.get_include()],
)

setup(
    name="Zcurve",
    version="0.0.1",
    description="DNA sequence encoding using Z-curve method",
    long_description=__doc__,
    author="Zetong Zhang",
    author_email="zhangzetong@tju.edu.cn",
    url="https://github.com/zetong-zhang/clovers",
    license="GPLv3",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Cython",
    ],
    keywords="bioinformatics genomics z-curve gene-prediction dna-sequence",
    packages=[],
    ext_modules=cythonize(
        [zcurve_extension],
        compiler_directives={
            "boundscheck": False,
            "wraparound": False,
            "cdivision": True,
            "language_level": 3,
        },
    ),
    cmdclass={"build_ext": BuildExt},
    zip_safe=False,
    python_requires=">=3.6",
    install_requires=["numpy>=1.26.0", "cython>=0.29"],
)

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    #name = 'hik_kernel',
    ext_modules = [
        Extension('c_kernels', ['c_kernels.pyx'])
    ],
    cmdclass = {'build_ext': build_ext}
)

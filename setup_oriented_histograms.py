from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    name = 'calculate_histograms',
    ext_modules = [
        Extension('c_oriented_histograms', ['c_oriented_histograms.pyx'])
    ],
    cmdclass = {'build_ext': build_ext}
)

import numpy as np
import os

# WARNING: run python script like this!!!!!!!!!!
# python REALGenomeSIM_LD_cgetter_setup.py build_ext --inplace

os.chdir('/home/greggj/GxE/REALGenomeSIM')
import pyximport
pyximport.install()
import distutils.core
import Cython.Build
distutils.core.setup(
    ext_modules = Cython.Build.cythonize('REALGenomeSIM_LD_cgetter.pyx'),
    include_dirs = [np.get_include()])
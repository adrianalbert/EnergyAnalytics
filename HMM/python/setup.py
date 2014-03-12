from distutils.core import setup, Extension

module1 = Extension('HMM_C', sources = ['HMM_C.c'])

setup (name = 'HMM_C',
                 version = '1.0',
                 description = 'HMM filtering / smoothing utilities',
                 ext_modules = [module1])

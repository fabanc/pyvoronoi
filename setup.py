import sys
import os
from setuptools import setup
from setuptools.extension import Extension
from pathlib import Path

version = '1.2.4'

"""
Note on using the setup.py:
setup.py operates in 2 modes that are based on the presence of the 'dev' file in the root of the project.
 - When 'dev' is present, Cython will be used to compile the .pyx sources. This is the development mode
   (as you get it in the git repository).
 - When 'dev' is absent, C/C++ compiler will be used to compile the .cpp sources (that were prepared in
   in the development mode). This is the distribution mode (as you get it on PyPI).

This way the package can be used without or with an incompatible version of Cython.

The idea comes from: https://github.com/MattShannon/bandmat
"""
dev_mode = os.path.exists('dev')

if dev_mode:
    from Cython.Distutils import build_ext

    print('Development mode: Compiling Cython modules from .pyx sources.')
    sources = ["pyvoronoi/pyvoronoi.pyx", "pyvoronoi/voronoi.cpp"]

else:
    from setuptools.command.build_ext import build_ext

    print('Distribution mode: Compiling from Cython generated .cpp sources.')
    sources = ["pyvoronoi/pyvoronoi.cpp", "pyvoronoi/voronoi.cpp"]


ext = Extension("pyvoronoi",
                sources=sources,
                include_dirs = ["pyvoronoi"],
                language="c++",
                optional=os.environ.get('CIBUILDWHEEL', '0') != '1'
                )


# This command has been borrowed from
# http://www.pydanny.com/python-dot-py-tricks.html

if sys.argv[-1] == 'tag':
    tag_command = "git tag -a v%s -m 'v%s'" % (version, version)
    print(tag_command)
    os.system(tag_command)
    os.system("git push --tags")
    sys.exit()

# Run generation of pyi files
pyi_command = f'{os.path.dirname(sys.executable)}{os.path.sep}Scripts{os.path.sep}cythonpeg pyvoronoi/*.pyx'
print(pyi_command)
os.system(pyi_command)

class build_ext_subclass( build_ext ):
    def build_extensions(self):
        print(f'Compiler type: {self.compiler.compiler_type} - Version: {sys.version_info.major}.{sys.version_info.minor}')
        # Starting from 3.11, the file longintrepr.h has moved. It is no longer under Python@3.XX\include but Python@3.XX\include\cpython        
        if sys.version_info.major == 3 and sys.version_info.minor >= 11:
            c = self.compiler.compiler_type
            if c == 'msvc':
                for e in self.extensions:
                    install_dir = os.path.dirname(sys.executable)
                    cpython_directory = os.path.join(install_dir, 'include', 'cpython')         
                    e.extra_compile_args = [cpython_directory]
        build_ext.build_extensions(self)

      

this_directory = Path(__file__).parent
setup(
    name='pyvoronoi',
    python_requires='>=3.8',
    version=version,
    description='Cython wrapper for the Boost Voronoi library (version 1.59.0)',
    long_description=(this_directory / "README.md").read_text(),
    long_description_content_type='text/markdown',
    author='Fabien Ancelin / Andrii Sydorchuk, Voxel8',
    author_email='',
    url='https://github.com/fabanc/pyvoronoi',
    keywords=['voronoi','Boost','polygon'],
    data_files=['pyvoronoi/pyvoronoi.pyi'],
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Cython",
        "Programming Language :: C++",
        "Environment :: Other Environment",
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Operating System :: OS Independent",
        "License :: OSI Approved",
        "License :: OSI Approved :: MIT License",
        "Topic :: Multimedia :: Graphics",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: GIS",
        "Topic :: Software Development :: Libraries :: Python Modules"
    ],
    ext_modules=[ext],
    cmdclass={
        'build_ext': build_ext_subclass
    },   
)

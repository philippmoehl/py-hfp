from glob import glob
from pathlib import Path
from setuptools import setup, find_packages

from pybind11.setup_helpers import Pybind11Extension, build_ext

__version__ = "0.1"

source_path = Path(r"hfp/src").glob('**/*.cpp')


ext_modules = [
        Pybind11Extension(
            "pyhfp.hfp", 
            sources=[
                str(file.as_posix())
                for file in source_path
            ],
            include_dirs=["hfp/src", "hfp/include"],
            define_macros = [('VERSION_INFO', __version__)],
            ),
    ]

setup(
    name="pyhfp",
    version=__version__,
    author="Philipp Moehl",
    author_email="philipp1994@gmail.com",
    url="https://github.com/philippmoehl/py-hfp",
    description="Hierarchical Fitting Primitives.",
    long_description="",
    packages=find_packages(),
    ext_modules=ext_modules,
    install_requires=[
        'numpy>=1.20.0',
        'trimesh>=3.23.5'
    ],
    extras_require={"test": "pytest"},
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    python_requires=">=3.7",
)

from glob import glob
from pathlib import Path
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

source_path = Path(r"hfp/src").glob('**/*.cpp')


ext_modules = [
        Pybind11Extension(
            "hfp", 
            sources=[
                str(file.as_posix())
                for file in source_path
            ],
            include_dirs=["hfp/src", "hfp/include"],
            ),
    ]

setup(
    name="hfp",
    version="0.1.0",
    author="Philipp Moehl",
    author_email="philipp1994@gmail.com",
    description="Hierarchical Fitting Primitives.",
    long_description="",
    ext_modules=ext_modules,
    extras_require={"test": "pytest"},
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
)

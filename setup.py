import glob
import pybind11

from setuptools import setup, Extension
from pybind11.setup_helpers import Pybind11Extension, build_ext


ALL_SOURCE_FILES = [f for f in glob.glob("src/*.cpp") if not any(x in f for x in ["testing", "interface", "library"])]

ext_modules = [
    Extension(
        "Qubos.Extension",
        sources=ALL_SOURCE_FILES,
        include_dirs=[pybind11.get_include(), "src"],
        language="c++",
        extra_compile_args=["-std=c++20"],
    ),
]

setup(
    name="qubos",
    version="0.0.1",
    author="Șerban",
    description="A library for quantum Stabilizer simulation",
    packages=["qubos"],
    package_dir={"": "src"},
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    python_requires=">=3.8",
)

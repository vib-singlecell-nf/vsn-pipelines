
import glob
from setuptools import setup, find_packages
import os
import setuptools
with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='singlecelltxbenchmark',
    version='0.0.1',
    author="Maxime De Waegeneer",
    author_email="",
    description="Single-Cell Pipelines in Nextflow",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/aertslab/SingleCellTxBenchmark",
    packages=setuptools.find_packages(where='src'),
    package_dir={'': 'src'},
    py_modules=[os.path.splitext(os.path.basename(path))[0]
                for path in glob.glob('src/*.py')],
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
         "License :: OSI Approved :: MIT License",
         "Operating System :: OS Independent",
    ],
)
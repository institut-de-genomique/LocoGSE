#!/usr/bin/env python
# -*- coding: utf-8 -*-

from os import path

try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, "README.md")) as f:
    long_description = f.read()


setup(
    name="LocoGSE",
    version="0.0.1",
    license="CeCILL",
    description="\n\n A Genome Size Estimation program. It is based on a linear relation between the depth and the genome size. \n A correction which depends of the family is added at this depth for a best prediction. \n For a question : https://github.com/institut-de-genomique/LocoGSE/issues or pierre.guenzi.tiberi@gmail.com",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Pierre Guenzi-Tiberi",
    author_email="pierre.guenzi.tiberi@gmail.com",
    packages=[
        "LocoGSE",
        "LocoGSE.lib",
    ],
    package_dir={"LocoGSE": "LocoGSE"},
    package_data={
        "LocoGSE": [
            "slopes/PlantFamilies.CoeffRegression.V1.txt",
            "slopes/PlantFamilies.CoeffRegression.V2.txt",
            "slopes/PlantFamilies.CoeffRegression.BUSCO.V2.txt",
            "db/OneKP.410genes.consensus.dmnd",
            "db/OneKP.410genes.consensus.fa",
            "db/BUSCO.ancestral.dmnd",
            "db/BUSCO.ancestral.fa",
            # "lib/",
        ],
    },
    entry_points={
        "console_scripts": [
            "LocoGSE=LocoGSE.LocoGSE:run",
            "check_seq_names=LocoGSE.check_seq_names:run",
        ],
    },
    include_package_data=True,
    zip_safe=False,
    keywords="LocoGSE",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: CEA CNRS Inria Logiciel Libre License, version 2.1 (CeCILL-2.1)",
        "Programming Language :: Python :: 3.7.5",
    ],
)

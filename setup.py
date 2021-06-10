"""Setup script for cressp"""

import os.path
from setuptools import setup, find_packages

# The directory containing this file
HERE = os.path.abspath(os.path.dirname(__file__))

# The text of the README file
with open(os.path.join(HERE, "README.md")) as fid:
    README = fid.read()

setup(
    name='cressp',
    version='0.1.0',
    author="Hyunsu An",
    author_email="ahs2202@gm.gist.ac.kr",
    description="a program to find cross-reactive epitopes with structural information from known protein structures.",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/ahs2202/cressp",
    license="GPLv3",
    packages=find_packages( ),
    include_package_data=True,
    install_requires=[
        'tensorflow==2.3.0',
        'biobookshelf>=0.1.18',
        'mhcflurry>=2.0.2',
    ],
    entry_points={
        "console_scripts": [
            "cressp=cressp.__main__:cressp",
        ]
    },
)

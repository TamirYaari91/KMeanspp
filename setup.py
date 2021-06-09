import os
from setuptools import setup, find_packages, Extension

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...

setup(
    name="mykmeanssp",
    version="0.0.4",
    author="Andrew Carter",
    author_email="andrewjcarter@gmail.com",
    description=("An demonstration of how to create, document, and publish "
                 "to the cheese shop a5 pypi.org."),
    license="BSD",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: Implementation :: CPython",
    ],
    ext_modules=[
        Extension(
            'mykmeanssp',
            ['kmeans.c'],
        )
    ]
)

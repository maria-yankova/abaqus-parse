"""Pip installation script."""

import os
import re
from setuptools import find_packages, setup


def get_version():

    ver_file = 'abaqus_parse/_version.py'
    with open(ver_file) as handle:
        ver_str_line = handle.read()

    ver_pattern = r'^__version__ = [\'"]([^\'"]*)[\'"]'
    match = re.search(ver_pattern, ver_str_line, re.M)
    if match:
        ver_str = match.group(1)
    else:
        msg = 'Unable to find version string in "{}"'.format(ver_file)
        raise RuntimeError(msg)

    return ver_str


def get_long_description():
    readme_file = 'README.md'
    with open(readme_file, encoding='utf-8') as handle:
        contents = handle.read()
    return contents


def get_changelog():
    changelog_file = 'CHANGELOG.md'
    with open(changelog_file, encoding='utf-8') as handle:
        contents = handle.read()
    return contents


setup(
    author="Maria S. Yankova",
    author_email='maria.yankova@manchester.ac.uk',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Operating System :: OS Independent',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Abaqus input and output file readers and writers.",
    install_requires=[
        'numpy',
        'fracture_fea_laf==0.1.1',
        'plotly',
    ],
    license="MIT license",
    long_description=get_long_description() + '\n\n' + get_changelog(),
    long_description_content_type='text/markdown',
    keywords='abaqus, fea, materials-science',
    name='abaqus-parse',
    packages=find_packages(),
    project_urls={
        'GitHub': 'https://github.com/maria-yankova/abaqus-parse'
    },
    version=get_version(),
)
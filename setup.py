# Copyright (c) 2019 J. Alvarez-Jarreta and C.J. Brasher
#
# This file is part of the LipidFinder software tool and governed by the
# 'MIT License'. Please see the LICENSE file that should have been
# included as part of this software.

import setuptools


with open('README.md', 'r') as readme:
    long_description = readme.read()

__version__ = 'unknown'
for line in open('LipidFinder/__init__.py'):
    if (line.startswith(('__version__', '__author__', '__email__'))):
        exec(line.strip())

setuptools.setup(
    name = "LipidFinder",
    version = __version__,
    author = __author__,
    author_email = __email__,
    description = "", # TODO: add description
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = "https://github.com/ODonnell-Lipidomics/LipidFinder",
    packages = setuptools.find_packages(),
    keywords = 'lipidomics LC/MS profile',
    python_requires = '>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, <4',
    classifiers = [
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires = [
        'cycler',
        'ipython',
        'ipywidgets',
        'matplotlib',
        'numpy',
        'pandas',
        'requests-toolbelt',
        'scipy',
        'setuptools',
        'xlsxwriter',
    ],
    entry_points = {
        'console_scripts': [
            'config_params.py=LipidFinder.config_params:main',
            'run_peakfilter.py=LipidFinder.run_peakfilter:main',
            'run_amalgamator.py=LipidFinder.run_amalgamator:main',
            'run_mssearch.py=LipidFinder.run_mssearch:main',
            'update_params.py=LipidFinder.update_params:main',
            ],
    },
    package_data = {
        'LipidFinder': ['Data/*']
    },
    include_package_data = True,
    #test_suite = "tests",
)

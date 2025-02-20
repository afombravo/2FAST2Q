from setuptools import setup, find_packages

from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

VERSION = '2.7.6' 
DESCRIPTION = """A Python3 program that counts sequence occurences in fastq files."""

setup(
        name="fast2q", 
        version=VERSION,
        author="Afonso M Bravo",
        author_email="<afonsombravo@hotmail.com>",
        url='https://github.com/afombravo/2FAST2Q',
        description=DESCRIPTION,
        long_description=long_description,
        long_description_content_type='text/markdown',
        packages=find_packages(),
        install_requires=["numpy >= 1.19.2",\
                           "matplotlib >= 3.3.4",\
                           "numba >= 0.53.1",\
                           "psutil",\
                           "argparse",\
                           "datetime",\
                           "tqdm >= 4.59.0",\
                           "dataclasses",\
                           "tk-tools >= 0.1",
                           "colorama"],
        entry_points={
        'console_scripts': [
            '2fast2q=fast2q.fast2q:main',  # Replace `2fast2q` with the command name you want to use
            ],
        },
        keywords=['CRISPRi', 'CRISPRi-Seq', 'essentiality', "gene fitness","fastq","screening"],
        classifiers= [
            "Development Status :: 6 - Mature",
            "Intended Audience :: Science/Research",
            "Programming Language :: Python :: 3",
            "Operating System :: OS Independent",
        ],
        include_package_data=True,
        package_data={'': ['data/*']},
)

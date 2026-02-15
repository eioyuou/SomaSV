from setuptools import setup, find_packages

setup(
    name='SomaSV',
    version='0.0.1',
    author='Rt G',
    description='Long-read Somatic Structural Variant Hunter',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/eioyuou/SomaSV',
    packages=find_packages(),
    python_requires='>=3.8',
    install_requires=[
        'numpy>=1.21.0',
        'pysam>=0.19.0',
        'numba>=0.56.0',
        'scipy>=1.7.0',
        'intervaltree>=3.1.0',
    ],
    entry_points={
        'console_scripts': [
            'somasv=main:main',
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX :: Linux',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)
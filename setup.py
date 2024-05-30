from setuptools import setup

setup(
    name='hicpep',
    version='0.0.1',
    description="Hi-C Pearson matrix's Estimated-Pc1-pattern is a toolkit \
                used for Hi-C AB compartment analysis.",
    author='Zhi-Rong Cheng',
    license='MIT',
    packages = ['hicpep'],
    install_requires=[
        'hic-straw>=1.3.1',
        'matplotlib>=3.8.1',
        'numpy>=1.26.1',
        'pandas>=1.5.3'
    ],
    python_requires='>=3.11.5',
)
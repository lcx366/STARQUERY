from setuptools import setup,find_packages 

setup(
    name='starcatalogquery',
    version='2.0.1',
    description='A package for building an offline star catalog query database',
    author='Chunxiao Li',
    author_email='lcx366@126.com',
    url='https://github.com/lcx366/STARQUERY',
    license='MIT',
    long_description_content_type='text/markdown',
    long_description=open('README.md', 'rb').read().decode('utf-8'),
    keywords = ['StarCatalog'],
    python_requires = '>=3.10',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
        'License :: OSI Approved :: MIT License',
        ],
    packages = find_packages(),
    include_package_data=True,
    install_requires=[
        'numpy',
        'scipy',
        'astropy',
        'pandas',
        'h5py',
        'healpy',
        'cartopy',
        'wget',
        'tqdm',
        'natsort',
        'duckdb'
        ],
)

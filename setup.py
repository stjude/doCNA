import setuptools

setuptools.setup(
    name='doCNA',
    entry_points = {
        'console_scripts': ['docna=doCNA.__main__:main',],
    },
    version='0.8.2',
    author='Karol Szlachta, Dennis Kennetz',
    description='Analysis of copy number profile in WGS.',
    packages=['doCNA'],
    package_dir={'doCNA': 'doCNA'},
    package_data={'doCNA': ['doCNA.ini']},
    install_requires=[
        'numpy==1.17.1',
        'pandas==1.1.5',
        'scipy==1.3.1',
        'scikit-learn==0.20.0',
    ],
    python_requires='>=3.7.0',
    zip_safe=True,
    license='Apache2.0',
    url='https://github.com/stjude/doCNA',
)

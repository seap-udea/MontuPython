import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    # ######################################################################
    # BASIC DESCRIPTION
    # ######################################################################
    name='montu',
    author='Jorge I. Zuluaga',
    author_email='jorge.zuluaga@udea.edu.co',
    description='Montu Python: astronomical ephemerides for the ancient world',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://pypi.org/project/montu',
    keywords='astronomy egypt history',
    license='MIT',

    # ######################################################################
    # CLASSIFIER
    # ######################################################################
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        ],
    version='0.20.5',

    # ######################################################################
    # FILES
    # ######################################################################
    package_dir={'': '.'},
    packages=setuptools.find_packages(where='.'),
    
    # ######################################################################
    # TESTS
    # ######################################################################
    test_suite='pytest',
    tests_require=['pytest'],

    # ######################################################################
    # DEPENDENCIES
    # ######################################################################
    install_requires=['scipy','ipython','matplotlib','tqdm','numpy','ephem',
                      'pymeeus','regex','pandas','tabulate',
                      'pyplanets','requests',
                      'dash','dash_bootstrap_components'],

    extras_require={
        'test': ['pytest'],
    },

    # ######################################################################
    # OPTIONS
    # ######################################################################
    include_package_data=True,
    package_data={'': ['data/*.*', 'tests/*.*']},
    scripts=['montu/bin/imontu', 'montu/bin/montu-gui'],
)

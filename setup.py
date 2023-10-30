import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    # ######################################################################
    # BASIC DESCRIPTION
    # ######################################################################
    name='montu',
    author='Jorge Zuluaga, Tito Vivas',
    author_email='jorge.zuluaga@gmail.com',
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
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        ],
    version='0.9.1',

    # ######################################################################
    # FILES
    # ######################################################################
    package_dir={'': '.'},
    packages=setuptools.find_packages(where='.'),
    
    # ######################################################################
    # ENTRY POINTS
    # ######################################################################
    entry_points={
        'console_scripts': ['install=montu.install:main'],
    },

    # ######################################################################
    # TESTS
    # ######################################################################
    test_suite='nose.collector',
    tests_require=['nose'],

    # ######################################################################
    # DEPENDENCIES
    # ######################################################################
    install_requires=['scipy','ipython','matplotlib','tqdm','numpy','ephem','pymeeus','regex','pandas','tabulate',
                      'astroquery','plotly','numpy','spiceypy','pyplanets','astropy','astroquery'],

    # ######################################################################
    # OPTIONS
    # ######################################################################
    include_package_data=True,
    package_data={'': ['data/*.*', 'tests/*.*']},
)

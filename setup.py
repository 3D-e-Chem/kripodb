from setuptools import setup, find_packages

exec(open('kripodb/version.py').read())

setup(
    name='kripodb',
    description='Library to interact with Kripo fragment, fingerprint and similarity data files.',
    version=__version__,
    packages=find_packages(),
    url='https://github.com/3D-e-Chem/kripodb',
    author='Stefan Verhoeven',
    author_email='s.verhoeven@esciencecenter.nl',
    install_requires=['pyroaring', 'blosc', 'tables', 'pandas', 'connexion', 'requests', 'scipy', 'progressbar2', 'six'],
    setup_requires=['pytest-runner'],
    package_data={
      'kripodb.webservice': ['swagger.yaml'],
    },
    include_package_data=True,
    zip_safe=False,
    entry_points={
        'console_scripts': [
            'kripodb=kripodb.script:main',
        ],
    },
    license='Apache',
    classifiers=[
        'License :: OSI Approved :: Apache Software License',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Framework :: Flask',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Environment :: Console',
        'Environment :: Web Environment',
    ]
)

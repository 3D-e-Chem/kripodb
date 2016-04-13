from setuptools import setup, find_packages

exec(open('kripodb/version.py').read())

setup(
    name='kripodb',
    description='Kripo database',
    version=__version__,
    packages=find_packages(),
    url='https://github.com/3D-e-Chem/kripodb',
    author='Stefan Verhoeven',
    author_email='s.verhoeven@esciencecenter.nl',
    install_requires=['intbitset', 'tables', 'pandas', 'connexion', 'requests'],
    package_data={
      'kripodb.webservice': ['swagger.json'],
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
    ]
)

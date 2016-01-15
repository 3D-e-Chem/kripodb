from setuptools import setup, find_packages

exec(open('modifiedtanimoto/version.py').read())

setup(
    name='modifiedtanimoto',
    description='Compute bit vector distance with modified tanimoto',
    version=__version__,
    packages=find_packages(),
    url='https://github.com/3D-e-Chem/python-modified-tanimoto',
    author='Stefan Verhoeven',
    author_email='s.verhoeven@esciencecenter.nl',
    install_requires=['intbitset', 'tables', 'nose', 'coverage', 'mock'],
    entry_points={
        'console_scripts': [
            'modified_tanimoto=modifiedtanimoto.script:main',
        ],
    },
    license='Apache',
    classifiers=[
        'License :: OSI Approved :: Apache Software License'
    ]
)

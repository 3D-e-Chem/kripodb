from setuptools import setup, find_packages
from setuptools.extension import Extension

exec(open('modifiedtanimoto/version.py').read())

# Only use Cython if it is available, else just use the pre-generated files
try:
    from Cython.Distutils import build_ext
    source_ext = '.pyx'
    cmdclass = {'build_ext': build_ext}
except ImportError:
    # If missing can be created with 'cython modifiedtanimoto/algorithm.pyx'
    source_ext = '.c'
    cmdclass = {}

ext_modules = [Extension('modifiedtanimoto.algorithm',
                         ['modifiedtanimoto/algorithm' + source_ext])]


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
    ],
    cmdclass=cmdclass,
    ext_modules=ext_modules
)

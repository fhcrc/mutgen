from setuptools import setup, find_packages

from mutgen.version import __version__

setup(name='mutgen',
        version=__version__,
        url='http://github.com/fhcrc/mutgen',
        description='Generate mutations under a motif based mutation model',
        author="Christopher Small",
        entry_points={
            'console_scripts': [
                'mutgen = mutgen.cli:main'
            ]},
        packages=find_packages(exclude=['tests']),
        install_requires=['biopython'])


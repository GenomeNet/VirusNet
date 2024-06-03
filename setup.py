# setup.py

from setuptools import setup, find_packages

setup(
    name='virusnet',
    version='0.1.0',
    packages=find_packages(),
    scripts=['virusnet/cli.py'],
    entry_points={
        'console_scripts': [
            'virusnet=virusnet.cli:main',
        ],
    },
)
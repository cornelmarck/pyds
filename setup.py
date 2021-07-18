from setuptools import setup, find_packages, Command
from os import system

class CleanCommand(Command):
    """Custom clean command to tidy up the project root."""
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        system('rm -vrf ./build ./dist ./*.pyc ./*.tgz ./*.egg-info')


setup(
    name='pyds',
    version='1.0.0',
    author='Cornel Marck',
    author_email='cm1520@ic.ac.uk',

    packages=find_packages(),
    cmdclass={'clean': CleanCommand}
)
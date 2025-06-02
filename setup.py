from setuptools import setup
setup(
    name='SMol_FIESTA',
    version='1.0',
    description='',
    author='Reyes-Lamothe Lab, McGill University',
    author_email='jose.rasconperez@mail.mcgill.ca',
    packages=['SMol_FIESTA'],
    install_requires=[
        "setuptools>=68.2.0",
        "numpy~=2.2.3",
        "matplotlib~=3.10.0",
        "pandas~=2.2.3",
        "scikit-image~=0.25.2",
        "scipy~=1.15.2",
        "tifffile~=2025.2.18",
        "seaborn~=0.13.2",
        "natsort~=8.4.0",
        "tqdm~=4.67.1",
        "imagecodecs~=2024.12.30"
    ],
    entry_points={
            'console_scripts': ['SMF=SMol_FIESTA.BatchRun_Rebind:run_scripts']
    }
)
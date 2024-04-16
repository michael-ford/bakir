from setuptools import setup, find_packages

setup(
    name="bakir",
    version="0.1",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    package_data={
        'bakir': ['data/*'],
    },
    entry_points={
        'console_scripts': [
            'bakir=bakir.main:run',
        ],
    },
)

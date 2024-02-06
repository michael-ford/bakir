from setuptools import setup, find_packages

setup(
    name="kir_annotator",
    version="0.1",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    package_data={
        'kir_annotator': ['data/*'],
    },
    entry_points={
        'console_scripts': [
            'kir-annotator=kir_annotator.main:run',
        ],
    },
)

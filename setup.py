import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

    setuptools.setup(
    name="piperna",
    version="0.6",
    author="Scott Furlan",
    author_email="scottfurlan@gmail.com",
    description="A python wrapper for fast and parallel processing of bulk RNA Seq data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/scfurl/piperna.git",
    packages=setuptools.find_packages(),
    package_data={'piperna': ['piperna/data/genomes.json']},
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=2.5',
    entry_points={'console_scripts': [
        'piperna = piperna.__main__:run_piperna',
    ]},
    )


"""
cd ~/Box\ Sync/PI_FurlanS/computation/develop/piperna/
rm dist/*
python3 setup.py sdist bdist_wheel
python3 -m twine upload dist/*



**At SCRI do the following**

module load python
python3 -m pip install --user pipx
python3 -m pipx ensurepath
pipx install --include-deps --pip-args '--trusted-host pypi.org --trusted-host files.pythonhosted.org' piperna
pipx install --spec git+https://github.com/scfurl/piperna@dev --include-deps piperna --pip-args '--trusted-host pypi.org --trusted-host files.pythonhosted.org'


**At the FHCRC do the following...**

module load Python/3.6.7-foss-2016b-fh1
python3 -m pip install --user pipx
python3 -m pipx ensurepath
pipx install --include-deps piperna
pipx install --spec git+https://github.com/scfurl/piperna --include-deps piperna 


"""
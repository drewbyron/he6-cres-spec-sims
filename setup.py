from setuptools import setup, find_packages

setup(
    name="he6_cres_spec_sims",
    version="0.1.1",
    description="A package for simulating cres data.",
    url="https://github.com/drewbyron/he6-cres-spec-sims",
    author="William (Drew) Byron",
    author_email="wbyron@uw.edu",
    license="BSD",
    packages=find_packages(),
    install_requires=[
        "ipykernel",
        "pathlib2",
        "numpy>=1.21",
        "pandas>=1.3",
        "scipy>=1.7",
        "pyyaml",
        "matplotlib>=3.4 ",
    ],
    include_package_data=True,
    classifiers=[
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent"
    ],
)

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="coatings",
    version="1.0.0",
    author="pacosalces",
    author_email="pacosalces@gmail.com",
    description="Package for the design and optimization of dielectric coatings",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://git.ligo.org/francisco.carcoba/Coatings/",
    license="LICENSE",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "numpy",
        "scipy",
        "matplotlib",
        "h5py",
        "pyyaml",
        "physunits",
        "corner",
        "emcee",
    ],
    python_requires=">=3.10",
)

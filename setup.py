import setuptools
from cryptanalysis import __version__

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="cryptanalysis",
    version=__version__,
    author="room2042",
    description="Package for simple cryptanalysis during CTFs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://gitlab.com/room2042/cryptanalysis",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)

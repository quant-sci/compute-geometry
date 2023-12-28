import setuptools

with open("README.md", "r", encoding = "utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name = "computational-geometry",
    version = "0.0.1",
    author = "Kleyton da Costa",
    author_email = "kleyton.vsc@gmail.com",
    description = "Computational Geometry Library",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    #url = "",
    #project_urls = {
    #    "Bug Tracker": "package issues URL",
    #},
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir = {"": "cgeom"},
    packages = setuptools.find_packages(where="cgeom"),
    python_requires = ">=3.6"
)
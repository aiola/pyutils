import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
        name="pyutils",
        packages=["pyutils"],
        version="v0.1-alpha",
        description="Python utilities",
        long_description=long_description,
        long_description_content_type="text/markdown",
        url="http://github.com/aiola/pyutils",
        author="Salvatore Aiola",
        author_email="salvatore.aiola@cern.ch",
        license="Apache 2.0",
        classifiers=[
                "Development Status :: 1 - Alpha",
                "Programming Language :: Python :: 2",
                "License :: OSI Approved :: Apache 2.0 License"
                ],
        install_requires=[
        "pandas",
        "numpy",
        "scipy",
        "matplotlib",
        "uproot"
        ],
        include_package_data=True)

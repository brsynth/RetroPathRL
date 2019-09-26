from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
        name="rp3",
        version="0.0",
        author="Mathilde Koch",
        author_email="mathilde.koch@inra.fr",
        description="Perform retrosynthesis with Monte-Carlo Tree Search algorithm",
	long_description=long_description,
    	long_description_content_type="text/markdown",
        url="https://github.com/brsynth/RetroPath3",
        packages=find_packages(),
	python_requires='>=3.6',
	include_package_data=True
)

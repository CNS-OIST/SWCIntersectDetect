import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="swc_intersect_detect",
    version="1.0.0",
    author="Weiliang Chen",
    author_email="w.chen@oist.jp",
    description="Detect branch intersections in a SWC/NEURON morphology and label the intersections for curation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/CNS-OIST/SWCIntersectDetect",
    packages=setuptools.find_packages(),
    install_requires=["numpy", "networkx", "rtree"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPL-3.0",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
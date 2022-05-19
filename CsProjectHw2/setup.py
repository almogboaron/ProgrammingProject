from setuptools import setup,find_packages, Extension
setup(
    name = "mykmeanssp",
    version="0.1.0",
    description="Module for ex2 SP",
    install_requires=["invoke"],
    packages=find_packages(),
    ext_modules=[
        ## (name of module, [name of c file])
        Extension("mykmeanssp",["kmeans.c"])
    ]
)
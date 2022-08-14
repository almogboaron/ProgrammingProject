from setuptools import setup,find_packages, Extension
setup(
    name = "spkmeansmodule",
    version="0.1.0",
    description="Module for final project SP",
    install_requires=["invoke"],
    packages=find_packages(),
    ext_modules=[
        ## (name of module, [name of c file])
        Extension("spkmeansmodule",["spkmeansmodule.c"])
    ]
)
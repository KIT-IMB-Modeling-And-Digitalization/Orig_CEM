from setuptools import setup, find_packages

setup(
    name="package-test",
    version="0.0.1",
    package_dir={"": "src"},              # use the src/ layout
    packages=find_packages(where="src"),  # find packages under src/
    python_requires=">=3.9",
    include_package_data=True,
    package_data={"package_test": ["_bin/*"]},
)

# This will try to import setuptools. If not here, it will reach for the
# embedded ez_setup (or the ez_setup package). If none, it fails with a message
try:
    from setuptools import setup
except ImportError:
    try:
        import ez_setup

        ez_setup.use_setuptools()
    except ImportError:
        raise ImportError(
            "DnaWeaver could not be installed, probably because"
            " neither setuptools nor ez_setup are installed on"
            "this computer. \nInstall ez_setup "
            "([sudo] pip install ez_setup) and try again."
        )

from setuptools import setup, find_packages

exec(open("dnaweaver/version.py").read())  # loads __version__

setup(
    name="DnaWeaver",
    version=__version__,
    author="Zulko",
    description="Make ordering and assembly plans for DNA sequences",
    url='https://github.com/Edinburgh-Genome-Foundry/DnaWeaver',
    long_description=open("pypi-readme.rst").read(),
    license="MIT",
    keywords="DNA optimization assembly ordering synthetic biology",
    packages=find_packages(exclude="docs"),
    include_package_data=True,
    install_requires=[
        "numpy",
        "Biopython",
        "proglog",
        "networkx",
        "flametree",
        "dna_features_viewer",
        "weasyprint",
        "pandas",
        "jinja2"
    ],
)

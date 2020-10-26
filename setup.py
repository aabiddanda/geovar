"""Setup for GeoVar package."""

from setuptools import setup

version = "1.0.0"

required = open("requirements.txt").read().split("\n")
with open("README.md", "r") as fh:
    long_description = fh.read()

data_files = [
    "data/integrated_call_samples_v3.20130502.1kg_pops.panel",
    "data/integrated_call_samples_v3.20130502.1kg_superpops.panel",
    "data/new_1kg_nyc_hg38_filt_total.biallelic_snps.superpops_amended2.ncat3x.filt_0.geodist_cnt.txt.gz",
    "data/new_1kg_nygc.chr22.biallelic_snps.filt.n5000.freq.txt",
    "data/new_1kg_nygc.chr22.biallelic_snps.filt.n5000.vcf.gz",
    "data/new_1kg_nygc.chr22.biallelic_snps.filt.n5000.vcf.gz.tbi",
]


setup(
    name="geovar",
    version=version,
    description="A package to provide functions to visualize joint allele frequency spectra",
    author="[aabiddanda]",
    author_email="aabiddanda@gmail.com",
    url="https://github.com/aabiddanda/geovar",
    packages=["geovar"],
    install_requires=required,
    long_description=long_description,
    long_description_content_type="text/markdown",
    include_package_data=True,
    package_data={"": data_files},
    license="MIT",
)

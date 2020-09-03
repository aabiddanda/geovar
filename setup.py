#!/usr/bin/env python

from setuptools import setup

version = '0.1.0'

required = open('requirements.txt').read().split('\n')
with open('README.md', 'r') as fh:
    long_description = fh.read()

setup(
    name='geovar',
    version=version,
    description='A package to provide functions to visualize joint allele frequency spectra[D[D[D[D[D[D[D[D[D[D[D[D[D',
    author='[aabiddanda]',
    author_email='aabiddanda@gmail.com',
    url='https://github.com/aabiddanda/geovar',
    packages=['geovar'],
    install_requires=required,
    long_description=long_description,
    long_description_content_type='text/markdown',
    include_package_data=True,
    package_data={"": ["data/integrated_call_samples_v3.20130502.1kg_pops.indiv_label",
                       "data/integrated_call_samples_v3.20130502.1kg_superpops.indiv_label",
                       "data/new_1kg_nyc_hg38_filt_total.biallelic_snps.superpops_amended2.ncat3x.filt_0.geodist_cnt.txt.gz",
                       "data/new_1kg_nygc.chr22.biallelic_snps.filt.n5000.freq.txt",
                       "data/new_1kg_nygc.chr22.biallelic_snps.filt.n5000.vcf.gz",
                       "data/new_1kg_nygc.chr22.biallelic_snps.filt.n5000.vcf.gz.tbi"]},
    license='MIT'
)

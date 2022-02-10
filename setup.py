# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
setup(
  name = 'firepy',
  packages = find_packages(), 
  version = '0.1.0',
  license = 'GPLv3',
  description = 'FIRE simulator',
  author = 'Kevin A.G. Smet',
  author_email = 'ksmet1977@gmail.com',
  url = 'https://github.com/ksmet1977/firepy',
  download_url = 'https://github.com/ksmet1977/firepy/archive/0.1.0.tar.gz',
  keywords = ['FIRE', 'Financial Independence', 'Retire Early', 'Financial', 'Independence', 'Retire','Early'], 
  install_requires=[
        'numpy',
		'matplotlib',
      ],
  package_data={'firepy': []},
  include_package_data = True,
  classifiers=[
    'Development Status :: 3 - Alpha',
    'Programming Language :: Python :: 3',
    ],  
  python_requires='>=3.5',
)
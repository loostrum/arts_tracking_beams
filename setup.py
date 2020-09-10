#!/usr/bin/env python

import os
from setuptools import setup, find_packages


with open('README.md') as f:
    readme = f.read()


with open(os.path.join('arts_tracking_beams', '__version__.py')) as version_file:
    version = {}
    exec(version_file.read(), version)
    project_version = version['__version__']


setup(name='arts_tracking_beams',
      version=project_version,
      description='Create a tracking beam from ARTS tied-array beam data',
      long_description=readme,
      long_description_content_type='text/markdown',
      url='http://github.com/loostrum/arts_tracking_beams',
      author='Leon Oostrum',
      author_email='oostrum@astron.nl',
      license='Apache Software License 2.0',
      packages=find_packages(),
      zip_safe=False,
      install_requires=['numpy',
                        'astropy',
                        'tqdm',
                        'joblib'],
      entry_points={'console_scripts':
                    ['arts_create_tracking_beam=arts_tracking_beams.create_tracking_beam:main_with_args']},
      classifiers=['License :: OSI Approved :: Apache Software License',
                   'Programming Language :: Python :: 3',
                   'Operating System :: OS Independent'],
      python_requires='>=3.6'
      )

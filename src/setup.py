from setuptools import setup,find_packages

setup(
  name = 'RT_prediction',         # How you named your package folder (MyLib)
  packages = ["RT_prediction"],   # Chose the same as "name"
  version = '0.0.1',      # Start with a small number and increase it with every change you make
  license='MIT',        # Chose a license from here: https://help.github.com/articles/licensing-a-repository
  description = 'Random Forest Based Retention Time Prediction',   # Give a short description about your library
  author = 'Ethan Stancliffe',                   # Type in your name
  author_email = 'estancl1234@gmail.com',      # Type in your E-Mail
  url = 'https://github.com/e-stan/RT_prediction',   # Provide either the link to your github or to your website
  download_url = 'https://github.com/e-stan/RT_prediction/archive/v0.0.1.tar.gz',    # I explain this later on
  keywords = ['Metabolomics', 'Retention Time',"Machine Learning"],   # Keywords that define your package best
  install_requires=[            # I get to this in a second
          'numpy',
          'sklearn',
          'pandas',
      ],
  classifiers=[
    'Development Status :: 3 - Alpha',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
    'Intended Audience :: Developers',      # Define that your audience are developers
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',   # Again, pick a license
    'Programming Language :: Python :: 3.7',      #Specify which pyhton versions that you want to support
  ],
  package_data={
    "":["molecular_properties.csv","compute_molecular_properties.R"]
    }
)
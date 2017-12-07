[//]: # ( [![PyPI package](https://img.shields.io/pypi/v/Kinetic-Datanator.svg)](https://pypi.python.org/pypi/Kinetic-Datanator) )
[![Documentation](https://readthedocs.org/projects/kinetic_datanator/badge/?version=latest)](http://kinetic_datanator.readthedocs.org)
[![Test results](https://circleci.com/gh/KarrLab/kinetic_datanator.svg?style=shield)](https://circleci.com/gh/KarrLab/kinetic_datanator)
[![Test coverage](https://coveralls.io/repos/github/KarrLab/kinetic_datanator/badge.svg)](https://coveralls.io/github/KarrLab/kinetic_datanator)
[![Code analysis](https://api.codeclimate.com/v1/badges/62e495c53a118f35afea/maintainability)](https://codeclimate.com/github/KarrLab/kinetic_datanator)
[![License](https://img.shields.io/github/license/KarrLab/kinetic_datanator.svg)](LICENSE)
![Analytics](https://ga-beacon.appspot.com/UA-86759801-1/kinetic_datanator/README.md?pixel)

# `kinetic_datanator`: tools for aggregating data for biochemical models

## Contents
* [Overview](#overview)
* [Installation](#installation)
* [Usage](#usage)
* [Documentation](#documentation)
* [Tests](#tests)
* [License](#license)
* [Development team](#development-team)
* [Questions and comments](#questions-and-comments)
* [References](#references)

## Overview
`kinetic_datanator` is a software tool for finding experimental data for building and calibrating dynamical models of cellular biochemistry such as metabolite, RNA, and protein abundances; protein complex compositions; transcription factor binding motifs; and kinetic parameters. ``kinetic_datanator`` is particularly useful for building large models, such as whole-cell models, that require large amounts of data to constrain large numbers of parameters. ``kinetic_datanator`` was motivated by the need for large amounts of data to constrain whole-cell models and the fact that this data is hard to utilize because it is scattered across numerous siloed repositories.

## Installation instructions, user instructions, and code documentation
Please see the [documentation](http://kinetic_datanator.readthedocs.org) at the Read the Docs for installation instructions, user instructions, and code documentation.

## Testing `kinetic_datanator`
To ensure `kinetic_datanator` works properly, we have developed extensive units tests of every aspect of `kinetic_datanator`. We recommend using `pytest` to run these tests as follows:

```
python -m pytest tests
```

## License
This software is released open-source under the [MIT license](LICENSE).

## Development team
The model was developed by the following individuals in the [Karr Lab](http://www.karrlab.org) at the Icahn School of Medicine at Mount Sinai in New York, USA.

* Yosef Roth
* Balazs Szigeti
* Saahith Pochiraju
* Jonathan Karr

## Questions and comments
Please contact the [Karr Lab](http://www.karrlab.org) with any questions or comments.

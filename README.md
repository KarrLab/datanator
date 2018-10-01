[![PyPI package](https://img.shields.io/pypi/v/datanator.svg)](https://pypi.python.org/pypi/datanator) 
[![Documentation](https://readthedocs.org/projects/datanator/badge/?version=latest)](http://docs.karrlab.org/datanator)
[![Test results](https://circleci.com/gh/KarrLab/datanator.svg?style=shield)](https://circleci.com/gh/KarrLab/datanator)
[![Test coverage](https://coveralls.io/repos/github/KarrLab/datanator/badge.svg)](https://coveralls.io/github/KarrLab/datanator)
[![Code analysis](https://api.codeclimate.com/v1/badges/e9b796130e29aee4672f/maintainability)](https://codeclimate.com/github/KarrLab/datanator)
[![License](https://img.shields.io/github/license/KarrLab/datanator.svg)](LICENSE)
![Analytics](https://ga-beacon.appspot.com/UA-86759801-1/datanator/README.md?pixel)

# Datanator: Toolkit for discovering and aggregating data for whole-cell modeling

## Contents
* [Overview](#overview)
* [Installation instructions and documentation](#installation-instructions-and-documentation)
* [Testing Datanator](#testing-datanator)
* [License](#license)
* [Development team](#development-team)
* [Questions and comments](#questions-and-comments)

## Overview
Extensive data is needed to build comprehensive predictive models of cells. Although the literature and public repositories contain extensive data about cells, this data is hard to utilize for modeling because it is scattered across a large number of sources; because it is described with inconsistent identifiers, units, and data models; and because there are few tools for finding relevant data for modeling specific species and environmental conditions. 

Datanator is a software tool for discovering, aggregating, and integrating the data needed for modeling cells. This includes metabolite, RNA, and protein abundances; protein complex compositions; transcription factor binding motifs; and kinetic parameters. Datanator is particularly useful for building large models, such as whole-cell models, that require large amounts of data to constrain large numbers of parameters.

This package contains the source code for Datanator. The data aggregated with Datanator is available at [https://quiltdata.com/package/karrlab/datanator](https://quiltdata.com/package/karrlab/datanator).

## Installation instructions and documentation
Please see the [documentation](http://docs.karrlab.org/datanator) for installation instructions, user instructions, and code documentation. 

Note, Datanator only supports Python 3. 

## Testing Datanator
To ensure Datanator works properly, we have developed extensive units tests of every aspect of `datanator`. We recommend using `pytest` to run these tests as follows:

```
python3 -m pytest tests
```

## License
This software is released open-source under the [MIT license](LICENSE).

## Development team
The model was developed by the [Karr Lab](https://www.karrlab.org) at the Icahn School of Medicine at Mount Sinai in New York, US.

* Saahith Pochiraju
* Yosef Roth
* Balazs Szigeti
* Jonathan Karr

## Questions and comments
Please contact the [Karr Lab](https://www.karrlab.org) with any questions or comments.

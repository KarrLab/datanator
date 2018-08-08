[![PyPI package](https://img.shields.io/pypi/v/Kinetic-Datanator.svg)](https://pypi.python.org/pypi/Kinetic-Datanator) 
[![Documentation](https://readthedocs.org/projects/kinetic-datanator/badge/?version=latest)](http://docs.karrlab.org/kinetic_datanator)
[![Test results](https://circleci.com/gh/KarrLab/kinetic_datanator.svg?style=shield)](https://circleci.com/gh/KarrLab/kinetic_datanator)
[![Test coverage](https://coveralls.io/repos/github/KarrLab/kinetic_datanator/badge.svg)](https://coveralls.io/github/KarrLab/kinetic_datanator)
[![Code analysis](https://api.codeclimate.com/v1/badges/62e495c53a118f35afea/maintainability)](https://codeclimate.com/github/KarrLab/kinetic_datanator)
[![License](https://img.shields.io/github/license/KarrLab/kinetic_datanator.svg)](LICENSE)
![Analytics](https://ga-beacon.appspot.com/UA-86759801-1/kinetic_datanator/README.md?pixel)

# Kinetic Datanator: tools for aggregating data for cell modeling

## Contents
* [Overview](#overview)
* [Installation instructions and documentation](#installation-instructions-and-documentation)
* [Testing Kinetic Datanator](#testing-kinetic_datanator)
* [License](#license)
* [Development team](#development-team)
* [Questions and comments](#questions-and-comments)

## Overview
Extensive data is needed to build comprehensive predictive models of cells. Although the literature and public repositories contain extensive data about cells, this data is hard to utilize for modeling because it is scattered across a large number of sources; because it is described with inconsistent identifiers, units, and data models; and because there are few tools for finding relevant data for modeling specific species and environmental conditions. 

Kinetic Datanator is a software tool for discovering, aggregating, and integrating the data needed for modeling cells. This includes metabolite, RNA, and protein abundances; protein complex compositions; transcription factor binding motifs; and kinetic parameters. Kinetic Datanator is particularly useful for building large models, such as whole-cell models, that require large amounts of data to constrain large numbers of parameters.

This package contains the source code for Kinetic Datanator. The data aggregated with Kinetic Datanator is available at [https://quiltdata.com/package/karrlab/kinetic_datanator](https://quiltdata.com/package/karrlab/kinetic_datanator).

## Installation instructions and documentation
Please see the [documentation](http://docs.karrlab.org/kinetic_datanator) for installation instructions, user instructions, and code documentation. 

Note, Kinetic Datanator only supports Python 3. 

## Testing Kinetic Datanator
To ensure Kinetic Datanator works properly, we have developed extensive units tests of every aspect of `kinetic_datanator`. We recommend using `pytest` to run these tests as follows:

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

import datanator.config.core
from datanator.util import mongo_util
from datanator.util import file_util
from datanator.util import molecule_util
import six
import requests
from xml import etree
import libsbml
import re
import datetime
import bs4
import html
import csv
import pubchempy
import sys
import Bio.Alphabet
import Bio.SeqUtils
import math
import logging
import pymongo
import hashlib
logging.basicConfig(filename='./logs/sabiork_parser.log', level=logging.WARNING, 
                    format='%(asctime)s %(levelname)s %(name)s %(message)s')
logger=logging.getLogger()


class SabioRk:

    def __init__(self, cache_dirname=None, MongoDB=None, replicaSet=None, db=None,
                 verbose=False, max_entries=float('inf'), username=None,
                 password=None, authSource='admin', webservice_batch_size=100,
                 excel_batch_size=100):
        self.cache_dirname = cache_dirname
        self.MongoDB = MongoDB
        self.replicaSet = replicaSet
        self.db = db
        self.verbose = verbose
        self.max_entries = max_entries
        self.client, self.db_obj, self.collection = mongo_util.MongoUtil(
            MongoDB=MongoDB, db=db, username=username, password=password,
            authSource=authSource).con_db('sabio_rk_new')
        self.client, self.db_obj, self.collection_compound = mongo_util.MongoUtil(
            MongoDB=MongoDB, db=db, username=username, password=password,
            authSource=authSource).con_db('sabio_compound')
        self.excel_batch_size = excel_batch_size
        self.ENDPOINT_DOMAINS = {
            'sabio_rk': 'http://sabiork.h-its.org',
            'uniprot': 'http://www.uniprot.org',
        }
        self.ENDPOINT_KINETIC_LAWS_SEARCH = self.ENDPOINT_DOMAINS['sabio_rk'] + \
            '/sabioRestWebServices/searchKineticLaws/entryIDs'
        self.ENDPOINT_WEBSERVICE = self.ENDPOINT_DOMAINS['sabio_rk'] + \
            '/sabioRestWebServices/kineticLaws'
        self.ENDPOINT_EXCEL_EXPORT = self.ENDPOINT_DOMAINS['sabio_rk'] + \
            '/entry/exportToExcelCustomizable'
        self.ENDPOINT_COMPOUNDS_PAGE = self.ENDPOINT_DOMAINS['sabio_rk'] + \
            '/compdetails.jsp'
        self.ENDPOINT_KINETIC_LAWS_PAGE = self.ENDPOINT_DOMAINS['sabio_rk'] + \
            '/kindatadirectiframe.jsp'
        self.SKIP_KINETIC_LAW_IDS = (51286,)
        self.PUBCHEM_MAX_TRIES = 10
        self.PUBCHEM_TRY_DELAY = 0.25
        self.webservice_batch_size = webservice_batch_size
        self.file_manager = file_util.FileUtil()

    def load_content(self):
        """ Download the content of SABIO-RK and store it to a remote mongoDB. """

        ##################################
        ##################################
        # determine ids of kinetic laws
        if self.verbose:
            print('Downloading the IDs of the kinetic laws ...')

        ids = self.load_kinetic_law_ids()

        if self.verbose:
            print('  Downloaded {} IDs'.format(len(ids)))

        ##################################
        ##################################
        # remove bad IDs
        ids = list(filter(lambda id: id not in self.SKIP_KINETIC_LAW_IDS, ids))

        # sort ids
        ids.sort()

        # load only `max_entries` IDs
        if len(ids) > self.max_entries:
            ids = ids[0:self.max_entries]

        ##################################
        ##################################
        # download kinetic laws
        exisitng_ids = self.collection.distinct('kinlaw_id')
        new_ids = list(set(ids).difference(set(exisitng_ids)))
        new_ids.sort()

        if self.verbose:
            print('Downloading {} kinetic laws ...'.format(len(new_ids)))

        missing_ids = self.load_kinetic_laws(new_ids)

        if self.verbose:
            print('  done')

        # fill in missing information from Excel export
        exisitng_ids = self.collection.distinct('kinlaw_id')
        loaded_new_ids = list(set(new_ids).intersection(exisitng_ids))
        loaded_new_ids.sort()

        if self.verbose:
            print('Updating {} kinetic laws ...'.format(len(loaded_new_ids)))

        self.load_missing_kinetic_law_information_from_tsv(loaded_new_ids)

        if self.verbose:
            print('  done')

        ##################################
        ##################################
        # fill in missing information from HTML pages
        if self.verbose:
            print('Updating {} kinetic laws ...'.format(len(loaded_new_ids)))

        self.load_missing_enzyme_information_from_html(loaded_new_ids)

        if self.verbose:
            print('  done')

    def load_kinetic_law_ids(self):
        """ Download the IDs of all of the kinetic laws stored in SABIO-RK

        Returns:
            :obj:`list` of :obj:`int`: list of kinetic law IDs

        """
        # create session
        response = requests.get(self.ENDPOINT_KINETIC_LAWS_SEARCH, params={
            'q': 'DateSubmitted:01/01/2000',
        })
        response.raise_for_status()

        # get IDs of kinetic laws
        root = etree.ElementTree.fromstring(response.text)
        ids = [int(float(node.text)) for node in root.findall('SabioEntryID')]

        # sort ids
        ids.sort()

        # return IDs
        return ids

    def load_kinetic_laws(self, ids):
        """ Download kinetic laws from SABIO-RK

        Args:
            ids (:obj:`list` of :obj:`int`): list of IDs of kinetic laws to download

        Raises:
            :obj:`Error`: if an HTTP request fails
        """
        # todo: scrape strain, recombinant, experiment type information from web pages
        session = requests

        batch_size = self.webservice_batch_size
        loaded_ids = []

        for i_batch in range(int(math.ceil(float(len(ids)) / batch_size))):
            if self.verbose and (i_batch % max(1, 100. / batch_size) == 0):
                print('  Downloading kinetic laws {}-{} of {} in SBML format'.format(
                    i_batch * batch_size + 1,
                    min(len(ids), i_batch * batch_size + max(100, batch_size)),
                    len(ids)))

            batch_ids = ids[i_batch *
                            batch_size:min((i_batch + 1) * batch_size, len(ids))]
            response = session.get(self.ENDPOINT_WEBSERVICE, params={
                'kinlawids': ','.join(str(id) for id in batch_ids),
            })

            response.raise_for_status()
            if not response.text:
                raise Exception('Unable to download kinetic laws with ids {}'.format(
                    ', '.join([str(id) for id in batch_ids])))

            loaded_ids += self.create_kinetic_laws_from_sbml(batch_ids,
                                                            response.content if six.PY2 else response.text)

        not_loaded_ids = list(set(ids).difference(loaded_ids))
        if not_loaded_ids:
            not_loaded_ids.sort()
            warning = 'Several kinetic laws were not found:\n  {}'.format(
                '\n  '.join([str(id) for id in not_loaded_ids]))
            logging.warning(warning)
        return not_loaded_ids

    def create_kinetic_laws_from_sbml(self, ids, sbml):
        """ Add kinetic laws defined in an SBML file to the local mongodb database

        Args:
            ids (:obj:`list` of :obj:`int`): list kinetic law IDs
            sbml (:obj:`str`): SBML representation of one or more kinetic laws (root)

        Returns:
            :obj:`tuple`:

                * :obj:`list` of :obj:`KineticLaw`: list of kinetic laws
                * :obj:`list` of :obj:`Compound` or :obj:`Enzyme`: list of species (compounds or enzymes)
                * :obj:`list` of :obj:`Compartment`: list of compartments
        """
        reader = libsbml.SBMLReader()
        doc = reader.readSBMLFromString(sbml)
        model = doc.getModel()

        functions = {}
        functions_sbml = model.getListOfFunctionDefinitions()
        for i_function in range(functions_sbml.size()):
            function_sbml = functions_sbml.get(i_function)
            math_sbml = function_sbml.getMath()
            if math_sbml.isLambda() and math_sbml.getNumChildren():
                eq = libsbml.formulaToL3String(
                    math_sbml.getChild(math_sbml.getNumChildren() - 1))
            else:
                eq = None
            if eq in ('', 'NaN'):
                eq = None
            functions[function_sbml.getId()] = eq

        units = {}
        units_sbml = model.getListOfUnitDefinitions()
        for i_unit in range(units_sbml.size()):
            unit_sbml = units_sbml.get(i_unit)
            units[unit_sbml.getId()] = unit_sbml.getName()

        # species
        specie_properties = {}
        species_sbml = model.getListOfSpecies()
        species = []
        for i_specie in range(species_sbml.size()):
            specie_sbml = species_sbml.get(i_specie)
            specie, properties = self.get_specie_from_sbml(specie_sbml)
            species.append(specie)
            specie_properties[specie_sbml.getId()] = properties

        # kinetic laws
        reactions_sbml = model.getListOfReactions()
        if reactions_sbml.size() != len(ids):
            raise ValueError('{} reactions {} is different from the expected {}'.format(
                reaction_sbml.size(), len(ids)))
        # kinetic_laws = []
        loaded_ids = []
        for i_reaction, _id in enumerate(ids):
            reaction_sbml = reactions_sbml.get(i_reaction)
            kinetic_law = self.create_kinetic_law_from_sbml(
                _id, reaction_sbml, species, specie_properties, functions, units)
            try:
                self.collection.update_one({'kinlaw_id': _id},
                                           {'$set': kinetic_law},
                                           upsert=True)
            except pymongo.errors.WriteError as err:
                logger.error(err)
            loaded_ids.append(_id)
        return loaded_ids

    def create_kinetic_law_from_sbml(self, id, sbml, root_species, specie_properties, functions, units):
        """ Make a kinetic law doc for mongoDB

        Args:
            id (:obj:`int`): identifier
            sbml (:obj:`libsbml.KineticLaw`): SBML-representation of a reaction (reaction_sbml)
            species (:obj:`list`): list of species in root sbml
            specie_properties (:obj:`dict`): additional properties of the compounds/enzymes

                * `is_wildtype` (:obj:`bool`): indicates if the enzyme is wildtype or mutant
                * `variant` (:obj:`str`): description of the variant of the eznyme
                * `modifier_type` (:obj:`str`): type of the enzyme (e.g. Modifier-Catalyst)

            functions (:obj:`dict` of :obj:`str`: :obj:`str`): dictionary of rate law equations (keys = IDs in SBML, values = equations)
            units (:obj:`dict` of :obj:`str`: :obj:`str`): dictionary of units (keys = IDs in SBML, values = names)

        Returns:
            :obj:`dictionary`: kinetic law

        Raises:
            :obj:`ValueError`: if the temperature is expressed in an unsupported unit
        """
        law = sbml.getKineticLaw()
        x_refs = self.create_cross_references_from_sbml(law)
        reaction_x_refs = self.create_cross_references_from_sbml(sbml)
        kinetic_law = {}
        # stop if kinetic law entry is empty
        if not law.getMetaId():
            return None

        """ participants """
        kinetic_law['reactants'] = []
        reactants = sbml.getListOfReactants()
        for i_part in range(reactants.size()):
            part_sbml = reactants.get(i_part)
            compound, compartment = self.get_specie_reference_from_sbml(
                part_sbml.getSpecies(), root_species)
            compound = self.load_compounds(compound)
            if 'structures' not in compound[0].keys():
                compound = self.infer_compound_structures_from_names(compound)
            part = {
                'compartment': compartment,
                'coefficient': part_sbml.getStoichiometry()}
            react = {**compound[0], **part}
            kinetic_law['reactants'].append(react)

        kinetic_law['products'] = []
        products = sbml.getListOfProducts()
        for i_part in range(products.size()):
            part_sbml = products.get(i_part)
            compound, compartment = self.get_specie_reference_from_sbml(
                part_sbml.getSpecies(), root_species)
            compound = self.load_compounds(compound)
            if 'structures' not in compound[0].keys():
                compound = self.infer_compound_structures_from_names(compound)
            part = {
                'compartment': compartment,
                'coefficient': part_sbml.getStoichiometry()}
            prod = {**compound[0], **part}
            kinetic_law['products'].append(prod)

        """ cross references """
        # Note: these are stored KineticLaws rather than under Reactions because this seems to how SABIO-RK stores this information.
        # For example, kinetic laws 16016 and 28003 are associated with reaction 9930, but they have different EC numbers 1.1.1.52 and
        # 1.1.1.50, respectively.
        kinetic_law['cross_references'] = list(
            filter(lambda x_ref: list(x_ref.keys())[0] not in ['taxonomy'], reaction_x_refs))

        # rate_law
        kinetic_law['equation'] = functions[law.getMetaId()[5:]]

        # parameters
        kinetic_law['parameters'] = []
        params = law.getListOfLocalParameters()
        for i_param in range(params.size()):
            param = params.get(i_param)

            match = re.match(
                r'^(.*?)_((SPC|ENZ)_([0-9]+)_(.*?))$', param.getId(), re.IGNORECASE)
            if match:
                observed_name = match.group(1)
                species, compartment = self.get_specie_reference_from_sbml(
                    match.group(2), root_species)
                if 'subunits' in species[0].keys():
                    compound = None
                    enzyme = species[0]
                else:
                    compound = species[0]
                    enzyme = None
            else:
                observed_name = param.getId()
                compound = None
                enzyme = None
                compartment = None
            observed_name = observed_name.replace('div', '/')

            observed_type = param.getSBOTerm()

            observed_units_id = param.getUnits()
            if observed_units_id:
                if observed_units_id in units:
                    observed_units = units[observed_units_id]
                else:
                    observed_units = observed_units_id
            else:
                observed_units = None

            observed_value = param.getValue()

            parameter = {
                'compound': compound,
                'enzyme': enzyme,
                'compartment': compartment,
                'observed_name': observed_name,
                'observed_type': observed_type,
                'observed_value': observed_value,
                'observed_units': observed_units,
                'modified': datetime.datetime.utcnow()
            }
            kinetic_law['parameters'].append(parameter)

        # modifiers to kinetic law
        kinetic_law['modifiers'] = []
        modifiers = sbml.getListOfModifiers()
        for i_modifier in range(modifiers.size()):
            modifier = modifiers.get(i_modifier)
            modifier_id = modifier.getSpecies()
            specie, compartment = self.get_specie_reference_from_sbml(
                modifier_id, root_species)
            type = specie_properties[modifier.getSpecies()]['modifier_type']
            if modifier_id[0:3] == 'SPC':
                part = {
                    'compartment': compartment,
                    'type': type
                }  # ReactionParticipant
                modif = {**specie[0], **part}
                kinetic_law['modifiers'].append(modif)
            elif modifier_id[0:3] == 'ENZ':
                kinetic_law['enzyme'], kinetic_law['enzyme_compartment'] = self.get_specie_reference_from_sbml(
                    modifier_id, root_species)
                kinetic_law['enzyme_type'] = specie_properties[modifier.getSpecies(
                )]['modifier_type']
                kinetic_law['taxon_wildtype'] = specie_properties[modifier_id]['is_wildtype']
                kinetic_law['taxon_variant'] = specie_properties[modifier_id]['variant']

        # taxon
        taxon = self.file_manager.search_dict_list(reaction_x_refs, 'taxonomy')
        if len(taxon) > 0:
            kinetic_law['taxon'] = int(taxon[0]['taxonomy'])
        else:
            kinetic_law['taxon'] = None

        """ conditions """
        conditions = law \
            .getAnnotation() \
            .getChild('sabiork') \
            .getChild('experimentalConditions')

        # temperature
        if conditions.hasChild('temperature'):
            temperature = conditions \
                .getChild('temperature') \
                .getChild('startValueTemperature') \
                .getChild(0) \
                .getCharacters()
            temperature = float(temperature)
            temperature_units = conditions \
                .getChild('temperature') \
                .getChild('temperatureUnit') \
                .getChild(0) \
                .getCharacters()
            if temperature_units not in ['°C', '��C']:
                raise ValueError(
                    'Unsupported temperature units: {}'.format(temperature_units))
            kinetic_law['temperature'] = temperature

        # pH
        if conditions.hasChild('pH'):
            ph = conditions \
                .getChild('pH') \
                .getChild('startValuepH') \
                .getChild(0) \
                .getCharacters()
            kinetic_law['ph'] = float(ph)

        # media
        if conditions.hasChild('buffer'):
            media = conditions \
                .getChild('buffer') \
                .getChild(0) \
                .getCharacters() \
                .strip()
            if six.PY2:
                media = unicode(media.decode('utf-8'))
            kinetic_law['media'] = media

        """ references """
        kinetic_law['references'] = list(
            filter(lambda x_ref: list(x_ref.keys())[0] != 'sabiork.kineticrecord', x_refs))

        """ updated """
        kinetic_law['modified'] = datetime.datetime.utcnow()

        return kinetic_law

    def get_compartment_from_sbml(self, sbml):
        """ get compartment from sbml

        Args:
            sbml (:obj:`libsbml.Compartment`): SBML-representation of a compartment

        Returns:
            :dictionary: compartment
        """
        name = sbml.getName()
        if name == 'Cell':
            return None

        modified = datetime.datetime.utcnow()
        compartment = {'name': name, 'modified': modified}

        return compartment

    def get_specie_reference_from_sbml(self, specie_id, species):
        """ Get the compound/enzyme associated with an SBML species by its ID

        Args:
            specie_id (:obj:`str`): ID of an SBML species

        Returns:
            :obj:`tuple`:

                * :obj:`Compound` or :obj:`Enzyme`: compound or enzyme
                * :obj:`Compartment`: compartment

        Raises:
            :obj:`ValueError`: if the species is not a compound or enzyme, no species
                with `id` = `specie_id` exists, or no compartment with `name` = `compartment_name`
                exists
        """
        tmp = specie_id.split('_')
        type = tmp[0]
        specie_id = int(float(tmp[1]))
        compartment_name = '_'.join(tmp[2:])

        if type == 'SPC':
            specie = self.file_manager.search_dict_list(
                species, '_id', value=specie_id)
            self.collection_compound.update_one({'_id': specie_id},
                                                {'$set': specie[0]}, upsert=True)
        elif type == 'ENZ':
            specie = self.file_manager.search_dict_list(
                species, '_id', value=specie_id)
        else:
            raise ValueError('Unsupported species type: {}'.format(type))

        if compartment_name != 'Cell':
            compartment = compartment_name
        else:
            compartment = None

        return (specie, compartment)

    def get_specie_from_sbml(self, sbml):
        """ get species information from sbml

        Args:
            sbml (:obj:`libsbml.Species`): SBML-representation of a compound or enzyme

        Returns:
            :obj:`tuple`:

                * :obj:`Compound`: or :obj:`Enzyme`: compound or enzyme
                * :obj:`dict`: additional properties of the compound/enzyme

                    * `is_wildtype` (:obj:`bool`): indicates if the enzyme is wildtype or mutant
                    * `variant` (:obj:`str`): description of the variant of the eznyme
                    * `modifier_type` (:obj:`str`): type of the enzyme (e.g. Modifier-Catalyst)

        Raises:
            :obj:`ValueError`: if a species is of an unsupported type (i.e. not a compound or enzyme)
        """
        # id, name
        type_id_compartment = sbml.getId().split('_')
        type = type_id_compartment[0]
        id = int(float(type_id_compartment[1]))

        # modifier type
        modifier_type = ''
        if sbml.isSetAnnotation() \
                and sbml.getAnnotation().hasChild('sabiork') \
                and sbml.getAnnotation().getChild('sabiork').hasChild('modifierType'):
            modifier_type = sbml \
                .getAnnotation() \
                .getChild('sabiork') \
                .getChild('modifierType') \
                .getChild(0) \
                .getCharacters()

        # create object or return existing object
        if type == 'SPC':
            name = sbml.getName()
            properties = {'modifier_type': modifier_type}
            specie = {'_id': id, 'name': name}
            self.collection_compound.update_one({'_id': id},
                                                {'$set': {'name': name}},
                                                upsert=True)
        elif type == 'ENZ':
            name, is_wildtype, variant = self.parse_enzyme_name(sbml.getName())
            if six.PY2:
                variant = unicode(variant.decode('utf-8'))
            properties = {'is_wildtype': is_wildtype,
                          'variant': variant, 'modifier_type': modifier_type}

            specie = {'_id': id, 'name': name}
        else:
            raise ValueError('Unsupported species type: {}'.format(type))

        # cross references
        cross_references = self.create_cross_references_from_sbml(sbml)
        if type == 'SPC':
            specie['cross_references'] = cross_references
            self.collection_compound.update_one({'_id': id},
                                                {'$set': {
                                                    'cross_references': cross_references}},
                                                upsert=True)
        elif type == 'ENZ':
            specie['subunits'] = []
            specie['cross_references'] = []
            for cross_reference in cross_references:
                if 'uniprot' in cross_reference:
                    specie['subunits'].append(cross_reference)
                else:
                    specie['cross_references'].append(cross_reference)

        # updated
        specie['modified'] = datetime.datetime.utcnow()

        return (specie, properties)

    def parse_enzyme_name(self, sbml):
        """ Parse the name of an enzyme in SBML for the enzyme name, wild type status, and variant
        description that it contains.

        Args:
            sbml (:obj:`str`): enzyme name in SBML

        Returns:
            :obj:`tuple`:

                * :obj:`str`: name
                * :obj:`bool`: if :obj:`True`, the enzyme is wild type
                * :obj:`str`: variant

        Raises:
            :obj:`ValueError`: if the enzyme name is formatted in an unsupport format
        """
        match = re.match(
            r'^(.*?)\(Enzyme\) (wildtype|mutant),?(.*?)$', sbml, re.IGNORECASE)
        if match:
            name = match.group(1)
            is_wildtype = match.group(2).lower() == 'wildtype'
            variant = match.group(3).strip()
            return (name, is_wildtype, variant)

        match = re.match(
            r'^Enzyme (wildtype|mutant),?( (.*?))*$', sbml, re.IGNORECASE)
        if match:
            if match.group(3):
                name = match.group(3).strip()
            else:
                name = None
            is_wildtype = match.group(1).lower() == 'wildtype'
            variant = None
            return (name, is_wildtype, variant)

        match = re.match(r'^Enzyme - *$', sbml, re.IGNORECASE)
        if match:
            name = None
            is_wildtype = True
            variant = None
            return (name, is_wildtype, variant)

        raise ValueError('Cannot parse enzyme name: {}'.format(sbml))

    def create_cross_references_from_sbml(self, sbml):
        """ Look up cross references from an SBML object to dictionary

        Args:
            sbml (:obj:`libsbml.SBase`): object in an SBML documentation

        Returns:
            :obj:`list` of dictionary: list of resources
        """
        if not sbml.isSetAnnotation():
            return []

        xml = sbml.getAnnotation().getChild('RDF').getChild('Description')

        attr = libsbml.XMLTriple(
            'resource', 'http://www.w3.org/1999/02/22-rdf-syntax-ns#', 'rdf')

        x_refs = []
        for i_child in range(xml.getNumChildren()):
            bag = xml.getChild(i_child).getChild('Bag')
            for i_li in range(bag.getNumChildren()):
                li = bag.getChild(i_li)

                val = li.getAttrValue(attr)
                if val.startswith('http://identifiers.org/ec-code/Complex '):
                    _, _, ids = val.partition(' ')
                    resources = [('ec-code', id) for id in ids.split('/')]
                else:
                    parsed_url = val.split('/')
                    if len(parsed_url) == 5:
                        resources = [(parsed_url[3], parsed_url[4])]
                    else:
                        resources = [(parsed_url[2], parsed_url[3])]

                for namespace, id in resources:
                    resource = {namespace: id}

                if resource not in x_refs:
                    x_refs.append(resource)

        return x_refs

    def load_compounds(self, compounds=None):
        """ Download information from SABIO-RK about all of the compounds stored sabio_compounds
                collection
        Args:
            compounds (:obj:`list` of :obj:`obj`): list of compounds to download

        Raises:
            :obj:`Error`: if an HTTP request fails
        """
        if compounds is None:
            compounds = self.collection_compound.find({})
            n_compounds = self.collection_compound.count_documents({})
        else:
            n_compounds = len(compounds)
        result = []

        for i_compound, c in enumerate(compounds):
            # print status
            # if self.verbose and (i_compound % 100 == 0):
            #     print('  Downloading compound {} of {}'.format(
            #         i_compound + 1, n_compounds))

            # download info
            response = requests.get(
                self.ENDPOINT_COMPOUNDS_PAGE, params={'cid': c['_id']})
            response.raise_for_status()

            # parse info
            doc = bs4.BeautifulSoup(response.text, 'html.parser')
            table = doc.find('table')

            # name
            node = table.find('span', id='commonName')
            if node:
                c['name'] = node.get_text()

            # synonyms
            c['synonyms'] = []
            synonym_label_node = table.find('b', text='Synonyms')
            if synonym_label_node:
                for node in list(synonym_label_node.parents)[1].find_all('span'):
                    name = node.get_text()
                    c['synonyms'].append(name)

            # structure
            c['structures'] = []

            inchi_label_node = table.find('b', text='InChI')
            if inchi_label_node:
                for node in list(inchi_label_node.parents)[1].find_all('span'):
                    value = node.get_text()
                    inchi = {'inchi': value}
                    norm = self.calc_inchi_formula_connectivity(inchi)
                    c['structures'].append({**inchi, **norm})

            smiles_label_node = table.find('b', text='SMILES')
            if smiles_label_node:
                for node in list(smiles_label_node.parents)[1].find_all('span'):
                    value = node.get_text()
                    smiles = {'smiles': value}
                    # norm = self.calc_inchi_formula_connectivity(smiles)
                    c['structures'].append(smiles)

            # cross references
            c['cross_references'] = []
            for node in table.find_all('a'):

                url = node.get('href')

                id = node.get_text()

                if url.startswith('http://www.genome.jp/dbget-bin/www_bget?cpd:'):
                    namespace = 'kegg.compound'
                elif url.startswith('http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?sid='):
                    namespace = 'pubchem.substance'
                elif url.startswith('http://www.ebi.ac.uk/chebi/searchId.do?chebiId='):
                    namespace = 'chebi'
                    id = 'CHEBI:' + id
                elif url.startswith('https://reactome.org/content/detail/'):
                    namespace = 'reactome'
                elif url.startswith('https://biocyc.org/compound?orgid=META&id='):
                    namespace = 'biocyc'
                elif url.startswith('https://www.metanetx.org/chem_info/'):
                    namespace = 'metanetx.chemical'
                elif url.startswith('http://www.chemspider.com/Chemical-Structure.'):
                    namespace = 'chemspider'
                elif url.startswith('http://sabiork.h-its.org/newSearch?q=sabiocompoundid:'):
                    continue
                else:
                    namespace = html.unescape(node.parent.parent.parent.find_all('td')[
                                              0].get_text()).strip()
                    warning = 'Compound {} has unkonwn cross reference type to namespace {}'.format(c['_id'],
                                                                                                        namespace)
                    logging.warning(warning)
                resource = {namespace: id}

                c['cross_references'].append(resource)

            # udated
            c['modified'] = datetime.datetime.utcnow()
            self.collection_compound.update_one({'_id': c['_id']},
                                                {'$set': {'synonyms': c['synonyms'],
                                                          'structures': c['structures'],
                                                          'cross_references': c['cross_references']}},
                                                upsert=True)
            result.append(c)
        return result 

    def load_missing_kinetic_law_information_from_tsv(self, ids):
        """ Update the properties of kinetic laws in mongodb based on content downloaded
        from SABIO in TSV format.

        Args:
            ids (:obj:`list` of :obj:`int`): list of IDs of kinetic laws to download
        """
        batch_size = self.excel_batch_size

        for i_batch in range(int(math.ceil(float(len(ids)) / batch_size))):
            if self.verbose:
                print('  Downloading kinetic laws {}-{} of {} in Excel format'.format(
                    i_batch * batch_size + 1,
                    min(len(ids), (i_batch + 1) * batch_size),
                    len(ids)))

            batch_ids = ids[i_batch *
                            batch_size:min((i_batch + 1) * batch_size, len(ids))]
            response = requests.get(self.ENDPOINT_EXCEL_EXPORT, params={
                'entryIDs[]': batch_ids,
                'fields[]': [
                    'EntryID',
                    'KineticMechanismType',
                    'Tissue',
                    'Parameter',
                ],
                'preview': False,
                'format': 'tsv',
                'distinctRows': 'false',
            })
            response.raise_for_status()
            if not response.text:
                cache = session.cache
                key = cache.create_key(response.request)
                cache.delete(key)
                raise Exception('Unable to download kinetic laws with ids {}'.format(
                    ', '.join([str(id) for id in batch_ids])))

            self.load_missing_kinetic_law_information_from_tsv_helper(
                response.text)

    def load_missing_kinetic_law_information_from_tsv_helper(self, tsv):
        """ Update the properties of kinetic laws in the mongodb based on content downloaded
        from SABIO in TSV format.

        Note: this method is necessary because neither of SABIO's SBML and Excel export methods provide
        all of the SABIO's content.

        Args:
            tsv (:obj:`str`): TSV-formatted table

        Raises:
            :obj:`ValueError`: if a kinetic law or compartment is not contained in the local sqlite database
        """
        # group properties
        tsv = tsv.split('\n')
        law_properties = {}
        for row in csv.DictReader(tsv, delimiter='\t'):
            entry_id = int(float(row['EntryID']))
            if entry_id not in law_properties:
                law_properties[entry_id] = {
                    'KineticMechanismType': row['KineticMechanismType'],
                    'Tissue': row['Tissue'],
                    'Parameters': [],
                }

            parameter = {}

            # type
            parameter['type'] = row['parameter.type']

            if row['parameter.type'] == 'kcat':
                parameter['type_code'] = 25
            elif row['parameter.type'] == 'Vmax':
                parameter['type_code'] = 186
            elif row['parameter.type'] == 'Km':
                parameter['type_code'] = 27
            elif row['parameter.type'] == 'Ki':
                parameter['type_code'] = 261
            else:
                parameter['type_code'] = None

            # associatated species
            if row['parameter.associatedSpecies'] in ['', '-']:
                parameter['associatedSpecies'] = None
            else:
                parameter['associatedSpecies'] = row['parameter.associatedSpecies']

            # start value
            if row['parameter.startValue'] in ['', '-']:
                parameter['startValue'] = None
            else:
                parameter['startValue'] = float(row['parameter.startValue'])

            # end value
            if row['parameter.endValue'] in ['', '-']:
                parameter['endValue'] = None
            else:
                parameter['endValue'] = float(row['parameter.endValue'])

            # error
            if row['parameter.standardDeviation'] in ['', '-']:
                parameter['standardDeviation'] = None
            else:
                parameter['standardDeviation'] = float(
                    row['parameter.standardDeviation'])

            # units
            if row['parameter.unit'] in ['', '-']:
                parameter['unit'] = None
            else:
                parameter['unit'] = row['parameter.unit']

            law_properties[entry_id]['Parameters'].append(parameter)

        # update properties
        for id, properties in law_properties.items():
            # get kinetic law
            projection = {'kinlaw_id': 1, 'mechanism': 1, 'tissue': 1, 'parameters': 1}
            q = self.collection.find_one({'kinlaw_id': id}, projection=projection)
            if q == None:
                raise ValueError('No Kinetic Law with id {}'.format(id))
            law = q

            # mechanism
            print('KineticMechanismType: ' + properties['KineticMechanismType'])
            print('Tissue: ' + properties['Tissue'])
            if properties['KineticMechanismType'] == 'unknown':
                law['mechanism'] = None
            else:
                law['mechanism'] = properties['KineticMechanismType']

            # tissue
            if properties['Tissue'] in ['', '-']:
                law['tissue'] = None
                properties.pop('Tissue')
            else:
                law['tissue'] = properties['Tissue']
                properties.pop('Tissue')

            # parameter
            parameters = []
            for param_properties in properties['Parameters']:
                param = self.get_parameter_by_properties(law, param_properties)
                if param is None:
                    if param_properties['type'] != 'concentration':
                        warning = 'Unable to find parameter `{}:{}` for law {}'.format(
                            param_properties['type'], param_properties['associatedSpecies'], law['kinlaw_id'])
                        logging.warning(warning)
                    continue

                param['observed_value'] = param_properties['startValue']
                param['observed_error'] = param_properties['standardDeviation']
                param['observed_units'] = param_properties['unit']
                parameters.append(param)

            # updated
            law['modified'] = datetime.datetime.utcnow()
            self.collection.update_one({'kinlaw_id': id},
                                       {'$set': law, '$set': {'parameters': parameters} },
                                       upsert=False)

    def get_parameter_by_properties(self, kinetic_law, parameter_properties):
        """ Get the parameter of :obj:`kinetic_law` whose attribute values are 
                equal to that of :obj:`parameter_properties`
        Args:
            kinetic_law (:obj:`KineticLaw`): kinetic law to find parameter of
            parameter_properties (:obj:`dict`): properties of parameter to find

        Returns:
            :obj:`Parameter`: parameter with attribute values equal to values of :obj:`parameter_properties`
        """
        if parameter_properties['type'] == 'concentration':
            return parameter_properties

        # match observed name and compound
        def func(parameter):
            return parameter['observed_type'] == parameter_properties['type_code'] and \
                ((parameter['compound'] is None and parameter_properties['associatedSpecies'] is None) or
                 (parameter['compound'] is not None and parameter['compound']['name'] == parameter_properties['associatedSpecies']))
        parameters = list(filter(func, kinetic_law['parameters']))
        if len(parameters) == 1:
            return parameters[0]

        # match observed name
        def func(parameter):
            return parameter['observed_type'] == parameter_properties['type_code']
        parameters = list(filter(func, kinetic_law['parameters']))
        if len(parameters) == 1:
            return parameters[0]

        # match compound
        def func(parameter):
            return (parameter['compound'] is None and parameter_properties['associatedSpecies'] is None) or \
                (parameter['compound'] is not None and parameter['compound']
                 ['name'] == parameter_properties['associatedSpecies'])
        parameters = list(filter(func, kinetic_law['parameters']))
        if len(parameters) == 1:
            return parameters[0]

        # match value
        def func(parameter):
            return parameter['observed_value'] == parameter_properties['startValue']
        parameters = list(filter(func, kinetic_law['parameters']))
        if len(parameters) == 1:
            return parameters[0]

    def infer_compound_structures_from_names(self, compounds):
        """ Try to use PubChem to infer the structure of compounds from their names

        Notes: we don't try look up structures from their cross references because SABIO has already gathered
        all structures from their cross references to ChEBI, KEGG, and PubChem

        Args:
            compounds (:obj:`list` of :obj:`dict`): list of compounds
        """
        result = []
        for i_compound, compound in enumerate(compounds):
            # if self.verbose and (i_compound % 100 == 0):
            #     print('  Trying to infer the structure of compound {} of {}'.format(
            #         i_compound + 1, len(compounds)))

            if compound['name'] == 'Unknown':
                continue

            for i_try in range(self.PUBCHEM_MAX_TRIES):
                try:
                    p_compounds = pubchempy.get_compounds(
                        compound['name'], 'name')
                    break
                except pubchempy.PubChemHTTPError:
                    if i_try < self.PUBCHEM_MAX_TRIES - 1:
                        # sleep to avoid getting overloading PubChem server and then try again
                        time.sleep(self.PUBCHEM_TRY_DELAY)
                    else:
                        raise

            for p_compound in p_compounds:
                namespace = 'pubchem.compound'
                id = str(p_compound.cid)
                q = self.file_manager.search_dict_list(
                    compound['cross_references'], namespace)
                if len(q) == 0:
                    resource = {namespace: id}
                    compound['cross_references'].append(resource)

                structure = {'inchi': p_compound.inchi}
                norm = self.calc_inchi_formula_connectivity(structure)
                tmp = compound.setdefault('structures', [])
                tmp.append(norm)
            self.collection_compound.update_one({'_id': compound['_id']},
                                                {'$set': compound})
            result.append(compound)

        return result

    def calc_inchi_formula_connectivity(self, structure):
        """ Calculate a searchable structures

        * InChI format
        * Core InChI format

            * Formula layer (without hydrogen)
            * Connectivity layer
        """

        # if necessary, convert structure to InChI
        if 'inchi' in structure:
            _value_inchi = structure['inchi']
        else:
            try:
                _value_inchi = molecule_util.Molecule(
                    structure=structure.get('smiles',)).to_inchi() or None
            except ValueError:
                _value_inchi = None

        # calculate formula (without hydrogen) and connectivity
        if _value_inchi:
            _value_inchi_formula_connectivity = molecule_util.InchiMolecule(_value_inchi) \
                .get_formula_and_connectivity()

        result = {'_value_inchi': _value_inchi,
                  '_value_inchi_formula_connectivity': _value_inchi_formula_connectivity}

        return result

    def load_missing_enzyme_information_from_html(self, ids):
        """ Loading enzyme subunit information from html

        Args:
            ids (:obj:`list` of :obj:`int`): list of IDs of kinetic laws to download
        """
        query = {'$and': [{'kinlaw_id': {'$in': ids}},
                          {'enzyme._id': {'$exists': True}}]}
        projection = {'enzyme': 1, 'kinlaw_id': 1}
        kinetic_laws = self.collection.find(
            filter=query, projection=projection)
        total_count = self.collection.count_documents(query)

        for i_kinetic_law, kinetic_law in enumerate(kinetic_laws):
            if self.verbose and (i_kinetic_law % 100 == 0):
                print('  Loading enzyme information for {} of {} kinetic laws'.format(
                    i_kinetic_law + 1, total_count))

            response = requests.get(self.ENDPOINT_KINETIC_LAWS_PAGE, params={
                                    'kinlawid': kinetic_law['kinlaw_id'], 'newinterface': 'true'})
            response.raise_for_status()

            enzyme = kinetic_law['enzyme'][0]
            subunits = enzyme['subunits']
            for subunit in subunits:
                subunit['coefficient'] = None

            doc = bs4.BeautifulSoup(response.text, 'html.parser')
            td = doc.find('td', text='Modifier-Catalyst')
            tr = td.parent
            td = tr.find_all('td')[-1]
            inner_html = td.decode_contents(formatter='html').strip() + ' '
            if inner_html == '- ':
                continue
            try:
                subunit_coefficients = self.parse_complex_subunit_structure(
                    inner_html)
            except Exception as error:
                six.reraise(
                    ValueError,
                    ValueError('Subunit structure for kinetic law {} could not be parsed: {}\n\t{}'.format(
                        kinetic_law['kinlaw_id'], inner_html, str(error).replace('\n', '\n\t'))),
                    sys.exc_info()[2])

            enzyme['subunits'] = []
            for subunit_id, coefficient in subunit_coefficients.items():
                xref = {'uniprot': subunit_id}
                coeff = {'coefficient': coefficient}
                subunit = {**xref, **coeff}
                enzyme['subunits'].append(subunit)
            enzyme = self.calc_enzyme_molecular_weights([enzyme], 1)
            self.collection.update_one({'kinlaw_id': kinetic_law['kinlaw_id']},
                                       {'$set': {'enzyme': enzyme}})

    def parse_complex_subunit_structure(self, text):
        """ Parse the subunit structure of complex into a dictionary of subunit coefficients

        Args:
            text (:obj:`str`): subunit structure described with nested parentheses

        Returns:
            :obj:`dict` of :obj:`str`, :obj:`int`: dictionary of subunit coefficients
        """
        # try adding missing parentheses
        n_open = text.count('(')
        n_close = text.count(')')
        if n_open > n_close:
            text += ')' * (n_open - n_close)
        elif n_open < n_close:
            text = '(' * (n_close - n_open) + text

        # for convenenice, add parenthesis at the beginning and end of the string and before and after each subunit
        if text[0] != '(':
            text = '(' + text + ')'
        text = text.replace('<a ', '(<a ').replace('</a>', '</a>)')
        text = text.replace('open(', "open ").replace('able=1\')', "able=1\' ")
        # parse the nested subunit structure
        i = 0
        stack = [{'subunits': {}}]
        while i < len(text):
            if text[i] == '(':
                stack.append({'start': i, 'subunits': {}})
                i += 1
            elif text[i] == ')':
                tmp = stack.pop()
                start = tmp['start']
                end = i

                subunits = tmp['subunits']
                if not subunits:
                    protein_pattern = re.compile('<a href="#" onclick="window\.open \'http://sabiork\.h-its\.org/proteindetails\.jsp\?enzymeUniprotID=(.*?)\',\'\',\'width=600,height=500,scrollbars=1,resizable=1\' ">.*?</a>')
                    matches = protein_pattern.findall(text[start+1:end])
                    for match in matches:
                        subunits[match] = 1

                match = re.match(r'^\)\*(\d+)', text[i:])
                if match:
                    i += len(match.group(0))
                    coefficient = int(float(match.group(1)))
                else:
                    i += 1
                    coefficient = 1

                for id in subunits.keys():
                    if id not in stack[-1]['subunits']:
                        stack[-1]['subunits'][id] = 0
                    stack[-1]['subunits'][id] += subunits[id] * coefficient

            else:
                i += 1

        # check that all subunits were extracted
        protein_pattern = re.compile('<a href="#" onclick="window\.open \'http://sabiork\.h-its\.org/proteindetails\.jsp\?enzymeUniprotID=(.*?)\',\'\',\'width=600,height=500,scrollbars=1,resizable=1\' ">.*?</a>')
        matches = protein_pattern.findall(text)
        if len(set(matches)) != len(stack[0]['subunits'].keys()):
            raise ValueError(
                'Subunit structure could not be parsed: {}'.format(text))

        return stack[0]['subunits']

    def calc_enzyme_molecular_weights(self, enzymes, length):
        """ Calculate the molecular weight of each enzyme

        Args:
            enzymes (:obj:`list` of :obj:`dict`): list of enzymes
        Returns:
            enzymes (:obj:`list` of :obj:`dict`): list of enzymes
        """
        letters = Bio.Alphabet.IUPAC.IUPACProtein.letters
        mean_aa_mw = Bio.SeqUtils.molecular_weight(letters, seq_type='protein') / len(letters)

        results = []
        for i_enzyme, enzyme in enumerate(enzymes):
            # if self.verbose and i_enzyme % 100 == 0:
            #     print('    Calculating molecular weight of enzyme {} of {}'.format(i_enzyme + 1, length))
            enzyme_molecular_weight = 0
            for subunit in enzyme['subunits']:
                if 'uniprot' in subunit:
                    response = requests.get(
                        self.ENDPOINT_DOMAINS['uniprot'] + '/uniprot/?query={}&columns=id,sequence&format=tab'.format(subunit['uniprot']))
                    response.raise_for_status()
                    seqs = list(csv.DictReader(response.text.split('\n'), delimiter='\t'))
                    if seqs:
                        subunit['sequence'] = next((seq['Sequence'] for seq in seqs if seq['Entry'] == subunit['uniprot']), seqs[0]['Sequence'])
                        iupac_seq = re.sub(r'[^' + Bio.Alphabet.IUPAC.IUPACProtein.letters + r']', '', subunit['sequence'])
                        subunit['molecular_weight'] = \
                            + Bio.SeqUtils.molecular_weight(iupac_seq, seq_type='protein') \
                            + (len(subunit['sequence']) - len(iupac_seq)) * mean_aa_mw
                        enzyme_molecular_weight += (subunit['coefficient'] or float('nan')) * subunit['molecular_weight']
                    else:
                        subunit['sequence'] = None
                        subunit['molecular_weight'] = None
                        enzyme_molecular_weight = float('nan')

            if not enzyme_molecular_weight or math.isnan(enzyme_molecular_weight):
                enzyme['molecular_weight'] = None
            else:
                enzyme['molecular_weight'] = enzyme_molecular_weight
            results.append(enzyme)
        return results


    def add_inchi_hash(self):
        query = {}
        projection = {'products': 1, 'reactants':1, 'kinlaw_id': 1}
        cursor = self.collection.find({}, projection = projection)

        def get_inchi_structure(chem):
            '''Given subsrate or product subdocument from sabio_rk
               find the corresponding inchi
               Args:
                    chem (:obj: `dict`)
                Returns:
                    (:obj: `str`): inchi string
            '''
            try:
                return chem['structures'][0].get('_value_inchi', None)
            except IndexError:
                return chem['structures'].append({})

        def iter_rxnp_subdoc(rxnp):
            '''Given a reactant or product array from sabio_rk
                append fields of hashed inchi
                Args:
                    rxnp (:obj: `list` of :obj: `dict`)
            '''
            for i in range(len(rxnp)):
                substrate_inchi = get_inchi_structure(rxnp[i])
                try:
                    hashed_inchi = hashlib.sha224(substrate_inchi.encode()).hexdigest()
                    rxnp[i]['structures'][0]['hashed_inchi'] = hashed_inchi
                except AttributeError:
                    rxnp[i]['structures'][0]['hashed_inchi'] = None

            return rxnp

        j = 0
        for doc in cursor:
            if j > self.max_entries:
                break
            if self.verbose == True and j % 100 == 0:
                print(j)
            substrates = doc['reactants']
            products = doc['products']
            new_subsrates = iter_rxnp_subdoc(substrates)
            new_products = iter_rxnp_subdoc(products)

            doc['reactants'] = new_subsrates
            doc['products'] = new_products

            self.collection.update_one({'kinlaw_id': doc['kinlaw_id']},
                            {'$set': {'reactants': doc['reactants'],
                                      'products': doc['products']} })
            j += 1

def main():
        db = 'datanator'
        username = datanator.config.core.get_config()[
            'datanator']['mongodb']['user']
        password = datanator.config.core.get_config(
        )['datanator']['mongodb']['password']
        MongoDB = datanator.config.core.get_config(
        )['datanator']['mongodb']['server']
        port = datanator.config.core.get_config(
        )['datanator']['mongodb']['port']
        replSet = datanator.config.core.get_config(
        )['datanator']['mongodb']['replSet']
        manager = SabioRk(MongoDB=MongoDB,  db=db,
                                 verbose=True, username=username,
                                 password=password)
        manager.load_content()

if __name__ == '__main__':
    main()
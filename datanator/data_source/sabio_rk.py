import datanator.config.core
from datanator.util import mongo_util
import six
import requests
from xml import etree
import libsbml
import re
import datetime

class SabioRk:

    def __init__(self, cache_dirname=None, MongoDB=None, replicaSet=None, db=None,
                 verbose=False, max_entries=float('inf'), username=None,
                 password=None, authSource='admin', webservice_batch_size=100):
        self.cache_dirname = cache_dirname
        self.MongoDB = MongoDB
        self.replicaSet = replicaSet
        self.db = db
        self.verbose = verbose
        self.max_entries = max_entries
        self.client, self.db_obj, self.collection = mongo_util.MongoUtil(
            MongoDB=MongoDB, db=db, username=username, password=password, authSource=authSource).con_db('sabio_rk')
        ENDPOINT_DOMAINS = {
            'sabio_rk': 'http://sabiork.h-its.org',
            'uniprot': 'http://www.uniprot.org',
        }
        self.ENDPOINT_KINETIC_LAWS_SEARCH = ENDPOINT_DOMAINS['sabio_rk'] + \
            '/sabioRestWebServices/searchKineticLaws/entryIDs'
        self.ENDPOINT_WEBSERVICE = ENDPOINT_DOMAINS['sabio_rk'] + \
            '/sabioRestWebServices/kineticLaws'
        self.ENDPOINT_EXCEL_EXPORT = ENDPOINT_DOMAINS['sabio_rk'] + \
            '/entry/exportToExcelCustomizable'
        self.ENDPOINT_COMPOUNDS_PAGE = ENDPOINT_DOMAINS['sabio_rk'] + \
            '/compdetails.jsp'
        self.ENDPOINT_KINETIC_LAWS_PAGE = ENDPOINT_DOMAINS['sabio_rk'] + \
            '/kindatadirectiframe.jsp'
        self.PUBCHEM_MAX_TRIES = 10
        self.PUBCHEM_TRY_DELAY = 0.25
        self.webservice_batch_size = webservice_batch_size

    def load_kinetic_law_ids(self):
        """ Download the IDs of all of the kinetic laws stored in SABIO-RK

        Returns:
            :obj:`list` of :obj:`int`: list of kinetic law IDs

        Raises:
            :obj:`Error`: if an HTTP request fails or the expected number of kinetic laws is not returned
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

        for i_batch in range(int(math.ceil(float(len(ids)) / batch_size))):
            if self.verbose and (i_batch % max(1, 100. / batch_size) == 0):
                print('  Downloading kinetic laws {}-{} of {} in SBML format'.format(
                    i_batch * batch_size + 1,
                    min(len(ids), i_batch * batch_size + max(100, batch_size)),
                    len(ids)))

            batch_ids = ids[i_batch * batch_size:min((i_batch + 1) * batch_size, len(ids))]
            response = session.get(self.ENDPOINT_WEBSERVICE, params={
                'kinlawids': ','.join(str(id) for id in batch_ids),
            })

            response.raise_for_status()
            if not response.text:
                raise Exception('Unable to download kinetic laws with ids {}'.format(', '.join([str(id) for id in batch_ids])))

            self.create_kinetic_laws_from_sbml(batch_ids, response.content if six.PY2 else response.text)

            if self.commit_intermediate_results:
                self.session.commit()

        # print warning with list of unidentified ids
        loaded_ids = [l.id for l in self.session.query(KineticLaw).order_by(KineticLaw.id)]

        not_loaded_ids = list(set(ids).difference(loaded_ids))
        if not_loaded_ids:
            not_loaded_ids.sort()
            warnings.warn('Several kinetic laws were not found:\n  {}'.format(
                '\n  '.join([str(id) for id in not_loaded_ids])), data_source.DataSourceWarning)

    def create_kinetic_laws_from_sbml(self, ids, sbml):
        """ Add kinetic laws defined in an SBML file to the local mongodb database

        Args:
            ids (:obj:`list` of :obj:`int`): list kinetic law IDs
            sbml (:obj:`str`): SBML representation of one or more kinetic laws

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
                eq = libsbml.formulaToL3String(math_sbml.getChild(math_sbml.getNumChildren() - 1))
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

        # compartments
        compartments_sbml = model.getListOfCompartments()
        compartments = []
        for i_compartment in range(compartments_sbml.size()):
            compartment_sbml = compartments_sbml.get(i_compartment)
            compartments.append(self.get_compartment_from_sbml(compartment_sbml))

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
            raise ValueError('{} reactions {} is different from the expected {}'.format(reaction_sbml.size(), len(ids)))
        kinetic_laws = []
        for i_reaction, id in enumerate(ids):
            reaction_sbml = reactions_sbml.get(i_reaction)
            kinetic_law = self.create_kinetic_law_from_sbml(
                id, reaction_sbml, specie_properties, functions, units)
            kinetic_laws.append(kinetic_law)

        return (kinetic_laws, species, compartments)

    def create_kinetic_law_from_sbml(self, id, sbml, specie_properties, functions, units):
        """ Add a kinetic law to the mongodb database

        Args:
            id (:obj:`int`): identifier
            sbml (:obj:`libsbml.KineticLaw`): SBML-representation of a reaction
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

        # stop if kinetic law entry is empty
        if not law.getMetaId():
            return None

        # ID
        annotated_id = next((int(float(x_ref['id'])) for x_ref in x_refs if x_ref['namespace'] == 'sabiork.kineticrecord'), None)
        if annotated_id is not None and annotated_id != id:
            raise ValueError('Annotated ID {} is different from expected ID {}'.format(annotated_id, id))
        id = annotated_id

        """ participants """
        kinetic_law['reactants'] = []
        reactants = sbml.getListOfReactants()
        for i_part in range(reactants.size()):
            part_sbml = reactants.get(i_part)
            compound, compartment = self.get_specie_reference_from_sbml(part_sbml.getSpecies())
            part = {
                'compound':compound,
                'compartment':compartment,
                'coefficient':part_sbml.getStoichiometry()}
            kinetic_law['reactants'].append(part)

        kinetic_law.products[:] = []
        products = sbml.getListOfProducts()
        for i_part in range(products.size()):
            part_sbml = products.get(i_part)
            compound, compartment = self.get_specie_reference_from_sbml(part_sbml.getSpecies())
            part = {
                'compound':compound,
                'compartment':compartment,
                'coefficient':part_sbml.getStoichiometry()}
            kinetic_law['products'].append(part)

        """ cross references """
        # Note: these are stored KineticLaws rather than under Reactions because this seems to how SABIO-RK stores this information.
        # For example, kinetic laws 16016 and 28003 are associated with reaction 9930, but they have different EC numbers 1.1.1.52 and
        # 1.1.1.50, respectively.
        kinetic_law.cross_references = list(filter(lambda x_ref: x_ref.namespace not in ['taxonomy'], reaction_x_refs))

        # rate_law
        kinetic_law.equation = functions[law.getMetaId()[5:]]

        # parameters
        kinetic_law.parameters = []
        params = law.getListOfLocalParameters()
        for i_param in range(params.size()):
            param = params.get(i_param)

            match = re.match(r'^(.*?)_((SPC|ENZ)_([0-9]+)_(.*?))$', param.getId(), re.IGNORECASE)
            if match:
                observed_name = match.group(1)
                species, compartment = self.get_specie_reference_from_sbml(match.group(2))
                if isinstance(species, Compound):
                    compound = species
                    enzyme = None
                else:
                    compound = None
                    enzyme = species
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

            parameter = Parameter(
                compound=compound,
                enzyme=enzyme,
                compartment=compartment,
                observed_name=observed_name,
                observed_type=observed_type,
                observed_value=observed_value,
                observed_units=observed_units,
                modified=datetime.datetime.utcnow(),
            )
            self.session.add(parameter)
            kinetic_law.parameters.append(parameter)

        # modifiers to kinetic law
        kinetic_law.modifiers[:] = []
        modifiers = sbml.getListOfModifiers()
        for i_modifier in range(modifiers.size()):
            modifier = modifiers.get(i_modifier)
            modifier_id = modifier.getSpecies()
            specie, compartment = self.get_specie_reference_from_sbml(modifier_id)
            type = specie_properties[modifier.getSpecies()]['modifier_type']
            if modifier.getSpecies()[0:3] == 'SPC':
                part = ReactionParticipant(
                    compound=specie,
                    compartment=compartment,
                    type=type,
                )
                self.session.add(part)
                kinetic_law.modifiers.append(part)
            elif modifier_id[0:3] == 'ENZ':
                kinetic_law.enzyme, kinetic_law.enzyme_compartment = self.get_specie_reference_from_sbml(modifier_id)
                kinetic_law.enzyme_type = specie_properties[modifier.getSpecies()]['modifier_type']
                kinetic_law.taxon_wildtype = specie_properties[modifier_id]['is_wildtype']
                kinetic_law.taxon_variant = specie_properties[modifier_id]['variant']

        # taxon
        kinetic_law.taxon = next((int(float(x_ref.id)) for x_ref in reaction_x_refs if x_ref.namespace == 'taxonomy'), None)

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
                raise ValueError('Unsupported temperature units: {}'.format(temperature_units))
            kinetic_law.temperature = temperature

        # pH
        if conditions.hasChild('pH'):
            ph = conditions \
                .getChild('pH') \
                .getChild('startValuepH') \
                .getChild(0) \
                .getCharacters()
            kinetic_law.ph = float(ph)

        # media
        if conditions.hasChild('buffer'):
            media = conditions \
                .getChild('buffer') \
                .getChild(0) \
                .getCharacters() \
                .strip()
            if six.PY2:
                media = unicode(media.decode('utf-8'))
            kinetic_law.media = media

        """ references """
        kinetic_law.references = list(filter(lambda x_ref: x_ref.namespace != 'sabiork.kineticrecord', x_refs))

        """ updated """
        kinetic_law.modified = datetime.datetime.utcnow()

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

    def get_specie_reference_from_sbml(self, specie_id):
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
            name = sbml.getName()
            properties = {'modifier_type': modifier_type}
            specie = {'_id': id, 'name': name}
        elif type == 'ENZ':
            name, is_wildtype, variant = self.parse_enzyme_name(sbml.getName())
            if six.PY2:
                variant = unicode(variant.decode('utf-8'))
            properties = {'is_wildtype': is_wildtype, 'variant': variant, 'modifier_type': modifier_type}

            specie = {'_id': id, 'molecular_weight': None, 'name': name}
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
        elif type == 'ENZ':
            name, is_wildtype, variant = self.parse_enzyme_name(sbml.getName())
            if six.PY2:
                variant = unicode(variant.decode('utf-8'))
            properties = {'is_wildtype': is_wildtype, 'variant': variant, 'modifier_type': modifier_type}

            specie = {'_id': id, 'molecular_weight': None, 'name': name}
        else:
            raise ValueError('Unsupported species type: {}'.format(type))

        # cross references
        cross_references = self.create_cross_references_from_sbml(sbml)
        if type == 'SPC':
            specie['cross_references'] = cross_references
        elif type == 'ENZ':
            specie['subunits'] = []
            specie['cross_references'] = []
            for cross_reference in cross_references:
                if cross_reference['namespace'] == 'uniprot':
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
        match = re.match(r'^(.*?)\(Enzyme\) (wildtype|mutant),?(.*?)$', sbml, re.IGNORECASE)
        if match:
            name = match.group(1)
            is_wildtype = match.group(2).lower() == 'wildtype'
            variant = match.group(3).strip()
            return (name, is_wildtype, variant)

        match = re.match(r'^Enzyme (wildtype|mutant),?( (.*?))*$', sbml, re.IGNORECASE)
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

        attr = libsbml.XMLTriple('resource', 'http://www.w3.org/1999/02/22-rdf-syntax-ns#', 'rdf')

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
                    resource = {'namespace': namespace, 'id': id}

                if resource not in x_refs:
                    x_refs.append(resource)

        return x_refs
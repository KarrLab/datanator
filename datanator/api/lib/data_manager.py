from datanator.core import data_model, common_schema


class BaseManager(object):

    def __init__(self, cache_dirname=None):
        self.data_source = common_schema.CommonSchema(cache_dirname=cache_dirname)

    def metadata_dump(self, component):
        """ Calculate a consensus statistical representation of the one or more observed values

        Args:
            component (:obj:`models.Observation`): model component dump data for

        Returns:
            :obj:`list` of :obj:`data_model.ObservedResultMetadata`: data model metadata object
        """
        genetics = None
        environment=None
        cross_references= []
        method=None
        synonym= []
        meta = component._metadata


        taxon = meta.taxon[0].name if meta.taxon else None
        variation = meta.cell_line[0].name if meta.cell_line else None
        genetics = data_model.Genetics(taxon=taxon, variation=variation)

        temperature = meta.conditions[0].temperature if meta.conditions else None
        ph = meta.conditions[0].ph if meta.conditions else None
        media = meta.conditions[0].media if meta.conditions else None
        growth_status = meta.conditions[0].growth_status if meta.conditions else None
        growth_system = meta.conditions[0].growth_system if meta.conditions else None
        environment = data_model.Environment(temperature=temperature, ph=ph, media=media, growth_status=growth_status, growth_system=growth_system)

        name = meta.method[0].name if meta.method else None
        description =  meta.method[0].comments if meta.method else None
        performer = meta.method[0].performer if meta.method else None
        hardware =  meta.method[0].hardware if meta.method else None
        software =  meta.method[0].software if meta.method else None
        method = data_model.Method(name=name, description=description, performer=performer, hardware=hardware, software=software)


        if meta.resource:
            for item in meta.resource:
                cross_references.append(data_model.Resource(namespace=item.namespace, id=item._id))
        else:
            cross_references.append(data_model.Resource(namespace=None, id=None))


        if meta.synonym:
            for item in meta.synonym:
                synonym.append(data_model.Synonym(name=item.name))
        else:
            synonym.append(data_model.Synonym(name=None))

        metadata_result = data_model.ObservedResultMetadata(genetics = genetics, environment=environment, cross_references=cross_references, method=method, synonym=synonym)

        return metadata_result

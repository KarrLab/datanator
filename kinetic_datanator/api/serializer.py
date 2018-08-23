from kinetic_datanator import ma
from kinetic_datanator.core import models


class CellCompartmentSerializer(ma.ModelSchema):

    class Meta:
        exclude = ["search_vector", '_metadata' , 'reaction']
        model = models.CellCompartment

class TaxonSerializer(ma.ModelSchema):

    class Meta:
        exclude = ["search_vector", '_metadata', '_experimentmetadata']
        model = models.Taxon

class MetadataSerializer(ma.ModelSchema):
    cell_compartment = ma.Nested(CellCompartmentSerializer, many=True)
    taxon = ma.Nested(TaxonSerializer, many=True)

    class Meta:
        fields = ['cell_compartment', "taxon"]
        model = models.Metadata

class StructureSerializer(ma.ModelSchema):

    class Meta:
        exclude = ['compound', 'name', 'type']
        model = models.Structure

class CompoundSerializer(ma.ModelSchema):

    _metadata = ma.Nested(MetadataSerializer)
    structure = ma.Nested(StructureSerializer)

    class Meta:
        exclude = ["search_vector", "parameter", "reaction", "_is_name_ambiguous", "concentration"]
        model = models.Compound

class ProteinComplexSerializer(ma.ModelSchema):
    _metadata = ma.Nested(MetadataSerializer)
    class Meta:
        exclude = ["search_vector"]
        model = models.ProteinComplex

class ProteinSubunitSerializer(ma.ModelSchema):
    _metadata = ma.Nested(MetadataSerializer)
    class Meta:
        exclude = ["search_vector"]
        model = models.ProteinSubunit




### ----------------------------------------------------- ######

class ResourceSerializer(ma.Schema):

    class Meta:
        fields = ['namespace', 'id']

class SpecieSerializer(ma.Schema):

    class Meta:
        fields = ['structure']


class EntityInteractionOrPropertySerializer(ma.Schema):

    cross_references = ma.Nested(ResourceSerializer, many=True)

    class Meta:
        fields = ['id', 'name', 'cross_references']

class InteractionSerializer(EntityInteractionOrPropertySerializer):


    class Meta(EntityInteractionOrPropertySerializer):
        # Fields to expose
        fields = EntityInteractionOrPropertySerializer.Meta.fields + ['position', 'score', 'confidence']


class ReactionParticipantSerializer(ma.Schema):

    specie = ma.Nested(SpecieSerializer)

    class Meta:
        fields = ['specie', 'coefficient', 'order']

class ReactionSerializer(InteractionSerializer):

    participants = ma.Nested(ReactionParticipantSerializer, many=True)

    class Meta(InteractionSerializer):
        fields = InteractionSerializer.Meta.fields+['kinetic_law_id', 'reversible', 'participants']


class ObservableSerializer(ma.Schema):

    interaction = ma.Nested(InteractionSerializer)
    specie = ma.Nested(SpecieSerializer)
    # compartment = ma.Nested(CompartmentSerializer)

    class Meta:
        # Fields to expose
        fields = ['property', 'interaction', 'specie']

class GeneticsSerializer(ma.Schema):

    class Meta:
        fields = ['taxon', 'variation']

class ObservedResultMetadataSerializer(ma.Schema):

    genetics = ma.Nested(GeneticsSerializer)

    class Meta:
        fields = ['genetics']


class ObservedResultSerializer(ma.Schema):

    metadata = ma.Nested(ObservedResultMetadataSerializer)

    class Meta:
        fields = ['metadata']

class ObservedValueSerializer(ObservedResultSerializer):

    observable = ma.Nested(ObservableSerializer)

    class Meta(ObservedResultSerializer):
        # Fields to expose
        fields = ObservedResultSerializer.Meta.fields + ['observable', 'value', 'error', 'units']

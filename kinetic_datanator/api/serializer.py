from kinetic_datanator.app import ma
from kinetic_datanator.core import models


class CompoundSerializer(ma.ModelSchema):
    class Meta:
        exclude = ["search_vector", "parameter", "reaction", "_is_name_ambiguous", "_metadata", "concentration"]
        model = models.Compound

class ProteinComplexSerializer(ma.ModelSchema):
    class Meta:
        exclude = ["search_vector"]
        model = models.ProteinComplex

class ProteinSubunitSerializer(ma.ModelSchema):
    class Meta:
        exclude = ["search_vector"]
        model = models.ProteinSubunit




### ----------------------------------------------------- ######




class ObservableSerializer(ma.Schema):

    # interaction = ma.Nested(InteractionSerializer)
    # specie = ma.Nested(SpecieSerializer)
    # compartment = ma.Nested(CompartmentSerializer)

    class Meta:
        # Fields to expose
        fields = ['property']

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

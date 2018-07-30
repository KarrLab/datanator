# """"
# TODO: Turn into marshmallow seraializers
# """"
#
#
# class SerializeClass(object):
#     """
#     Mixin class provides data object special functions such as serialization
#     """
#     ##TODO: Make more succinct and use a stack instead of recursion
#     def serialize_relationships(self):
#         result_json = {}
#         for relation in self.__mapper__.relationships.keys():
#             if relation == '_metadata':
#                 continue
#             objs = getattr(self, relation)
#
#             if objs == None:
#                 result_json[relation] = None
#             elif isinstance(objs, list):
#                 result_json[relation] ={i: {c.name: getattr(meta, c.name) for c in meta.__table__.columns} for i,meta in enumerate(objs)}
#             else:
#                 result_json[relation] = {c.name: getattr(objs, c.name) for c in objs.__table__.columns}
#
#         return(result_json)
#
#     def serialize_metadata(self, obj):
#         result_json = {}
#         for relation in obj.__mapper__.relationships.keys():
#             if relation == 'observation':
#                 continue
#             objs = getattr(obj, relation)
#             result_json[relation] ={i: {c.name: getattr(meta, c.name) for c in meta.__table__.columns} for i,meta in enumerate(objs)} if objs else None
#         return(result_json)
#
#     def serialize(self, metadata=False, relationships=False):
#         """
#         Attributes:
#             metadata (:obj:`bool`): Switch for gathering metadata serialization data
#             relationships (:obj:`bool`): Switch for gathering relationship serialization data
#
#         """
#
#         json = {c.name: getattr(self, c.name) for c in self.__table__.columns}
#         if relationships:
#             json['relationships'] = self.serialize_relationships()
#         if metadata:
#             json['metadata'] = self.serialize_metadata(self._metadata)
#         return json

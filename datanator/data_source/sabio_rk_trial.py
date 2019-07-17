from suds.client import Client
from suds.transport import TransportError

class SabioRk:

    def suds_trial(self):
        url = "http://sabio.villa-bosch.de/sabiork?wsdl"
        try:
            client = Client(url)
            print(client.service.getPathwayNames("ABC"))
        except TransportError as e:
            return str(e)

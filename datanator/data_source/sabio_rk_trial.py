from suds.client import Client


class SabioRk:

    def suds_trial(self):
        url = "http://sabiork.h-its.org/sabiork?wsdl"
        client = Client(url)
        print(client.service.getPathwayNames("ABC"))

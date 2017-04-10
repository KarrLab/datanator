import requests
from jxmlease import parse
import urllib
from xml.etree import ElementTree

#url = 'http://sabio.villa-bosch.de/compdetails.jsp?cid=35'

file = open("InchiToCompound2.txt", "w")
file.write("")
errorFile = open("datasetErrors2.txt", "w")
errorFile.write("")
errorFile = open("datasetErrors2.txt", "a")
file = open("InchiToCompound2.txt", "a")
i = 0
while i < 30000:
	try:

		answer = requests.get("http://sabio.villa-bosch.de/compdetails.jsp?cid={}".format(i)).text
		if len(answer) > 1600:
		#print(answer)
		#xmldata = '<root>' + answer.text + '</root>'
		#print(xmldata)
		#tree = ElementTree.fromstring(xmldata)
		#print(answer.text)
			start =  answer.index("common")
			string =  answer[start+12:start+400]
			name = string[:string.find("<")]


		#print(parse(tree))
			try:
				start =  answer.index("Inchi")
				string =  answer[start:start+400]
				inchi = string[:string.find("<")]
			except:
				inchi = "No Inchi Found"

			file.write("{}".format(i) + " - " +inchi + " - " + name + "\n")
	except:
		errorFile.write("{}".format(i) + "\n")
	i = i+1





"""
answer = requests.get("http://sabio.villa-bosch.de/compdetails.jsp?cid=35")
print(answer.text)
parsed = parse(answer.text)
print(parsed)
"""
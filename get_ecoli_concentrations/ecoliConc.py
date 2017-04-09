import json
import requests 
from jxmlease import Parser

myparser = Parser()



def trimInchi(inchi):
	if "/h" in inchi:
		end = inchi.index("/h")
		inchi = inchi[:end]
	return inchi


#def getMetabolomicInfo(inchi)

inchi = """InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/t4-,6-,7-,10-/m1/s1"""

with open('ecmdb.json') as data_file:    
    data = json.load(data_file)



for entry in data:
	if trimInchi(entry["moldb_inchi"]) == trimInchi(inchi):
		print entry['name']
		compId = entry["m2m_id"]
		print compId

		response = requests.get("""http://ecmdb.ca/compounds/{}.xml""".format(compId)).text

		#print response
		data = myparser(response)
		print data["compound"]["concentrations"]

#print data[0]
#print data[1]["moldb_inchi"]

"""
i = 0
for thing in data:
	answer = thing["moldb_inchi"]
	if answer[:5] != "InChI":
		print answer
	i = i+1
print i
"""
""" 
:Author: Yosef Roth <yosefdroth@gmail.com>
:Date: 2017-04-06
:Copyright: 2017, Karr Lab
:License: MIT
"""

import requests
import logging
from kinetic_datanator.util import compound_util
from requests.packages.urllib3.exceptions import InsecureRequestWarning
requests.packages.urllib3.disable_warnings(InsecureRequestWarning)


def easy_find_ec_number(substrates, products):
    #Unlike this method, get_ec_number requires a precise matchup of substrates to products
    #therefore, this method - easy_find_ec_number - tries all the combinations in get_ec_number
    #its basically an add-on to get_ec_number that makes using it much easier

    logging.info("ec_number_finder: Attempting to find EC number")

    ec_num = ""

    #first check to make sure that each compound has a smiles or an inchi
    #if they do not, return an empty string
    every_compound_has_structural_info = True

    for compound in substrates:
        if len(compound.inchi_smiles)==0:
            every_compound_has_structural_info = False
    for compound in products:
        if len(compound.inchi_smiles)==0:
            every_compound_has_structural_info = False
    #print("Every Compound Has Stuctural Info: {}".format(every_compound_has_structural_info))

    if every_compound_has_structural_info==False:
        logging.info("ec_number_finder: Not every compound has structural information, and therefore no EC number was found")

    if every_compound_has_structural_info==True:
        logging.info("ec_number_finder: Every compound has structural information, so an EC number should be found")
        sub_inchi_smiles = []
        prod_inchi_smiles = []

        for compound in substrates:
            sub_inchi_smiles.append(compound.inchi_smiles)
            #print("Substrate: {}".format(compound.id))
            logging.info("ec_number_finder: Substrate: {}".format(compound.inchi_smiles))
        for compound in products:
            prod_inchi_smiles.append(compound.inchi_smiles)
            logging.info("ec_number_finder: Product: {}".format(compound.inchi_smiles))
            #print("Product: {}".format(compound.id))


        i = 0
        while i<len(sub_inchi_smiles) and len(ec_num) == 0:
            ec_num = get_ec_number(sub_inchi_smiles, prod_inchi_smiles)
            sub_inchi_smiles = sub_inchi_smiles[-1:] + sub_inchi_smiles[:-1]
            i += 1

        if len(products)>len(substrates) and len(ec_num) == 0:
            prod_inchi_smiles = prod_inchi_smiles[-1:] + prod_inchi_smiles[:-1]
            ec_num = get_ec_number(sub_inchi_smiles, prod_inchi_smiles)
            i += 1

        if len(ec_num)==0:
            logging.error("ec_number_finder: No EC Found, but all structural information is present. This is troubling. Look into this")
        else:
            logging.info("ec_number_finder: EC Number Found - {}".format(ec_num))
    #print("ec_num: {}".format(ec_num))
    return ec_num

def get_ec_number(sub_array, prod_array):
    #input an array of substrates and an array of products
    #outputs the predicted EC number based on KEGG's ezyme 
    #the subsrates must be precisely matched with  their corresponding product

    s1 = 0
    s2 = 0
    s3 = 0

    p1 = 0
    p2 = 0
    p3 = 0
    
    if len(sub_array) > 0:
        s1 = compound_util.Compound(sub_array[0]).to_mol()
    if len(prod_array) >0:
        p1 = compound_util.Compound(prod_array[0]).to_mol()
    if len(sub_array) > 1:
        s2 = compound_util.Compound(sub_array[1]).to_mol()
    if len(prod_array) > 1:
        p2 = compound_util.Compound(prod_array[1]).to_mol()
    if len(sub_array) > 2:
        s3 = compound_util.Compound(sub_array[2]).to_mol()
    if len(prod_array) > 2:
        p3 = compound_util.Compound(prod_array[2]).to_mol()

    headers = {
        'Origin': 'http://www.genome.jp',
        'Accept-Encoding': 'gzip, deflate',
        'Accept-Language': 'en-US,en;q=0.8',
        'Upgrade-Insecure-Requests': '1',
        'User-Agent': 'Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/55.0.2883.87 Safari/537.36',
        'Content-Type': 'multipart/form-data; boundary=----WebKitFormBoundarykniU1HTxfbLnE5QK',
        'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
        'Cache-Control': 'max-age=0',
        'Referer': 'http://www.genome.jp/tools-bin/predict_reaction',
        'Connection': 'keep-alive',
    }

    data = """$------WebKitFormBoundarykniU1HTxfbLnE5QK\r\nContent-Disposition: form-data; name="QUERY_MODE"\r\n\r\nMULTI\r\n------WebKitFormBoundarykniU1HTxfbLnE5QK\r\nContent-Disposition: form-data; name="S_ENTRY"\r\n\r\n\r\n------WebKitFormBoundarykniU1HTxfbLnE5QK\r\nContent-Disposition: form-data; name="P_ENTRY"\r\n\r\n\r\n------WebKitFormBoundarykniU1HTxfbLnE5QK\r\nContent-Disposition: form-data; name="S_MOLFILE"; filename=""\r\nContent-Type: application/octet-stream\r\n\r\n\r\n------WebKitFormBoundarykniU1HTxfbLnE5QK\r\nContent-Disposition: form-data; name="P_MOLFILE"; filename=""\r\nContent-Type: application/octet-stream\r\n\r\n\r\n------WebKitFormBoundarykniU1HTxfbLnE5QK\r\nContent-Disposition: form-data; name="S_MOLTEXT"\r\n\r\n{}\r\n------WebKitFormBoundarykniU1HTxfbLnE5QK\r\nContent-Disposition: form-data; name="P_MOLTEXT"\r\n\r\n{}\r\n------WebKitFormBoundarykniU1HTxfbLnE5QK\r\nContent-Disposition: form-data; name="S_ENTRY"\r\n\r\n\r\n------WebKitFormBoundarykniU1HTxfbLnE5QK\r\nContent-Disposition: form-data; name="P_ENTRY"\r\n\r\n\r\n------WebKitFormBoundarykniU1HTxfbLnE5QK\r\nContent-Disposition: form-data; name="S_MOLFILE"; filename=""\r\nContent-Type: application/octet-stream\r\n\r\n\r\n------WebKitFormBoundarykniU1HTxfbLnE5QK\r\nContent-Disposition: form-data; name="P_MOLFILE"; filename=""\r\nContent-Type: application/octet-stream\r\n\r\n\r\n------WebKitFormBoundarykniU1HTxfbLnE5QK\r\nContent-Disposition: form-data; name="S_MOLTEXT"\r\n\r\n{}\r\n------WebKitFormBoundarykniU1HTxfbLnE5QK\r\nContent-Disposition: form-data; name="P_MOLTEXT"\r\n\r\n{}\r\n------WebKitFormBoundarykniU1HTxfbLnE5QK\r\nContent-Disposition: form-data; name="S_ENTRY"\r\n\r\n\r\n------WebKitFormBoundarykniU1HTxfbLnE5QK\r\nContent-Disposition: form-data; name="P_ENTRY"\r\n\r\n\r\n------WebKitFormBoundarykniU1HTxfbLnE5QK\r\nContent-Disposition: form-data; name="S_MOLFILE"; filename=""\r\nContent-Type: application/octet-stream\r\n\r\n\r\n------WebKitFormBoundarykniU1HTxfbLnE5QK\r\nContent-Disposition: form-data; name="P_MOLFILE"; filename=""\r\nContent-Type: application/octet-stream\r\n\r\n{}\r\n------WebKitFormBoundarykniU1HTxfbLnE5QK\r\nContent-Disposition: form-data; name="S_MOLTEXT"\r\n\r\n{}\r\n------WebKitFormBoundarykniU1HTxfbLnE5QK\r\nContent-Disposition: form-data; name="P_MOLTEXT"\r\n\r\n\r\n------WebKitFormBoundarykniU1HTxfbLnE5QK--\r\n""".format(s1,p1,s2,p2,s3,p3)
    response = requests.post('http://www.genome.jp/tools-bin/predict_view', headers=headers, data=data)
    

    #now we get the EC number
    string = response.text
    string = string[string.find("""<input type=hidden name=file value="""):]
    string = string[string.find("e-zyme"): string.find(""">""")-1]
    ezyme = string


    data = {'elup':'c',
    'file': ezyme,
    'name':':C00002:C00020:C00144:C00044::',
    'pview':'gif'}
    url2 = "http://www.genome.jp/tools-bin/e-zyme2/result.cgi"
    request = requests.post(url2, data=data)
    request.raise_for_status()
    #print(request.text)

    string = request.text
    string = string[string.find("""<td align=center><a href="http://www.genome.jp/kegg-bin/get_htext?htext=ko01000.keg&query="""):]
    return string[90:95]

def format_ec_for_sabio(ec_number):
    string = ""
    if len(ec_number) > 0:
        #string = "ECNumber: ("
        i = 0
        while i < 101:
            string = string + " OR " + '"{}.{}"'.format(ec_number, i)
            #print(string)
            i += 1
        #remove first "or"
        string = string[4:]
        string = "ECNumber: (" + string + ")"

    #print(string)
    return string

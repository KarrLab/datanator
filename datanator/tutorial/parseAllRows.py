import csv
#read .tab file 
in_stream = open('Computed-Gram_negative_without_outer_membrane-PSORTdb-3.00.tab','r') 
lines_reader = csv.reader(in_stream,delimiter='\t')

#row1 is the first row of the file.  Each column is a separate element in the list.
#next() removes the first line from lines_reader
row1 = next(lines_reader)

#jsonText stores text of JSON file that will be created
jsonText="{\n\t"

#for each line in the file, create a JSON file.
for line in lines_reader:
    for i in range(len(row1)):
        jsonText+='"'+row1[i]+'": '  #variable name in JSON file
        #if the value is a decimal/number
        if (line[i].replace(".","").isdecimal()):
            jsonText+=line[i]+",\n\t" #put value without quotation marks into the JSON file
        else:
            jsonText+='"'+line[i]+'",\n\t' #put value between quotation marks into the JSON file
    
    jsonText+="}"

    #create JSON file
    out_stream=open(line[0][line[0].find("W"):(line[0][line[0].find("W"):]).find("|")+line[0].find("W")]+'.json','w')
    out_stream.write(jsonText)
    out_stream.close()

    
        

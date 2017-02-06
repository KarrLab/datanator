#import openpyxl
from openpyxl import Workbook


#takes in FormattedData as an input. Saves a file 

def createExcelSheet(FormattedDataList, species):
	topRow = [['Name', 'Sabio Reaction IDs', 
	'Vmax Median Value', 'Vmax Median Entry', 'Vmax Min Entry', 'Vmax Max Entry', 'Proximity']]#, 'Sabio Entry IDs', 'Vmax Values',
	#'Km Median Value', 'Km Median Entry', 'Km Min Enry', 'Km Max Entry', 'Proximity', 'Sabio Entry IDs', 'Km Values']

	rows = []
	for formattedData in FormattedDataList:
		row = []
		row = [formattedData.id, "{}".format(formattedData.reactionIDs), "{}".format(formattedData.VmaxData.medianEntry.vmax),
		"{}".format(formattedData.VmaxData.medianEntry.__dict__), "{}".format(formattedData.VmaxData.minEntry.__dict__), "{}".format(formattedData.VmaxData.maxEntry.__dict__),
		"{}".format(formattedData.VmaxData.closestEntries[0].proximity)]
		rows.append(row)


	rows = topRow + rows

	wb = Workbook()
	dest_filename = 'Blueberry2.xlsx'
	ws1 = wb.active
	ws1.title = "Kinetics"


	row = 1
	while row < len(rows)+1:
		col = 1
		while col < len(rows[0])+1:
			_ = ws1.cell(column=col, row=row, value=rows[row-1][col-1])
			col = col+1
		row = row+1

	wb.save(filename = dest_filename)

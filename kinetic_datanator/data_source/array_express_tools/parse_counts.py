

class Measurement():
	def __init__(self, gene_name, count_value):
		self.gene_name = gene_name
		self.value = float(count_value)
		self.value_relative_to_total_counts = None

class Sample():

	def __init__(self, sample_name):
		self.sample_name = sample_name
		self.measurements = []
		self.whole_genome = False
		total_value_all_counts = None


def get_data_from_file(text_file):
	file = open('counts.txt')

	samples = []

	lines = file.readlines()
	content = [x.strip().split(' ') for x in lines]
	for sample_name in content[0]:#.split(' '):
		samples.append(Sample(sample_name=sample_name))
	for num, sample in enumerate(samples):
		measurements = []
		for row in content[1:]:
			measurements.append(Measurement(row[0], row[num+1]))
		sample.measurements = measurements

	for sample in samples:
		total = 0
		for m in sample.measurements:
			total = total+m.value
		sample.total_value_all_counts = total

		for measurement in sample.measurements:
			measurement.value_relative_to_total_counts = measurement.value/sample.total_value_all_counts
			#print measurement.value_relative_to_total_counts


	return samples



if __name__ == '__main__':  
	file = open('counts.txt')
	samples = get_data_from_file(file)

	for sample in samples:
		print sample.total_value_all_counts
		print len(sample.measurements)
		"""
		total = 0
		for m in sample.measurements:
			total = total+m.value
		print total
		"""


	

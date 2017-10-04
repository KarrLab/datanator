import unittest
from kinetic_datanator.data_source.array_express_tools import parse_counts



class QuickTest(unittest.TestCase):

	def test_get_data_from_file(self):
		file = open('tests/data_source/array_express_tools/counts.txt')
		samples = parse_counts.get_data_from_file(file)

		values = []
		for sample in samples:
			if sample.sample_name == 'aHSC':
				for m in sample.measurements:
					if m.gene_name == 'ENSMUSG00000000154':
						values.append(m.value)
		self.assertEqual(values, [2780, 2025, 4141])

		values = []
		for sample in samples:
			if sample.sample_name == 'MPP1':
				for m in sample.measurements:
					if m.gene_name == 'ENSMUSG00000061848':
						values.append(m.value)
		self.assertEqual(values, [114, 64, 86])

		total_counts = [128332069.0,
						27900638.0,
						103846727.0,
						123455535.0,
						85306936.0,
						125538695.0,
						132944642.0,
						105929159.0,
						141834239.0,
						]
		test_counts = []
		for sample in samples:
			test_counts.append(sample.total_value_all_counts)
		self.assertEqual(test_counts, total_counts)


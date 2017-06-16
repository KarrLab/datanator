from kinetic_datanator.kinetic_datanator.data_source import array_express
import unittest
import array_express
import requests.exceptions
import tempfile
import os
import shutil

"""
a = array_express.ArrayExpress()
db_session = a.session

accesssion_nums = [
	"E-MTAB-5207",
	"E-MTAB-5219",
	"E-MTAB-3145",
	"E-MTAB-5661",
	"E-MTAB-5385",
	"E-MTAB-5207",
	]

a.load_experiments(experiment_ids=accesssion_nums)
experiments = db_session.query(array_express.Experiment.id).all()
db_session.commit()


for m in db_session.query(array_express.Experiment).all():
	print [getattr(m, x.__str__().split('.')[1]) for x in array_express.Experiment.__table__.columns]
"""

class TestArrayExpress(unittest.TestCase):

	@classmethod
	def setUpClass(cls):
		cls.cache_dirname = tempfile.mkdtemp()
		src = array_express.ArrayExpress(cache_dirname=cls.cache_dirname, download_backup=False, load_content=False)
		src.session.close()
		src.engine.dispose()

	@classmethod
	def tearDownClass(cls):
		shutil.rmtree(cls.cache_dirname)


	def setUp(self):
		self.src = array_express.ArrayExpress(cache_dirname=self.cache_dirname, download_backup=False, load_content=False)

	def tearDown(self):
		src = self.src
		src.session.close()
		src.engine.dispose()

	def test_load_experiments(self):
		src = self.src
		session = src.session
		accesssion_nums = [
			"E-MTAB-5219",
			"E-MTAB-3145",
			"E-MTAB-5661",
			"E-MTAB-5385",
			"E-MTAB-5207",
		]
		src.load_experiments(experiment_ids=accesssion_nums)
		q = session.query(array_express.Experiment)

		"""
		print session.query(array_express.Experiment.id).all()
		
		for thing in q:
			print thing.id
		for m in session.query(array_express.Experiment).all():
			print [getattr(m, x.__str__().split('.')[1]) for x in array_express.Experiment.__table__.columns]

		#print session.query(array_express.Experiment.id).all()
		for thing in session.query(array_express.Experiment).\
		filter(array_express.Experiment.id=="E-MTAB-5207"):
			print thing.experiment_type
		"""
		
			
		self.assertEqual(q.count(), 5)
		#print [c.description for c in q.all()][0:5]
		self.assertEqual([c.id for c in q.all()], [
			'E-MTAB-5207', 
			'E-MTAB-3145', 
			'E-MTAB-5385', 
			'E-MTAB-5219', 
			'E-MTAB-5661'
		])

		self.assertEqual([c.name for c in q.all()], [
			'Dynamic regulation of VEGF-inducible genes by an ERK-ERG-p300 transcriptional network', 
			'Microarray analysis of palatal shelves from wild-type versus p63-null mouse embryos', 
			'Transcription profiling by array of neonatal and children dermal fibroblasts isolated from cleft lip and adult dermal fibroblasts of breast and face', 
			'Global expression profiling time series of early human embryonic stem cell differentiation towards the mesoderm lineage', 
			'Deciphering the relationship between polycomb repression and stochastic gene expression from single-cell RNA-seq data'
		])

		self.assertEqual([c.experiment_type for c in q.all()], [
			'transcription profiling by array',
			'transcription profiling by array', 
			'transcription profiling by array', 
			'transcription profiling by array', 
			'RNA-seq of coding RNA from single cells'
		])

		self.assertEqual([c.organism for c in q.all()], [
			'Homo sapiens', 
			'Mus musculus', 
			'Homo sapiens', 
			'Homo sapiens', 
			'Mus musculus'
		])

		self.assertEqual([c.description for c in q.all()], [
			'Effect of ERG knockdown, with and without VEGF stimulation, on vascular endothelial cell (HUVEC) gene expression', 
			'Mutations in the transcription factor p63 underlie of a series of human malformation syndromes which are defined by a combination of epidermal, limb and craniofacial abnormalities including cleft lip and palate. Transcription profiling was performed to determine the role of p63 in vivo mouse palatal shelves. Microarray analysis was done of palatal shelves dissected from E14.0 wild-type versus p63-null mouse embryos.', 
			'We measured expression profiles of skin fibroblasts harvested during lip cleft surgery of neonates (until 16 days after birth) and children (until 10 years after birth) and compared them to fibroblast isolated from adult skin of breast, scar, and face. We observed clear differences between the cleft and facial fibroblasts and the fibroblasts of breast, as the former originate from the neural crest. Further, we observed differences between adult face fibroblasts and children fibroblasts.', 
			'Human pluripotent stem cells (hPSCs) is a promising cell type with capacities of self-renewal and pluripotency, and they can differentiate into in principal any cell types in the body. To fully make use of this unique cell type and develop efficient and reproducible differentiation protocols, more knowledge is needed about regulatory mechanisms and molecular pathways important to direct the differentiation towards specific functional cell types. The aim of this experiment was to study the events and regulators governing early differentiation human embryonic stem cells (hESCs) towards the mesoderm lineage. To this end, global expression profiling was performed for 11 time points on the hESC line SA121 (Takara-Clontech). Three biological replicates (A, B and C) were taken at day 0 (before differentiation) and on each day thereafter for ten days.', 
			'Polycomb repressive complexes are important histone modifiers, which silence gene expression, yet there exists a subset of polycomb-bound genes actively transcribed by RNA polymerase II. To investigate the switching between polycomb-repressed and active states, we sequence mRNA from OS25 mouse embryonic stem cells cultured in serum/LIF. To validate our finding that polycomb modulates stochastic gene expression and transcriptional bursting, we perform knockout experiments and we sequence mRNA from Ring1A knockout (untreated) and Ring1A/B double knockout cells with constitutive Ring1A knockout and tamoxifen-inducible conditional Ring1B knockout.'])


		q = session.query(array_express.ExperimentDesign)

		self.assertEqual([c.experiment_id for c in q.all()], [2, 3, 3, 3, 3, 3, 4])

		self.assertEqual([c.name for c in q.all()], [
			'genetic modification design', 
			'case control design', 
			'cell type comparison design', 
			'clinical history design', 
			'development or differentiation design', 
			'disease state design', 
			'development or differentiation design'
			])


	def test_load_samples(self):

		src = self.src
		session = src.session
		accesssion_nums = [
			"E-MTAB-5219",
			"E-MTAB-3145",
			#"E-MTAB-5661",
			#"E-MTAB-5385",
			#"E-MTAB-5207",
		]
		src.load_experiments(experiment_ids=accesssion_nums)
		exp = session.query(array_express.Experiment)
		src.load_samples(exp)

		q = session.query(array_express.Sample)
		#session.commit()
		#print [c.experiment_id for c in q.all()][0:10]

		self.assertEqual([c.index for c in q.all()[0:10]], [
			'MUT_41', 
			'MUT_42', 
			'MUT_43', 
			'WT_40', 
			'WT_45', 
			'WT_49', 
			'01_d0_a', 
			'02_d0_b', 
			'03_d0_c', 
			'04_d1_a'
		])

		self.assertEqual([c.index for c in q.all()[0:10]], [c.name for c in q.all()[0:10]])
		self.assertEqual([c.experiment_id for c in q.all()[0:10]], [1, 1, 1, 1, 1, 1, 2, 2, 2, 2])


		q = session.query(array_express.Characteristic)

		self.assertEqual([c.name for c in q.all()[0:10]], [
			'organism', 
			'genotype', 
			'organism part', 
			'developmental stage', 
			'organism', 
			'genotype', 
			'organism part', 
			'developmental stage', 
			'organism', 
			'genotype'
			])


		self.assertEqual([c.value for c in q.all()[0:10]], [
			'Mus musculus', 
			'p63-/-', 
			'palate', 
			'embryonic day 14.5', 
			'Mus musculus', 
			'p63-/-', 
			'palate', 
			'embryonic day 14.5', 
			'Mus musculus', 
			'p63-/-'])


		q = session.query(array_express.Variable)
		#print [c.unit for c in q.all()][0:10]

		self.assertEqual([c.name for c in q.all()[0:10]], [
			'genotype', 
			'genotype', 
			'genotype', 
			'genotype', 
			'genotype', 
			'genotype', 
			'time', 
			'time', 
			'time', 
			'time'
		])

		self.assertEqual([c.value for c in q.all()[0:10]], [
			'p63-/-', 
			'p63-/-', 
			'p63-/-', 
			'wild type genotype', 
			'wild type genotype', 
			'wild type genotype', 
			'0', 
			'0', 
			'0', 
			'1'
		])


		self.assertEqual([c.unit for c in q.all()[0:10]], [None, None, None, None, None, None, 'day', 'day', 'day', 'day'])


		#for m in session.query(array_express.Variable).all()[0:10]:
		#	print [getattr(m, x.__str__().split('.')[1]) for x in array_express.Variable.__table__.columns]

		#q = session.query(array_express.Sample)
		#print [c.name for c in q.all()][0:5]

""" Tests of the command line utility

:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-04-12
:Copyright: 2017, Karr Lab
:License: MIT
"""

from capturer import CaptureOutput
from kinetic_datanator.__main__ import App
import unittest


class TestCli(unittest.TestCase):

    @unittest.skip('implement me')
    def test_GetKineticsController(self):
        pass

    @unittest.skip('implement me')
    def test_GenerateTemplateController(self):
        pass

    def test_taxonomy_get_rank(self):
        with App(argv=['taxonomy', 'get-rank', '2097']) as app:
            with CaptureOutput() as capturer:
                app.run()
                self.assertEqual(capturer.get_text(), "species")

        with App(argv=['taxonomy', 'get-rank', 'Mycoplasma genitalium']) as app:
            with CaptureOutput() as capturer:
                app.run()
                self.assertEqual(capturer.get_text(), "species")

        with App(argv=['taxonomy', 'get-rank', 'Mycoplasma genitalium XXX']) as app:
            self.assertRaises(ValueError, lambda: app.run())

    def test_taxonomy_get_parents(self):
        with App(argv=['taxonomy', 'get-parents', 'bacteria']) as app:
            with CaptureOutput() as capturer:
                app.run()
                self.assertEqual(capturer.get_text(), "root\ncellular organisms")

        with App(argv=['taxonomy', 'get-parents', 'XXX']) as app:
            self.assertRaises(ValueError, lambda: app.run())

    def test_taxonomy_get_common_ancestor(self):
        with App(argv=['taxonomy', 'get-common-ancestor', 'Mycoplasma genitalium', 'Mycoplasma pneumoniae']) as app:
            with CaptureOutput() as capturer:
                app.run()
                self.assertEqual(capturer.get_text(), 'Mycoplasma')

        with App(argv=['taxonomy', 'get-common-ancestor', 'Mycoplasma genitalium', 'XXX']) as app:
            with CaptureOutput() as capturer:
                self.assertRaises(ValueError, lambda: app.run())

    def test_taxonomy_get_distance_to_common_ancestor(self):
        with App(argv=['taxonomy', 'get-distance-to-common-ancestor', 'Mycoplasma genitalium', 'Mycoplasma pneumoniae']) as app:
            with CaptureOutput() as capturer:
                app.run()
                self.assertEqual(float(capturer.get_text()), 1.)

        with App(argv=['taxonomy', 'get-distance-to-common-ancestor', 'Mycoplasma genitalium', 'XXX']) as app:
            self.assertRaises(ValueError, lambda: app.run())

    def test_taxonomy_get_distance_to_root(self):
        with App(argv=['taxonomy', 'get-distance-to-root', 'bacteria']) as app:
            with CaptureOutput() as capturer:
                app.run()
                self.assertEqual(float(capturer.get_text()), 2.)

        with App(argv=['taxonomy', 'get-distance-to-root', 'XXX']) as app:
            self.assertRaises(ValueError, lambda: app.run())

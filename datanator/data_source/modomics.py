""" Use BpForms to gather rRNA and tRNA modification information from 
`MODOMICS <https://iimcb.genesilico.pl/modomics/>`.

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2020-04-23
:Copyright: 2019, Karr Lab
:License: MIT
"""

import bpforms
import bs4
import csv
import os
import requests
import requests_cache


class Modomics(object):
    """ Use BpForms to gather rRNA and tRNA modification information from 
    `MODOMICS <https://iimcb.genesilico.pl/modomics/>`.
    """
    ENDPOINT = 'https://iimcb.genesilico.pl/modomics/sequences/'

    def run(self):
        # create dict of MODOMICS single character monomer codes
        modomics_short_code_to_monomer = {}
        for monomer in bpforms.rna_alphabet.monomers.values():
            for identifier in monomer.identifiers:
                if identifier.ns == 'modomics.short_name':
                    modomics_short_code_to_monomer[identifier.id] = monomer
        modomics_short_code_to_monomer['a'] = bpforms.rna_alphabet.monomers.get('A')
        modomics_short_code_to_monomer['c'] = bpforms.rna_alphabet.monomers.get('C')
        modomics_short_code_to_monomer['g'] = bpforms.rna_alphabet.monomers.get('G')
        modomics_short_code_to_monomer['u'] = bpforms.rna_alphabet.monomers.get('U')
        modomics_short_code_to_monomer['b'] = bpforms.rna_alphabet.monomers.get('0522U')
        modomics_short_code_to_monomer['B'] = bpforms.rna_alphabet.monomers.get('0C')
        modomics_short_code_to_monomer['E'] = bpforms.rna_alphabet.monomers.get('662A')
        modomics_short_code_to_monomer['h'] = bpforms.rna_alphabet.monomers.get('21511U')
        # modomics_short_code_to_monomer['H'] = bpforms.rna_alphabet.monomers.get('0C')
        modomics_short_code_to_monomer['J'] = bpforms.rna_alphabet.monomers.get('0U')
        modomics_short_code_to_monomer['l'] = bpforms.rna_alphabet.monomers.get('253U')
        modomics_short_code_to_monomer['L'] = bpforms.rna_alphabet.monomers.get('2G')
        modomics_short_code_to_monomer['K'] = bpforms.rna_alphabet.monomers.get('1G')
        modomics_short_code_to_monomer['M'] = bpforms.rna_alphabet.monomers.get('42C')
        # modomics_short_code_to_monomer['N'] = bpforms.rna_alphabet.monomers.get('?U')
        modomics_short_code_to_monomer['P'] = bpforms.rna_alphabet.monomers.get('9U')
        modomics_short_code_to_monomer['R'] = bpforms.rna_alphabet.monomers.get('22G')
        modomics_short_code_to_monomer['T'] = bpforms.rna_alphabet.monomers.get('5U')
        modomics_short_code_to_monomer['Z'] = bpforms.rna_alphabet.monomers.get('09U')
        modomics_short_code_to_monomer['7'] = bpforms.rna_alphabet.monomers.get('7G')
        modomics_short_code_to_monomer['#'] = bpforms.rna_alphabet.monomers.get('0G')
        modomics_short_code_to_monomer[':'] = bpforms.rna_alphabet.monomers.get('0A')
        modomics_short_code_to_monomer['='] = bpforms.rna_alphabet.monomers.get('6A')
        modomics_short_code_to_monomer['?'] = bpforms.rna_alphabet.monomers.get('5C')
        modomics_short_code_to_monomer['λ'] = bpforms.rna_alphabet.monomers.get('04C')
        modomics_short_code_to_monomer['"'] = bpforms.rna_alphabet.monomers.get('1A')
        modomics_short_code_to_monomer["'"] = bpforms.rna_alphabet.monomers.get('3C')
        modomics_short_code_to_monomer[','] = bpforms.rna_alphabet.monomers.get('522U')
        modomics_short_code_to_monomer['\\'] = bpforms.rna_alphabet.monomers.get('05U')
        modomics_short_code_to_monomer['ℑ'] = bpforms.rna_alphabet.monomers.get('00G')
        modomics_short_code_to_monomer[']'] = bpforms.rna_alphabet.monomers.get('19U')
        modomics_short_code_to_monomer['ˆ'] = bpforms.rna_alphabet.monomers.get('00A')
        modomics_short_code_to_monomer['gluQtRNA'] = bpforms.rna_alphabet.monomers.get('105G')
        modomics_short_code_to_monomer['m22G'] = bpforms.rna_alphabet.monomers.get('22G')
        # modomics_short_code_to_monomer[';'] = bpforms.rna_alphabet.monomers.get('?G')
        # modomics_short_code_to_monomer['<'] = bpforms.rna_alphabet.monomers.get('?C')

        # output directory
        out_dir = os.path.dirname(__file__)

        # create cache for web queries
        cache_name = os.path.join(out_dir, 'modomics')
        session = requests_cache.core.CachedSession(cache_name, backend='sqlite', expire_after=None)
        session.mount(self.ENDPOINT, requests.adapters.HTTPAdapter(max_retries=5))

        # parse rRNA and tRNA data
        monomer_codes = {}
        rrna_forms = self._run_rrna(session, modomics_short_code_to_monomer, monomer_codes,
                                    os.path.join(out_dir, 'modomics.rrna.tsv'))
        trna_forms = self._run_trna(
            session, modomics_short_code_to_monomer, monomer_codes, os.path.join(out_dir, 'modomics.trna.tsv'))

        # return results
        return rrna_forms, trna_forms

    def _run_rrna(self, session, modomics_short_code_to_monomer, monomer_codes, out_filename):
        response = session.get(self.ENDPOINT, params={
            'RNA_type': 'rRNA',
            'RNA_subtype': 'all',
            'organism': 'all species',
            'vis_type': 'Modomics symbols',
        })

        response.raise_for_status()

        doc = bs4.BeautifulSoup(response.text, 'lxml')
        table = doc.find('table', {'id': 'tseq'})
        tbody = table.find('tbody')
        rows = tbody.find_all('tr')
        rna_forms = []
        for row in rows:
            if not isinstance(row, bs4.element.Tag):
                continue

            cells = row.find_all('td')

            rna_form = bpforms.RnaForm()
            unsupported_codes = set()
            for child in cells[5].children:
                if child.name is None or child.name == 'span':
                    if child.name is None:
                        text = str(child)
                    else:
                        text = child.text

                    for code in text.strip().replace('-', '').replace('_', ''):
                        monomer = modomics_short_code_to_monomer.get(code, None)
                        if monomer is None:
                            unsupported_codes.add(code)
                            monomer = bpforms.Monomer(id=code)
                        else:
                            monomer_codes[code] = monomer
                        rna_form.seq.append(monomer)
                elif child.name == 'a':
                    code = child.get('href').replace('/modomics/modifications/', '')
                    monomer = modomics_short_code_to_monomer.get(code, None)
                    if monomer is None:
                        unsupported_codes.add(code)
                        monomer = bpforms.Monomer(id=code)
                    else:
                        monomer_codes[code] = monomer
                    rna_form.seq.append(monomer)
                else:
                    raise Exception('Unsupported child {}'.format(child.name))

            rna_forms.append({
                'GenBank': cells[0].find('a').text.strip(),
                'Organism': cells[3].text.strip(),
                'Organellum': cells[4].text.strip(),
                'Type': cells[2].text.strip(),
                'Sequence (MODOMICS)': cells[5].text.strip().replace('-', '').replace('_', ''),
            })
            self._analyze_form(rna_form, unsupported_codes, rna_forms[-1])

        # save results to tab-separated file
        self._save_results(rna_forms, ['GenBank', 'Type'], out_filename)

        return rna_forms

    def _run_trna(self, session, modomics_short_code_to_monomer, monomer_codes, out_filename):
        response = session.get(self.ENDPOINT, params={
            'RNA_type': 'tRNA',
            'RNA_subtype': 'all',
            'organism': 'all species',
            'vis_type': 'Modomics symbols',
        })
        response.raise_for_status()

        doc = bs4.BeautifulSoup(response.text, 'lxml')
        table = doc.find('table', {'id': 'tseq'})
        tbody = table.find('tbody')
        rows = tbody.find_all('tr')
        rna_forms = []

        for row in rows:
            cells = row.find_all('td')

            rna_form = bpforms.RnaForm()
            unsupported_codes = set()
            for child in cells[5].children:
                if child.name is None or child.name == 'span':
                    if child.name is None:
                        text = str(child)
                    else:
                        text = child.text

                    for code in text.strip().replace('-', '').replace('_', ''):
                        monomer = modomics_short_code_to_monomer.get(code, None)
                        if monomer is None:
                            unsupported_codes.add(code)
                            monomer = bpforms.Monomer(id=code)
                        else:
                            monomer_codes[code] = monomer
                        rna_form.seq.append(monomer)
                elif child.name == 'a':
                    code = child.get('href').replace('/modomics/modifications/', '')
                    monomer = modomics_short_code_to_monomer.get(code, None)
                    if monomer is None:
                        unsupported_codes.add(code)
                        monomer = bpforms.Monomer(id=code)
                    else:
                        monomer_codes[code] = monomer
                    rna_form.seq.append(monomer)
                else:
                    raise Exception('Unsupported child {}'.format(child.name))

            rna_forms.append({
                'Amino acid type': cells[1].text.strip(),
                'Anticodon': cells[2].text.strip(),
                'Organism': cells[3].text.strip(),
                'Organellum': cells[4].text.strip(),
                'Sequence (MODOMICS)': cells[5].text.strip().replace('-', '').replace('_', ''),
            })
            self._analyze_form(rna_form, unsupported_codes, rna_forms[-1])

        # save results to tab-separated file
        self._save_results(rna_forms, ['Amino acid type', 'Anticodon'], out_filename)

        return rna_forms

    def _analyze_form(self, rna_form, unsupported_codes, results_dict):
        results_dict['Sequence (BpForms)'] = str(rna_form)
        results_dict['Sequence (IUPAC)'] = canonical_seq = rna_form.get_canonical_seq()
        results_dict['Length'] = len(rna_form.seq)

        results_dict['Number of modifications'] = len(rna_form.seq) \
            - rna_form.seq.count(bpforms.rna_alphabet.monomers.A) \
            - rna_form.seq.count(bpforms.rna_alphabet.monomers.C) \
            - rna_form.seq.count(bpforms.rna_alphabet.monomers.G) \
            - rna_form.seq.count(bpforms.rna_alphabet.monomers.U)
        results_dict['Number of modified A'] = canonical_seq.count('A') - rna_form.seq.count(bpforms.rna_alphabet.monomers.A)
        results_dict['Number of modified C'] = canonical_seq.count('C') - rna_form.seq.count(bpforms.rna_alphabet.monomers.C)
        results_dict['Number of modified G'] = canonical_seq.count('G') - rna_form.seq.count(bpforms.rna_alphabet.monomers.G)
        results_dict['Number of modified U'] = canonical_seq.count('U') - rna_form.seq.count(bpforms.rna_alphabet.monomers.U)

        if unsupported_codes:
            results_dict['BpForms errors'] = 'MODOMICS sequence uses monomeric forms {}'.format(
                ', '.join(unsupported_codes))
        else:
            results_dict['Formula'] = str(rna_form.get_formula())
            results_dict['Molecular weight'] = rna_form.get_mol_wt()
            results_dict['Charge'] = rna_form.get_charge()

            canonical_form = bpforms.RnaForm().from_str(canonical_seq)
            results_dict['Canonical formula'] = str(canonical_form.get_formula())
            results_dict['Canonical molecular weight'] = canonical_form.get_mol_wt()
            results_dict['Canonical charge'] = canonical_form.get_charge()

            results_dict['Extra formula'] = str(rna_form.get_formula() - canonical_form.get_formula())
            results_dict['Extra molecular weight'] = rna_form.get_mol_wt() - canonical_form.get_mol_wt()
            results_dict['Extra charge'] = rna_form.get_charge() - canonical_form.get_charge()

            results_dict['BpForms errors'] = ' '.join(rna_form.validate())

    def _save_results(self, rna_forms, variable_field_names, out_filename):
        with open(out_filename, 'w') as file:
            writer = csv.DictWriter(file,
                                    fieldnames=variable_field_names + [
                                        'Organism', 'Organellum',
                                        'Sequence (MODOMICS)',
                                        'Sequence (BpForms)',
                                        'Sequence (IUPAC)',
                                        'Length',
                                        'Number of modifications',
                                        'Number of modified A',
                                        'Number of modified C',
                                        'Number of modified G',
                                        'Number of modified U',
                                        'Formula', 'Molecular weight', 'Charge',
                                        'Canonical formula', 'Canonical molecular weight', 'Canonical charge',
                                        'Extra formula', 'Extra molecular weight', 'Extra charge',
                                        'BpForms errors'],
                                    dialect='excel-tab')
            writer.writeheader()

            for rna_form in rna_forms:
                writer.writerow(rna_form)

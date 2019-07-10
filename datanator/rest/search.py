import functools
from datanator.rest.query import front_end_query
from flask import jsonify
from flask import (
    Blueprint, flash, g, redirect, render_template, request, session, url_for
)
import json
import re


bp = Blueprint('/search', __name__)
bp_r = Blueprint('/results', __name__)

#have not used
@bp.route('/search', methods=('GET', 'POST'))
def search():
    if request.method == 'POST':
        molecule_name = request.form['molecule_name']
        organism_name = request.form['organism_name']
        abstract_default = request.form['abstract_default']
        session['molecule_name'] = molecule_name
        session['organism_name'] = organism_name
        session['abstract_default'] = abstract_default
        return redirect(url_for('/search.results'))


    return render_template('/search.html')
@bp_r.route('/results/<path:molecule_name>/<string:organism_name>', methods=('GET', 'POST'))
@bp_r.route('/results/<path:molecule_name>/<string:organism_name>/<string:abstract_default>', methods=('GET', 'POST'))
def results(molecule_name, organism_name, abstract_default=False):
    #q = front_end_query.QueryFrontEnd()
    print(molecule_name)
  
    print("{} {}".format(molecule_name, organism_name))
    if abstract_default:
        abstract_default = (abstract_default.lower() == "true")

    list_jsons = front_end_query.QueryFrontEnd().molecule_name_query(molecule_name, organism_name, abstract_default=abstract_default)

    return(jsonify(list_jsons))
    #return render_template('/results.html', results=[the_json_1,the_json_2])

    #return render_template('/results.html', results=json.dumps(the_json))
if __name__ == '__main__':
    bp_r.run(debug=True, host='0.0.0.0')
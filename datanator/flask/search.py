import functools
from . import front_end_query
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
        session['molecule_name'] = molecule_name
        session['organism_name'] = organism_name
        return redirect(url_for('/search.results'))


    return render_template('/search.html')



@bp.route('/results/<path:molecule_name>/<string:organism_name>', methods=('GET', 'POST'))
def results(molecule_name, organism_name):
    #q = front_end_query.QueryFrontEnd()
    print(molecule_name)
  
    print("{} {}".format(molecule_name, organism_name))

    list_jsons = front_end_query.molecule_name_query(molecule_name, organism_name)

    return(jsonify(list_jsons))
    #return render_template('/results.html', results=[the_json_1,the_json_2])

    #return render_template('/results.html', results=json.dumps(the_json))

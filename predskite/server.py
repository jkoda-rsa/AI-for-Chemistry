from flask import Flask, render_template, request, jsonify
from molview import display_3D_structure
from pcepredicted import estimate_etl_pce, estimate_htl_pce
from getsmiles import get_smiles

import os
import shutil
import sys

app = Flask(__name__)
#app.config.update(SERVER_NAME='127.0.0.1:5000')


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/results/', methods=['POST'])
def display_request():
    smiles = request.form["input"]
    action = request.form["action"]
    pdb_block = display_3D_structure(smiles)

    if action == "ETL":
        result = estimate_etl_pce(smiles)
    elif action == "HTL":
        result = estimate_htl_pce(smiles)
    else:
        result = None  # handle unexpected action value

    smiles_str = get_smiles(smiles)  # Generate the SMILES string

    return render_template('results.html', res={"result": result, "pdb_block": pdb_block, "input": smiles, "smiles_str": smiles_str})



if __name__ == '__main__':  
    app.run(debug=True)







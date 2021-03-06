from flask import Flask, render_template, request
from utils import utils

app = Flask(__name__)



@app.route("/get_pca", methods=["POST"])
def get_pca():
    smiles_list = request.json.get("smiles_list")
    pca = utils.get_pca(smiles_list)
    return {"pca": pca}


@app.route("/get_similarities", methods=["POST"])
def get_similarities():
    try:
        r = request.get_json()
        similarities = utils.get_similarities(smiles0=r['smiles0'], smiles_list=r['smiles_list'])
        return {"msg": "success", "similarities": similarities}
    except:
        return {"msg": "fail", "similarities": []}


@app.route("/draw_smiles", methods=["POST"])
def draw_smiles():
    """Draw a smiles string as an SVG string"""
    try:
        smiles = request.json.get("smiles")
        drawing = utils.smiles_to_svg(smiles)
        return {"msg": "success", "svg": drawing, "smiles": smiles}
    except Exception as e:
        print(e)
        return {"msg": "fail"}

@app.route("/")
def index():
    return render_template("index.html")

@app.route("/about")
def about():
    return render_template("about.html")

@app.errorhandler(404)
def not_found_error(error):
    """This page is returned when application experiences a 404 error"""
    return render_template("error.html", error_num=404), 404

@app.errorhandler(500)
def internal_error(error):
    """This page is returned when application experiences a 500 error"""
    return render_template("error.html", error_num=500), 500

if __name__ == "__main__":
    app.run(debug=True)

from flask import Flask, render_template, request
from utils import utils

app = Flask(__name__)


@app.route("/")
def index():
    return render_template("index.html")

@app.route("/get_similarities", methods=["POST"])
def get_similarities():
    r = request.get_json()
    similarities = utils.get_similarities(
        smiles0=r['smiles0'],
        smiles_list=r['smiles_list'],
    )
    return {"similarities": similarities}



if __name__ == "__main__":
    app.run(debug=True)

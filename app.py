
from flask import Flask, render_template, request
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import base64
import requests
import os

app = Flask(__name__)

def get_smiles_from_name(name):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/IsomericSMILES/JSON"
        res = requests.get(url)
        data = res.json()
        smiles = data['PropertyTable']['Properties'][0]['IsomericSMILES']
        return smiles
    except:
        return None

@app.route("/", methods=["GET", "POST"])
def home():
    output = ""
    img = ""
    molblock = ""
    name = ""

    if request.method == "POST":
        name = request.form.get("name")
        smiles = request.form.get("smiles")

        # 🔥 AUTO FETCH FROM PUBCHEM
        if not smiles and name:
            smiles = get_smiles_from_name(name)

        if not smiles:
            output = "❌ Could not fetch SMILES"
            return render_template("index.html", output=output, img=img, molblock=molblock, name=name)

        mol = Chem.MolFromSmiles(smiles)

        if mol is None:
            output = "❌ Invalid molecule"
            return render_template("index.html", output=output, img=img, molblock=molblock, name=name)

        mol = Chem.AddHs(mol)
        Chem.AssignStereochemistry(mol, force=True)

        # 🔬 CHIRALITY
        output += "🔬 Chiral Centers:\n"
        centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)

        if not centers:
            output += "None\n"
        else:
            for i, c in centers:
                output += f"Atom {i} → {c}\n"

        # 🧪 GROUPS
        benzyl = Chem.MolFromSmarts("[CH2]c1ccccc1")
        phenyl = Chem.MolFromSmarts("c1ccccc1")
        benzoyl = Chem.MolFromSmarts("C(=O)c1ccccc1")

        output += f"\nBenzyl: {len(mol.GetSubstructMatches(benzyl))}"
        output += f"\nPhenyl: {len(mol.GetSubstructMatches(phenyl))}"
        output += f"\nBenzoyl: {len(mol.GetSubstructMatches(benzoyl))}"

        # 🖼️ 2D IMAGE
        img_data = Draw.MolToImage(mol)
        import io
        buf = io.BytesIO()
        img_data.save(buf, format="PNG")
        img = base64.b64encode(buf.getvalue()).decode()

        # 🧊 3D STRUCTURE
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol)
        molblock = Chem.MolToMolBlock(mol)

    return render_template("index.html", output=output, img=img, molblock=molblock, name=name)

if __name__ == "__main__":
    port = int(os.environ.get("PORT", 10000))
    app.run(host="0.0.0.0", port=port)

from flask import Flask, render_template, request
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import base64
from io import BytesIO

app = Flask(__name__)

def analyze(smiles):
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return "Invalid SMILES!", None, None

    mol = Chem.AddHs(mol)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)

    result = ""

    # =========================
    # 🔬 CHIRALITY
    # =========================
    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    result += "Chiral Centers:\n"

    if len(centers) == 0:
        result += "None\n"
    else:
        for idx, config in centers:
            result += f"Atom {idx} → {config}\n"

    # =========================
    # 🧪 BENZYL GROUP
    # =========================
    benzyl_pattern = Chem.MolFromSmarts("[CH2]c1ccccc1")
    benzyl_matches = mol.GetSubstructMatches(benzyl_pattern)

    result += f"\nBenzyl groups: {len(benzyl_matches)}\n"

    # =========================
    # 🧪 PHENYL RINGS
    # =========================
    phenyl_pattern = Chem.MolFromSmarts("c1ccccc1")
    phenyl_matches = mol.GetSubstructMatches(phenyl_pattern)

    result += f"Phenyl rings: {len(phenyl_matches)}\n"

    # =========================
    # 🧪 BENZOYL GROUP
    # =========================
    benzoyl_pattern = Chem.MolFromSmarts("C(=O)c1ccccc1")
    benzoyl_matches = mol.GetSubstructMatches(benzoyl_pattern)

    result += f"Benzoyl groups: {len(benzoyl_matches)}\n"

    # =========================
    # 🧪 2D IMAGE
    # =========================
    img = Draw.MolToImage(mol)
    buffer = BytesIO()
    img.save(buffer, format="PNG")
    img_str = base64.b64encode(buffer.getvalue()).decode()

    # =========================
    # 🧊 3D STRUCTURE (CLEAN)
    # =========================
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)

    mol = Chem.RemoveHs(mol)  # clean look
    mol_block = Chem.MolToMolBlock(mol)

    return result, img_str, mol_block


@app.route("/", methods=["GET", "POST"])
def home():
    output = ""
    img = None
    molblock = None
    name = ""

    if request.method == "POST":
        name = request.form.get("name", "Unknown Compound")
        smiles = request.form.get("smiles", "")
        output, img, molblock = analyze(smiles)

    return render_template(
        "index.html",
        output=output,
        img=img,
        molblock=molblock,
        name=name
    )


if __name__ == "__main__":
    app.run(debug=True)
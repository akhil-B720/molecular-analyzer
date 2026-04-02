from flask import Flask, render_template, request
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, Descriptors
import base64
import requests
import io
import os

app = Flask(__name__)

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

def get_compound_data(name):
    """Fetch rich compound data from PubChem by name."""
    result = {
        "smiles": None,
        "iupac_name": None,
        "molecular_formula": None,
        "molecular_weight": None,
        "synonyms": [],
        "cid": None,
        "inchikey": None,
        "xlogp": None,
        "hbd": None,
        "hba": None,
        "tpsa": None,
        "rotatable_bonds": None,
        "charge": None,
        "complexity": None,
        "canonical_smiles": None,
    }

    try:
        # Step 1: Resolve name → CID
        cid_url = f"{PUBCHEM_BASE}/compound/name/{requests.utils.quote(name)}/cids/JSON"
        r = requests.get(cid_url, timeout=8)
        r.raise_for_status()
        cid = r.json()["IdentifierList"]["CID"][0]
        result["cid"] = cid
    except Exception:
        return result

    try:
        # Step 2: Fetch all properties in one call
        props = [
            "IsomericSMILES", "CanonicalSMILES", "IUPACName",
            "MolecularFormula", "MolecularWeight", "InChIKey",
            "XLogP", "HBondDonorCount", "HBondAcceptorCount",
            "TPSA", "RotatableBondCount", "Charge", "Complexity"
        ]
        prop_url = f"{PUBCHEM_BASE}/compound/cid/{cid}/property/{','.join(props)}/JSON"
        r = requests.get(prop_url, timeout=8)
        r.raise_for_status()
        p = r.json()["PropertyTable"]["Properties"][0]

        result["smiles"]            = p.get("IsomericSMILES")
        result["canonical_smiles"]  = p.get("CanonicalSMILES")
        result["iupac_name"]        = p.get("IUPACName")
        result["molecular_formula"] = p.get("MolecularFormula")
        result["molecular_weight"]  = p.get("MolecularWeight")
        result["inchikey"]          = p.get("InChIKey")
        result["xlogp"]             = p.get("XLogP")
        result["hbd"]               = p.get("HBondDonorCount")
        result["hba"]               = p.get("HBondAcceptorCount")
        result["tpsa"]              = p.get("TPSA")
        result["rotatable_bonds"]   = p.get("RotatableBondCount")
        result["charge"]            = p.get("Charge")
        result["complexity"]        = p.get("Complexity")
    except Exception:
        pass

    try:
        # Step 3: Fetch synonyms (top 5)
        syn_url = f"{PUBCHEM_BASE}/compound/cid/{cid}/synonyms/JSON"
        r = requests.get(syn_url, timeout=8)
        r.raise_for_status()
        syns = r.json()["InformationList"]["Information"][0].get("Synonym", [])
        result["synonyms"] = syns[:5]
    except Exception:
        pass

    return result


@app.route("/", methods=["GET", "POST"])
def home():
    output = ""
    img = ""
    molblock = ""
    name = ""
    pubchem = {}
    rdkit_props = {}
    error = ""

    if request.method == "POST":
        name = request.form.get("name", "").strip()
        smiles_input = request.form.get("smiles", "").strip()

        # Fetch from PubChem if name given
        if name:
            pubchem = get_compound_data(name)
            smiles = smiles_input or pubchem.get("smiles")
        else:
            smiles = smiles_input

        if not smiles:
            error = "Could not find compound. Check the name or enter a SMILES string."
            return render_template("index.html", error=error, pubchem=pubchem,
                                   rdkit_props=rdkit_props, output=output,
                                   img=img, molblock=molblock, name=name)

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            error = "Invalid SMILES string — RDKit could not parse the molecule."
            return render_template("index.html", error=error, pubchem=pubchem,
                                   rdkit_props=rdkit_props, output=output,
                                   img=img, molblock=molblock, name=name)

        mol_h = Chem.AddHs(mol)
        Chem.AssignStereochemistry(mol_h, force=True)

        # ── Chirality ──
        centers = Chem.FindMolChiralCenters(mol_h, includeUnassigned=True)
        chiral_lines = []
        if not centers:
            chiral_lines.append("None")
        else:
            for i, c in centers:
                chiral_lines.append(f"Atom {i} → {c}")

        # ── Functional group matches ──
        benzyl  = Chem.MolFromSmarts("[CH2]c1ccccc1")
        phenyl  = Chem.MolFromSmarts("c1ccccc1")
        benzoyl = Chem.MolFromSmarts("C(=O)c1ccccc1")
        hydroxyl = Chem.MolFromSmarts("[OX2H]")
        amine    = Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(NC=O)]")
        carbonyl = Chem.MolFromSmarts("[CX3]=[OX1]")
        carboxyl = Chem.MolFromSmarts("C(=O)[OH]")
        ester    = Chem.MolFromSmarts("C(=O)OC")
        amide    = Chem.MolFromSmarts("C(=O)N")

        groups = {
            "Benzyl":   len(mol_h.GetSubstructMatches(benzyl)),
            "Phenyl":   len(mol_h.GetSubstructMatches(phenyl)),
            "Benzoyl":  len(mol_h.GetSubstructMatches(benzoyl)),
            "Hydroxyl": len(mol_h.GetSubstructMatches(hydroxyl)),
            "Amine":    len(mol_h.GetSubstructMatches(amine)),
            "Carbonyl": len(mol_h.GetSubstructMatches(carbonyl)),
            "Carboxyl": len(mol_h.GetSubstructMatches(carboxyl)),
            "Ester":    len(mol_h.GetSubstructMatches(ester)),
            "Amide":    len(mol_h.GetSubstructMatches(amide)),
        }

        # ── RDKit computed properties ──
        rdkit_props = {
            "Mol. Weight":      round(Descriptors.MolWt(mol), 3),
            "Exact Mass":       round(Descriptors.ExactMolWt(mol), 5),
            "LogP (RDKit)":     round(Descriptors.MolLogP(mol), 3),
            "Heavy Atoms":      mol.GetNumHeavyAtoms(),
            "Ring Count":       Descriptors.RingCount(mol),
            "Aromatic Rings":   Descriptors.NumAromaticRings(mol),
            "Fraction Csp3":    round(Descriptors.FractionCSP3(mol), 3),
            "TPSA (RDKit)":     round(Descriptors.TPSA(mol), 2),
            "H Donors":         Descriptors.NumHDonors(mol),
            "H Acceptors":      Descriptors.NumHAcceptors(mol),
            "Rot. Bonds":       Descriptors.NumRotatableBonds(mol),
            "Stereo Centers":   len(centers),
        }

        # Lipinski rule of five
        lipinski_pass = (
            rdkit_props["Mol. Weight"] <= 500 and
            rdkit_props["LogP (RDKit)"] <= 5 and
            rdkit_props["H Donors"] <= 5 and
            rdkit_props["H Acceptors"] <= 10
        )
        rdkit_props["Lipinski (Ro5)"] = "Pass" if lipinski_pass else "Fail"

        output = {
            "chiral": chiral_lines,
            "groups": groups,
        }

        # ── 2D image (remove explicit H for cleaner drawing) ──
        img_mol = Chem.RemoveHs(mol_h)
        img_data = Draw.MolToImage(img_mol, size=(400, 300))
        buf = io.BytesIO()
        img_data.save(buf, format="PNG")
        img = base64.b64encode(buf.getvalue()).decode()

        # ── 3D structure ──
        AllChem.EmbedMolecule(mol_h, AllChem.ETKDGv3())
        AllChem.UFFOptimizeMolecule(mol_h)
        molblock = Chem.MolToMolBlock(mol_h)

    return render_template("index.html",
                           output=output,
                           img=img,
                           molblock=molblock,
                           name=name,
                           pubchem=pubchem,
                           rdkit_props=rdkit_props,
                           error=error)


if __name__ == "__main__":
    port = int(os.environ.get("PORT", 10000))
    app.run(host="0.0.0.0", port=port)

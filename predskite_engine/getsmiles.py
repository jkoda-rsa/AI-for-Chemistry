import requests

def get_pubchem_smiles(compound_name):
    compound_name = compound_name.replace(' ', '+')
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound_name}/property/CanonicalSMILES/JSON"
    response = requests.get(url)

    if response.status_code == 200:
        data = response.json()
        if "PropertyTable" in data and "Properties" in data["PropertyTable"]:
            smiles = data["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
            return smiles
        else:
            print("No properties found for this compound.")
            return None
    else:
        print("Failed to retrieve data from PubChem.")
        return None


def get_smiles(molecule):
    smiles_dict = {
        'X59': 'COC1=CC=C(C=C1)N(C1=CC2=C(C=C1)C1=CC=C(C=C1C21C2=CC=CC=C2OC=2C=CC=CC12)N(C1=CC=C(C=C1)OC)C1=CC=C(C=C1)OC)C1=CC=C(C=C1)OC',
        'EDOT-OMeTPA': 'N(=CC=1SC(C=NC2=CC=C(C=C2)N(C3=CC=C(OC)C=C3)C4=CC=C(OC)C=C4)=C5OCCOC15)C6=CC=C(C=C6)N(C7=CC=C(OC)C=C7)C8=CC=C(OC)C=C8',
        'BTT-1': 'O(C1=CC=C(C=C1)N(C=2SC=3C(C2)=C4SC(=CC4=C5SC(=CC35)N(C6=CC=C(OC)C=C6)C7=CC=C(OC)C=C7)N(C8=CC=C(OC)C=C8)C9=CC=C(OC)C=C9)C%10=CC=C(OC)C=C%10)C',
        'EDTA': 'O=C(O)CN(CC(=O)O)CCN(CC(=O)O)CC(=O)O',
        'HATNA-F6': 'FC1=CC=2N=C3C(=NC2C=C1F)C=4N=C5C=C(F)C(F)=CC5=NC4C=6N=C7C=C(F)C(F)=CC7=NC36',
        'NDI-P': 'O=C(O)CCN1C(=O)C2=CC=C3C(=O)N(C(=O)C4=CC=C(C1=O)C2=C34)CCC(=O)O'
    }

    if molecule in smiles_dict:
        return smiles_dict[molecule]
    else:
        return get_pubchem_smiles(molecule)
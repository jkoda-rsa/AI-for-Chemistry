from rdkit import Chem
import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem
from getsmiles import get_smiles

def molecule_to_smiles(mol):
    smiles = Chem.MolToSmiles(mol)
    return smiles

def display_3D_structure(molecule):
    smiles = get_smiles(molecule)
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    conf_id = AllChem.EmbedMolecule(mol, AllChem.ETKDG())

    if conf_id != -1:
        AllChem.UFFOptimizeMolecule(mol, confId=conf_id)

        pdb_block = Chem.MolToPDBBlock(mol)

        # Set up the viewer
#        viewer = py3Dmol.view(width=600, height=600)
#        viewer.addModel(pdb_block,'pdb')
#        viewer.setStyle({'stick':{}})
#        viewer.setBackgroundColor('0xeeeeee')
#        viewer.zoomTo()
#        viewer.rotate(720, {'x':1, 'y':1, 'z':1}) # increase rotation speed
#        viewer.spin('x', 1)
#        viewer.show()
    else:
        print(f"Failed to generate a conformer for the molecule: {smiles}")

    return pdb_block



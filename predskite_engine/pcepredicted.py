import numpy as np
from deepchem.feat import RDKitDescriptors
from rdkit import Chem
import pickle
from getsmiles import get_smiles


def process_smiles(smiles_list, labels):
    unique_smiles = []
    unique_labels = []
    
    for i, smiles in enumerate(smiles_list):
        if isinstance(smiles, str):
            molecule = Chem.MolFromSmiles(smiles)
            if molecule is not None:  
                canonical_smiles = Chem.MolToSmiles(molecule, canonical=True)
                if canonical_smiles not in unique_smiles:
                    unique_smiles.append(canonical_smiles)
                    unique_labels.append(labels[i])
    
    return unique_smiles, unique_labels


def featurize_smiles(smiles_list):
    featurizer = RDKitDescriptors()
    features = featurizer.featurize(smiles_list)
    print(f"Number of generated molecular descriptors: {features.shape[1]}")

    features = features[:, ~np.isnan(features).any(axis=0)]
    print(f"Number of molecular descriptors without invalid values: {features.shape[1]}")

    return features


def adjust_features(features, target_features):
    if features.shape[1] > target_features:
        # Use some criteria to discard features
        # Here, just for simplicity, we take the first `target_features` number of features
        features = features[:,:target_features]
    elif features.shape[1] < target_features:
        # In case there are less features, you can pad with zeros (or some other method)
        padding = np.zeros((features.shape[0], target_features - features.shape[1]))
        features = np.hstack((features, padding))

    return features


def estimate_etl_pce(molecule):
    model = pickle.load(open('model/ranf_etl_model.sav', 'rb'))

    smiles_str = get_smiles(molecule)
    processed_smiles = process_smiles([smiles_str], [None])[0]
    features = featurize_smiles(processed_smiles)
    features = adjust_features(features, 197)

    feature_vector = features.reshape(1, -1)

    pce = model.predict(feature_vector)[0]
    pce = round(pce, 3)

    return pce


def estimate_htl_pce(molecule):
    model = pickle.load(open('model/ranf_htl_model.sav', 'rb'))

    smiles_str = get_smiles(molecule)
    processed_smiles = process_smiles([smiles_str], [None])[0]
    features = featurize_smiles(processed_smiles)
    features = adjust_features(features, 197)

    feature_vector = features.reshape(1, -1)

    pce = model.predict(feature_vector)[0]
    pce = round(pce, 3)

    return pce
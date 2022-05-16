from rdkit import Chem
import rdkit.Chem.Draw

# Ejemplos de códigos SMILES
molecules = ['Oc1ccccc1', 'Cc1ccccc1', 'c1c2ccccc2ccc1', 'Cc1c(C)cccc1 ']

for i, mol in enumerate(molecules):

    # Creación de la molécula
    m = Chem.MolFromSmiles(mol)

    # Creación de la imagen
    Chem.Draw.MolToFile(m, f'molecule{i}.png')

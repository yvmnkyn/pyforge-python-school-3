from rdkit import Chem
from rdkit.Chem import Draw

def substructure_search(mols, mol):
    
    substructure = Chem.MolFromSmiles(mol)  # Creating an RDKit molecule object for the substructure
    matching_mols = []  # Creating an empty list to store molecules that match the substructure

    # Iterate over the list of molecule SMILES strings (mols)
    for smiles in mols:
        molecule = Chem.MolFromSmiles(smiles)  # Creating a RDKit molecule object for the current molecule

        if molecule.HasSubstructMatch(substructure):  # Checking if the molecule contains the substructure
            matching_mols.append(smiles)  # If it does, adding the SMILES string to the list of matching molecules
    
    
    return matching_mols # The list matching molecules will be returned



if __name__ == "__main__":
    
    # 2. Basic Molecular Representations
    # Create a molecule from a SMILES string
    smiles = "CCO"  # Ethanol
    molecule = Chem.MolFromSmiles(smiles)

    # Draw the molecule
    Draw.MolToImage(molecule)

    # 3. Substructure Search Example
    # Define the molecule and the substructure to search for
    benzene = Chem.MolFromSmiles("c1ccccc1")
    ethanol = Chem.MolFromSmiles("CCO")

    # Perform the substructure search
    match = ethanol.HasSubstructMatch(benzene)
    print("Benzene ring found in ethanol:", match)

    # 4. Molecular Visualization
    # Create a list of molecules
    smiles_list = ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]
    molecules = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]

    # Draw the molecules in a grid
    img = Draw.MolsToGridImage(molecules, molsPerRow=2, subImgSize=(200, 200), returnPNG=False)
    img.show()
    
    # Example usage of substructure search function
    substructure_smiles = "c1ccccc1"
    matching_molecules = substructure_search(smiles_list, substructure_smiles)
    print(matching_molecules)
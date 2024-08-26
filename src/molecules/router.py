from fastapi import APIRouter, Depends, HTTPException
from molecules.dao import MoleculeDAO
from src.molecules.request_body import RBMolecule, MoleculeUpdate
from src.molecules.schema import MoleculeResponse, MoleculeAdd

router = APIRouter(prefix="/molecules", tags=["molecules"])

@router.get("/", summary="Get all molecules")
async def get_all_molecules(request_body: RBMolecule = Depends()) -> list[MoleculeResponse]:
    return await MoleculeDAO.find_all_molecules(**request_body.to_dict())

@router.get("/{molecule_id}", summary="Get a molecule with an id")
async def get_molecule_by_id(molecule_id: int) -> MoleculeResponse | dict:
    rez = await MoleculeDAO.find_one_or_none_by_id(data_id=molecule_id)
    if rez is None:
        return {"message": f"The molecule with id {molecule_id} is not found!"}
    return rez

@router.post("/add/")
async def add_molecule(molecule: MoleculeAdd) -> dict:
    check = await MoleculeDAO.add_molecule(**molecule.model_dump())
    if check:
        return {"message": "The molecule is added!", "molecule": molecule}
    else:
        return {"message": "Error adding a molecule!"}

@router.put("/{molecule_id}", summary="Update a molecule")
async def update_molecule(molecule_id: int, molecule_update: MoleculeUpdate) -> dict:
    check = await MoleculeDAO.update_molecule(molecule_id, **molecule_update.model_dump())
    if check:
        return {"message": "The molecule is updated!", "molecule": molecule_update}
    else:
        raise HTTPException(status_code=400, detail="Error updating the molecule!")

@router.delete("/delete/{molecule_id}")
async def delete_molecule_by_id(molecule_id: int) -> dict:
    check = await MoleculeDAO.delete_molecule_by_id(molecule_id=molecule_id)
    if check:
        return {"message": f"The molecule with id {molecule_id} is deleted!"}
    else:
        return {"message": "Error deleting a molecule!"}

@router.get("/substructure_search/", summary="Search molecules by substructure")
async def search_substructure(substructure: str) -> list[MoleculeResponse]:
    try:
        matching_molecules = await MoleculeDAO.substructure_search(substructure)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))

    return [MoleculeResponse(identifier=molecule.identifier, smiles=molecule.smiles) for molecule in matching_molecules]

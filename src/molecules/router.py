import logging
from fastapi import APIRouter, Depends, HTTPException, Query
from molecules.dao import MoleculeDAO
from molecules.request_body import RBMolecule, MoleculeUpdate
from molecules.schema import MoleculeResponse, MoleculeAdd

# Initialize logging
logger = logging.getLogger(__name__)

router = APIRouter(prefix="/molecules", tags=["molecules"])

@router.get("/", summary="Get all molecules")
async def get_all_molecules(limit: int = Query(default=100)) -> list[MoleculeResponse]:
    logger.info("Fetching all molecules with a limit of %s", limit)
    molecules = [MoleculeResponse(**molecule) for molecule in await MoleculeDAO.find_all_molecules(limit=limit)]
    logger.info("Fetched %s molecules", len(molecules))
    return molecules

@router.get("/{molecule_id}", summary="Get a molecule with an id")
async def get_molecule_by_id(molecule_id: int) -> MoleculeResponse | dict:
    logger.info("Fetching molecule with ID: %s", molecule_id)
    rez = await MoleculeDAO.find_one_or_none_by_id(data_id=molecule_id)
    if rez is None:
        logger.warning("Molecule with ID %s not found", molecule_id)
        return {"message": f"The molecule with id {molecule_id} is not found!"}
    logger.info("Molecule with ID %s fetched successfully", molecule_id)
    return rez

@router.post("/add/", summary="Add a new molecule")
async def add_molecule(molecule: MoleculeAdd) -> dict:
    logger.info("Adding a new molecule: %s", molecule.identifier)
    check = await MoleculeDAO.add_molecule(**molecule.model_dump())
    if check:
        logger.info("Molecule added successfully: %s", molecule.identifier)
        return {"message": "The molecule is added!", "molecule": molecule}
    else:
        logger.error("Error adding molecule: %s", molecule.identifier)
        return {"message": "Error adding a molecule!"}

@router.put("/{molecule_id}", summary="Update a molecule")
async def update_molecule(molecule_id: int, molecule_update: MoleculeUpdate) -> dict:
    logger.info("Updating molecule with ID: %s", molecule_id)
    check = await MoleculeDAO.update_molecule(molecule_id, **molecule_update.model_dump())
    if check:
        logger.info("Molecule with ID %s updated successfully", molecule_id)
        return {"message": "The molecule is updated!", "molecule": molecule_update}
    else:
        logger.error("Error updating molecule with ID %s", molecule_id)
        raise HTTPException(status_code=400, detail="Error updating the molecule!")

@router.delete("/delete/{molecule_id}", summary="Delete a molecule")
async def delete_molecule_by_id(molecule_id: int) -> dict:
    logger.info("Deleting molecule with ID: %s", molecule_id)
    check = await MoleculeDAO.delete_molecule_by_id(molecule_id=molecule_id)
    if check:
        logger.info("Molecule with ID %s deleted successfully", molecule_id)
        return {"message": f"The molecule with id {molecule_id} is deleted!"}
    else:
        logger.error("Error deleting molecule with ID %s", molecule_id)
        return {"message": "Error deleting a molecule!"}

@router.get("/substructure_search/", summary="Search molecules by substructure")
async def search_substructure(substructure: str) -> list[MoleculeResponse]:
    logger.info("Searching for molecules with substructure: %s", substructure)
    try:
        matching_molecules = await MoleculeDAO.substructure_search(substructure)
    except ValueError as e:
        logger.error("Error during substructure search: %s", e)
        raise HTTPException(status_code=400, detail=str(e))
    
    logger.info("Found %s matching molecules", len(matching_molecules))
    return [MoleculeResponse(identifier=molecule.identifier, smiles=molecule.smiles) for molecule in matching_molecules]

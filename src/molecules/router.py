import logging
import json
from fastapi import APIRouter, Depends, HTTPException, Query, Request
from sqlalchemy.ext.asyncio import AsyncSession
from molecules.dao import MoleculeDAO
from molecules.schema import MoleculeResponse, MoleculeAdd
from molecules.request_body import MoleculeUpdate
from database import get_session

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/molecules", tags=["molecules"])

@router.get("/", summary="Get all molecules")
async def get_molecules(limit: int = 100, session: AsyncSession = Depends(get_session)):
    molecules = await MoleculeDAO.find_all_molecules(session=session, limit=limit)
    return molecules

@router.get("/{molecule_id}", summary="Get a molecule by ID")
async def get_molecule_by_id(molecule_id: int, session: AsyncSession = Depends(get_session)):
    molecule = await MoleculeDAO.find_one_or_none_by_id(session=session, data_id=molecule_id)
    if not molecule:
        raise HTTPException(status_code=404, detail=f"Molecule with ID {molecule_id} not found")
    return molecule

@router.post("/add/", summary="Add a new molecule")
async def add_molecule(molecule: MoleculeAdd, session: AsyncSession = Depends(get_session)) -> dict:
    logger.info("Adding a new molecule: %s", molecule.identifier)
    check = await MoleculeDAO.add_molecule(session=session, **molecule.model_dump())
    if check:
        logger.info("Molecule added successfully: %s", molecule.identifier)
        return {"message": "The molecule is added!", "molecule": molecule}
    else:
        logger.error("Error adding molecule: %s", molecule.identifier)
        return {"message": "Error adding a molecule!"}

@router.put("/{molecule_id}", summary="Update a molecule")
async def update_molecule(molecule_id: int, molecule_update: MoleculeUpdate, session: AsyncSession = Depends(get_session)) -> dict:
    logger.info("Updating molecule with ID: %s", molecule_id)
    check = await MoleculeDAO.update_molecule(session=session, molecule_id=molecule_id, **molecule_update.model_dump())
    if check:
        logger.info("Molecule with ID %s updated successfully", molecule_id)
        return {"message": "The molecule is updated!", "molecule": molecule_update}
    else:
        logger.error("Error updating molecule with ID %s", molecule_id)
        raise HTTPException(status_code=400, detail="Error updating the molecule!")

@router.delete("/delete/{molecule_id}", summary="Delete a molecule")
async def delete_molecule_by_id(molecule_id: int, session: AsyncSession = Depends(get_session)) -> dict:
    logger.info("Deleting molecule with ID: %s", molecule_id)
    check = await MoleculeDAO.delete_molecule_by_id(session=session, molecule_id=molecule_id)
    if check:
        logger.info("Molecule with ID %s deleted successfully", molecule_id)
        return {"message": f"The molecule with id {molecule_id} is deleted!"}
    else:
        logger.error("Error deleting molecule with ID %s", molecule_id)
        return {"message": "Error deleting a molecule!"}

@router.get("/substructure_search/", summary="Search molecules by substructure")
async def search_substructure(substructure: str, request: Request, session: AsyncSession = Depends(get_session)) -> list[MoleculeResponse]:
    logger.info("Searching for molecules with substructure: %s", substructure)
    
    # Get Redis client from the FastAPI app state
    redis_client = request.app.state.redis_client

    # Check Redis cache for cached result
    cache_key = f"substructure_search:{substructure}"
    cached_result = await redis_client.get(cache_key)

    if cached_result:
        logger.info("Returning cached result for substructure: %s", substructure)
        return json.loads(cached_result)

    try:
        matching_molecules = await MoleculeDAO.substructure_search(session=session, substructure=substructure)
    except ValueError as e:
        logger.error("Error during substructure search: %s", e)
        raise HTTPException(status_code=400, detail=str(e))

    # Cache the search result for 60 seconds
    await redis_client.setex(cache_key, 60, json.dumps([m.dict() for m in matching_molecules]))

    logger.info("Found %s matching molecules", len(matching_molecules))
    return [MoleculeResponse(identifier=molecule.identifier, smiles=molecule.smiles) for molecule in matching_molecules]

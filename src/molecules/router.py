import logging
import json
from fastapi import APIRouter, Depends, HTTPException, Query, Request
from sqlalchemy.ext.asyncio import AsyncSession
from celery.result import AsyncResult
from molecules.dao import MoleculeDAO
from molecules.schema import MoleculeResponse, MoleculeAdd
from molecules.request_body import MoleculeUpdate
from database import get_session
from tasks import substructure_search_task

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

# Endpoint to start the substructure search task (Celery task)
@router.post("/substructure_search/", summary="Start a substructure search task")
async def start_substructure_search(substructure_smiles: str):
    task = substructure_search_task.delay(substructure_smiles)
    return {"task_id": task.id}

# Endpoint to check task status and get result
@router.get("/substructure_search/{task_id}/", summary="Get substructure search result")
async def get_substructure_search_result(task_id: str):
    task_result = AsyncResult(task_id)

    if task_result.state == 'PENDING':
        return {"status": "Pending", "result": None}
    elif task_result.state == 'SUCCESS':
        return {"status": "Completed", "result": task_result.result}
    else:
        return {"status": task_result.state, "result": None}

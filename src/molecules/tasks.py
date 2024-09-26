from celery import Celery
from time import sleep
from molecules.dao import MoleculeDAO 

# Initialize Celery app
celery_app = Celery('tasks', broker='redis://localhost:6379/0', backend='redis://localhost:6379/0')

# Celery task for substructure search
@celery_app.task
def substructure_search_task(substructure_smiles: str):
    sleep(5)

    try:
        result = MoleculeDAO.substructure_search(substructure_smiles)
    except ValueError as e:
        result = f"Error in substructure search: {e}"

    return result

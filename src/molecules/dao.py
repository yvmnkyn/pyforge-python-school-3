from sqlalchemy import select
from src.database import async_session_maker
from src.molecules.models import MoleculeModel
from src.dao.base import BaseDAO
from rdkit import Chem

class MoleculeDAO(BaseDAO):
    model = MoleculeModel

    @classmethod
    async def find_all_molecules(cls):
        async with async_session_maker() as session:
            query = select(cls.model)
            result = await session.execute(query)
            return result.scalars().all()

    @classmethod
    async def add_molecule(cls, **molecule_data):
        return await cls.add(**molecule_data)

    @classmethod
    async def delete_molecule_by_id(cls, molecule_id: int):
        return await cls.delete(id=molecule_id)

    @classmethod
    async def update_molecule(cls, molecule_id: int, **molecule_data):
        return await cls.update({"id": molecule_id}, **molecule_data)

    @classmethod
    async def substructure_search(cls, substructure_smiles: str):
        substructure = Chem.MolFromSmiles(substructure_smiles)
        if not substructure:
            raise ValueError("Invalid substructure SMILES")

        async with async_session_maker() as session:
            query = select(cls.model)
            result = await session.execute(query)
            molecules = result.scalars().all()

        matching_ids = []
        for molecule in molecules:
            mol = Chem.MolFromSmiles(molecule.smiles)
            if mol and mol.HasSubstructMatch(substructure):
                matching_ids.append(molecule)

        return matching_ids

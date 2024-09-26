from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession
from molecules.models import MoleculeModel
from dao.base import BaseDAO

class MoleculeDAO:
    
    @classmethod
    async def find_all_molecules(cls, session: AsyncSession, limit: int = 100):
        query = select(MoleculeModel).limit(limit)
        result = await session.execute(query)
        return result.scalars().all()

    @classmethod
    async def find_one_or_none_by_id(cls, session: AsyncSession, data_id: int):
        query = select(MoleculeModel).filter_by(id=data_id)
        result = await session.execute(query)
        return result.scalars().first()

    @classmethod
    async def add_molecule(cls, session: AsyncSession, **molecule_data):
        new_molecule = MoleculeModel(**molecule_data)
        session.add(new_molecule)
        await session.commit()
        await session.refresh(new_molecule)
        return new_molecule

    @classmethod
    async def update_molecule(cls, session: AsyncSession, molecule_id: int, **update_data):
        query = select(MoleculeModel).filter_by(id=molecule_id)
        result = await session.execute(query)
        molecule = result.scalars().first()
        if molecule:
            for key, value in update_data.items():
                setattr(molecule, key, value)
            await session.commit()
            await session.refresh(molecule)
            return molecule
        return None

    @classmethod
    async def delete_molecule_by_id(cls, session: AsyncSession, molecule_id: int):
        query = select(MoleculeModel).filter_by(id=molecule_id)
        result = await session.execute(query)
        molecule = result.scalars().first()
        if molecule:
            await session.delete(molecule)
            await session.commit()
            return True
        return False

    @staticmethod
    async def substructure_search(session: AsyncSession, substructure: str):
        """
        Searches for molecules that contain the given substructure.
        """
        query = select(MoleculeModel).where(MoleculeModel.smiles.contains(substructure))
        result = await session.execute(query)
        return result.scalars().all()
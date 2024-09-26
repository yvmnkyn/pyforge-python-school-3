from sqlalchemy.orm import Mapped
from database import Base, str_uniq, int_pk

class MoleculeModel(Base):
    id: Mapped[int_pk]
    identifier: Mapped[str_uniq]
    smiles: Mapped[str_uniq]

    def __str__(self):
        return (
            f"{self.__class__.__name__}(id={self.id}, "
            f"identifier={self.identifier!r},"
            f"smiles={self.smiles!r})"
        )

    def __repr__(self):
        return str(self)
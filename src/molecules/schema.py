import re
from pydantic import BaseModel, Field, field_validator

class MoleculeResponse(BaseModel):
    id: int
    identifier: str = Field(..., min_length=1, max_length=100, description="Molecule identifier")
    smiles: str = Field(..., min_length=1, max_length=100, description="SMILES string")

    @field_validator("smiles")
    def validate_smiles(cls, value):
        if not re.match(
            r"^(?:[A-Z][a-z]?|[a-z])(?:(?:[1-9]\d*)?(?:\[(?:(?:[A-Z][a-z]?(?:@[@]?)?)"
            r"|[#+-]|\d+)?\])?|(?:[-=#$:/\\])?(?:[A-Z][a-z]?|[a-z])|[().\[\]])*((?:[1-9]\d*)?)$",
            value,
        ):
            raise ValueError("This SMILES has an invalid structure")
        return value

class MoleculeAdd(BaseModel):
    identifier: str = Field(..., min_length=1, max_length=100, description="Molecule identifier")
    smiles: str = Field(..., min_length=1, max_length=100, description="SMILES string")

    @field_validator("smiles")
    def validate_smiles(cls, value):
        if not re.match(
            r"^(?:[A-Z][a-z]?|[a-z])(?:(?:[1-9]\d*)?(?:\[(?:(?:[A-Z][a-z]?(?:@[@]?)?)"
            r"|[#+-]|\d+)?\])?|(?:[-=#$:/\\])?(?:[A-Z][a-z]?|[a-z])|[().\[\]])*((?:[1-9]\d*)?)$",
            value,
        ):
            raise ValueError("This SMILES has an invalid structure")
        return value
    

from pydantic import BaseModel, Field

class RBMolecule(BaseModel):
    molecule_id: int | None = None
    identifier: str | None = None
    smiles: str | None = None

    def to_dict(self) -> dict:
        data = {
            "id": self.molecule_id,
            "identifier": self.identifier,
            "smiles": self.smiles,
        }
        filtered_data = {key: value for key, value in data.items() if value is not None}
        return filtered_data
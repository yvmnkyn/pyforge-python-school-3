from pydantic import BaseModel, Field, field_validator
import re

class RBMolecule:
    def __init__(
        self,
        molecule_id: int | None = None,
        identifier: str | None = None,
        smiles: str | None = None,
    ):
        self.id = molecule_id
        self.identifier = identifier
        self.smiles = smiles

    def to_dict(self) -> dict:
        data = {
            "id": self.id,
            "identifier": self.identifier,
            "smiles": self.smiles,
        }
        filtered_data = {key: value for key, value in data.items() if value is not None}
        return filtered_data


class MoleculeUpdate(BaseModel):
    smiles: str = Field(..., min_length=1, max_length=100, description="Updated SMILES string")

    @field_validator("smiles")
    def validate_smiles(cls, value):
        if not re.match(
            r"^(?:[A-Z][a-z]?|[a-z])(?:(?:[1-9]\d*)?(?:\[(?:(?:[A-Z][a-z]?(?:@[@]?)?)"
            r"|[#+-]|\d+)?\])?|(?:[-=#$:/\\])?(?:[A-Z][a-z]?|[a-z])|[().\[\]])*((?:[1-9]\d*)?)$",
            value,
        ):
            raise ValueError("This SMILES has an invalid structure")
        return value

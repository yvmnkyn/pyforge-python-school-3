from fastapi.testclient import TestClient
from .main import app, molecules

client = TestClient(app)

def test_add_molecule():
    response = client.post("/molecules/", json={"identifier": "mol1", "smiles": "CCO"})
    assert response.status_code == 201
    assert response.json() == {"identifier": "mol1", "smiles": "CCO"}

def test_get_molecule():
    molecules["mol1"] = "CCO"  # Pre-load the molecule
    response = client.get("/molecules/mol1")
    assert response.status_code == 200
    assert response.json() == {"identifier": "mol1", "smiles": "CCO"}

def test_update_molecule():
    molecules["mol1"] = "CCO"
    response = client.put("/molecules/mol1", json={"smiles": "CCN"})
    assert response.status_code == 200
    assert response.json() == {"identifier": "mol1", "smiles": "CCN"}

def test_delete_molecule():
    molecules["mol1"] = "CCO"
    response = client.delete("/molecules/mol1")
    assert response.status_code == 204
    assert "mol1" not in molecules

def test_list_molecules():
    molecules["mol1"] = "CCO"
    response = client.get("/molecules/")
    assert response.status_code == 200
    assert response.json() == [{"identifier": "mol1", "smiles": "CCO"}]

def test_substructure_search():
    molecules["mol1"] = "CCO"
    molecules["mol2"] = "CCN"
    response = client.get("/substructure_search/", params={"substructure": "CC"})
    assert response.status_code == 200
    assert len(response.json()) == 2  # Both molecules should match

def test_invalid_substructure_search():
    response = client.get("/substructure_search/", params={"substructure": "invalid_smiles"})
    assert response.status_code == 400
    assert response.json() == {"detail": "Invalid substructure SMILES"}

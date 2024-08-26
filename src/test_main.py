from fastapi.testclient import TestClient
from main import app

client = TestClient(app)

def test_add_molecule():
    response = client.post(
        "/molecules/", json={"identifier": "mol1", "smiles": "CCO"}
    )
    assert response.status_code == 201
    assert response.json() == {"identifier": "mol1", "smiles": "CCO"}

def test_get_molecule():
    # Assuming your DAO or database setup is correctly populated
    response = client.get("/molecules/1")
    assert response.status_code == 200
    assert response.json() == {"id": 1, "identifier": "mol1", "smiles": "CCO"}

def test_update_molecule():
    response = client.put("/molecules/1", json={"smiles": "CCN"})
    assert response.status_code == 200
    assert response.json() == {"id": 1, "identifier": "mol1", "smiles": "CCN"}

def test_delete_molecule():
    response = client.delete("/molecules/1")
    assert response.status_code == 204

def test_list_molecules():
    response = client.get("/molecules/")
    assert response.status_code == 200
    assert isinstance(response.json(), list)

def test_substructure_search():
    response = client.get(
        "/substructure_search/", params={"substructure": "CC"}
    )
    assert response.status_code == 200
    assert len(response.json()) > 0

def test_invalid_substructure_search():
    response = client.get(
        "/substructure_search/", params={"substructure": "invalid_smiles"}
    )
    assert response.status_code == 400
    assert response.json() == {"detail": "Invalid substructure SMILES"}

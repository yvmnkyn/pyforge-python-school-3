from fastapi.testclient import TestClient
from main import app

client = TestClient(app)

def test_add_molecule():
    response = client.post(
        "/molecules/add/", json={"identifier": "mol1", "smiles": "CCO"}
    )
    assert response.status_code == 200
    assert response.json()["molecule"]["identifier"] == "mol1"

def test_get_molecule():
    response = client.get("/molecules/1")
    assert response.status_code == 200
    assert response.json()["identifier"] == "mol1"

def test_update_molecule():
    response = client.put("/molecules/1", json={"smiles": "CCN"})
    assert response.status_code == 200
    assert response.json()["molecule"]["smiles"] == "CCN"

def test_delete_molecule():
    response = client.delete("/molecules/delete/1")
    assert response.status_code == 200
    assert response.json()["message"] == "The molecule with id 1 is deleted!"

def test_list_molecules():
    response = client.get("/molecules/")
    assert response.status_code == 200
    assert isinstance(response.json(), list)

def test_substructure_search():
    response = client.get(
        "/molecules/substructure_search/", params={"substructure": "CC"}
    )
    assert response.status_code == 200
    assert len(response.json()) > 0

def test_invalid_substructure_search():
    response = client.get(
        "/molecules/substructure_search/", params={"substructure": "invalid_smiles"}
    )
    assert response.status_code == 400

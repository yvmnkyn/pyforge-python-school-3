import pytest
from httpx import AsyncClient
from redis.asyncio import Redis
from redis.exceptions import ConnectionError
from main import app
import fakeredis

@pytest.fixture(scope="module")
async def async_client():
    async with AsyncClient(app=app, base_url="http://test") as ac:
        yield ac

@pytest.fixture(scope="module")
async def mock_redis():
    redis_client = fakeredis.aioredis.FakeRedis()
    yield redis_client
    await redis_client.close()

@pytest.mark.asyncio
async def test_add_molecule(async_client):
    response = await async_client.post("/molecules/add/", json={"identifier": "mol1", "smiles": "CCO"})
    assert response.status_code == 200
    assert response.json()["molecule"]["identifier"] == "mol1"

@pytest.mark.asyncio
async def test_get_molecule(async_client):
    response = await async_client.get("/molecules/1")
    assert response.status_code == 200
    assert response.json()["identifier"] == "mol1"

@pytest.mark.asyncio
async def test_update_molecule(async_client):
    response = await async_client.put("/molecules/1", json={"smiles": "CCN"})
    assert response.status_code == 200
    assert response.json()["molecule"]["smiles"] == "CCN"

@pytest.mark.asyncio
async def test_delete_molecule(async_client):
    response = await async_client.delete("/molecules/delete/1")
    assert response.status_code == 200
    assert response.json()["message"] == "The molecule with id 1 is deleted!"

@pytest.mark.asyncio
async def test_list_molecules(async_client):
    response = await async_client.get("/molecules/")
    assert response.status_code == 200
    assert isinstance(response.json(), list)

@pytest.mark.asyncio
async def test_substructure_search(async_client, mock_redis):
    # First request should not be cached
    response = await async_client.get("/molecules/substructure_search/", params={"substructure": "CC"})
    assert response.status_code == 200
    assert len(response.json()) > 0

    # Ensure result is cached
    response = await async_client.get("/molecules/substructure_search/", params={"substructure": "CC"})
    assert response.status_code == 200
    assert len(response.json()) > 0

@pytest.mark.asyncio
async def test_invalid_substructure_search(async_client):
    response = await async_client.get("/molecules/substructure_search/", params={"substructure": "invalid_smiles"})
    assert response.status_code == 400

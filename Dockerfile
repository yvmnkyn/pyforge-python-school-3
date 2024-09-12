# Use Miniconda base image
FROM continuumio/miniconda3

# Set the working directory
WORKDIR /app

# Add PYTHONPATH to ensure that /app is recognized as a module path
ENV PYTHONPATH=/app

# Install necessary dependencies
RUN pip install fastapi uvicorn[standard] pydantic sqlalchemy alembic asyncpg pydantic-settings redis
RUN conda install -c conda-forge rdkit
RUN pip install pytest httpx
RUN pip install celery redis

# Copy application code to /app
COPY ./src /app

# Command to run the FastAPI application
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]

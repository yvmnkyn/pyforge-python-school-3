# Use Miniconda base image
FROM continuumio/miniconda3

# Set the working directory
WORKDIR /app

# Add PYTHONPATH to ensure that /app is recognized as a module path
ENV PYTHONPATH=/app

# Install necessary dependencies (does not change often)
RUN pip install fastapi uvicorn[standard] pydantic sqlalchemy alembic asyncpg
RUN conda install -c conda-forge rdkit

# Copy application code (changes more frequently)
COPY ./src /app

# Debug: List contents of /app to verify structure
RUN ls -la /app

# Command to run the FastAPI application
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]

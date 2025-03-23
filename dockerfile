# Use a base Python image with RDKit pre-installed
FROM continuumio/miniconda3:latest

# Set working directory
WORKDIR /app

# Copy environment.yml to container
COPY environment.yml .

# Create conda environment
RUN conda env create -f environment.yml

# Make RUN commands use the new environment
SHELL ["conda", "run", "-n", "feplanner", "/bin/bash", "-c"]

# Copy application code and templates
COPY app.py .
COPY templates/ templates/

# Expose port for Flask
EXPOSE 5001

# Run Flask application
CMD ["conda", "run", "--no-capture-output", "-n", "feplanner", "python", "app.py"]
# FEPlanner

A Flask-based web service that processes SDF files and plans FEP+ maps using the OpenFE library with Lomap atom mapping.

## Features

- Upload SDF files through a REST API or web interface
- Generate FEP+ maps using Lomap atom mapper
- Visualize the resulting network with D3.js
- Configure mapping parameters (3D information, distance cutoff, element changes)
- Dockerized application for easy deployment

## Prerequisites

To run the application directly, you'll need:

- Python 3.9+
- RDKit
- OpenFE
- Flask

Alternatively, you can use Docker to run the application without installing the dependencies.

## Installation

### Using Docker (recommended)

1. Build the Docker image:

```bash
docker build -t fep-mapping-service .
```

2. Run the container:

```bash
docker run -p 5000:5000 fep-mapping-service
```

The web service will be available at http://localhost:5000

### Manual Installation

1. Create a Conda environment with the required dependencies:

```bash
conda create -n openfe -c conda-forge python=3.9 rdkit openmmforcefields openmm openfe
conda activate openfe
pip install flask
```

2. Create the necessary directory structure:

```bash
mkdir -p templates
```

3. Save the Flask application as `app.py`
4. Save the HTML template as `templates/index.html`
5. Run the application:

```bash
python app.py
```

## API Usage

### Health Check

```
GET /health
```

Returns a simple health check response to verify the service is running.

### Generate FEP+ Map

```
POST /plan-fep-map
```

Upload an SDF file to generate an FEP+ map.

**Parameters:**

- `file`: SDF file containing molecules (required)
- `threed`: Boolean to use 3D information for mapping (default: true)
- `max3d`: Float for maximum 3D distance for atom mapping (default: 1.0)
- `element_change`: Boolean to allow changes in atom elements (default: false)

**Example using curl:**

```bash
curl -X POST -F "file=@your_molecules.sdf" -F "threed=true" -F "max3d=1.0" -F "element_change=false" http://localhost:5000/plan-fep-map
```

**Response:**

```json
{
  "status": "success",
  "network": {
    "nodes": [
      {
        "name": "Mol1",
        "num_atoms": 23,
        "formula": "C14H9NO"
      },
      ...
    ],
    "edges": [
      {
        "molecule_a": "Mol1",
        "molecule_b": "Mol2",
        "mapping": {"0": 0, "1": 1, ...},
        "score": 0.95
      },
      ...
    ]
  }
}
```

## Web Interface

The web interface is available at the root URL (http://localhost:5000). It provides a user-friendly way to:

1. Upload SDF files
2. Configure mapping parameters
3. Visualize the resulting network
4. Inspect edge details by clicking on connections

## Development

The application structure is as follows:

- `app.py`: The main Flask application with API endpoints and processing logic
- `templates/index.html`: The web interface with D3.js visualization
- `Dockerfile`: Container configuration for deployment

To extend the functionality:

1. Modify `app.py` to add new API endpoints or processing features
2. Update `templates/index.html` to enhance the web interface
3. Rebuild the Docker image if using containerized deployment

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- [OpenFE](https://docs.openfree.energy/): Free energy calculation library
- [RDKit](https://www.rdkit.org/): Cheminformatics and machine learning toolkit
- [Flask](https://flask.palletsprojects.com/): Web framework
- [D3.js](https://d3js.org/): Data visualization library

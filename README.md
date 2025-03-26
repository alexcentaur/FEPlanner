# FEP+ Mapping Service

A web application for generating FEP+ maps from SDF files. This service helps plan free energy perturbation calculations by creating networks of molecules with optimal atom mappings.

## Features

- Upload SDF files containing multiple molecules
- Generate different types of FEP networks:
  - Minimal Spanning Network
  - Minimal Redundant Network
  - Radial Network (with user-selected center ligand)
- Interactive network visualization with:
  - Drag-and-drop node positioning
  - Zoom and pan capabilities
  - Edge score visualization
  - Atom mapping details on edge click
- Manual node connection capability
- Export network data in JSON format
- Compound structure visualization
- Network summary statistics

## Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/FEPlanner.git
cd FEPlanner
```

2. Create a virtual environment and activate it:
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install the required packages:
```bash
pip install -r requirements.txt
```

## Usage

1. Start the web server:
```bash
python app.py
```

2. Open your web browser and navigate to:
```
http://localhost:5000
```

3. Upload your SDF file and configure the mapping parameters:
   - Use 3D information for mapping (optional)
   - Set maximum 3D distance for atom mapping
   - Allow changes in atom elements (optional)
   - Choose network type (Minimal Spanning, Minimal Redundant, or Radial)
   - For Radial networks, select a center ligand

4. Click "Generate FEP+ Map" to create the network

5. Interact with the results:
   - View the network graph
   - Examine compound structures
   - Check network summary
   - Download network data as JSON
   - Create manual connections between nodes

## Network Types

### Minimal Spanning Network
- Creates a network with the minimum number of edges needed to connect all molecules
- Optimizes for efficiency and reduces computational cost
- Best for initial screening or when computational resources are limited

### Minimal Redundant Network
- Includes additional edges beyond the minimal spanning network
- Provides alternative paths between molecules
- More robust to failures and allows for cross-validation
- Better for production runs where accuracy is crucial

### Radial Network
- Creates a hub-and-spoke topology with a user-selected center ligand
- Useful for analyzing relationships between a central molecule and others
- Helps identify key molecules in the network
- Ideal for lead optimization or when studying a specific molecule's interactions

## Requirements

- Python 3.8 or higher
- OpenFE library
- RDKit
- Flask
- vis.js (included via CDN)


import os
import tempfile
import json
import base64
from flask import Flask, request, jsonify, send_file, render_template, session
from werkzeug.utils import secure_filename
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import openfe
from openfe.setup import LomapAtomMapper
from openfe.setup.ligand_network_planning import generate_minimal_spanning_network

app = Flask(__name__)
app.secret_key = os.urandom(24)  # For session management

# Configure upload folder for temporary storage
UPLOAD_FOLDER = tempfile.mkdtemp()
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024  # 16 megabytes

# Configure allowed file extensions
ALLOWED_EXTENSIONS = {'sdf'}

# Global variable to store SDF data temporarily for the current session
sdf_cache = {}

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/health', methods=['GET'])
def health_check():
    """Simple health check endpoint to verify service is running."""
    return jsonify({
        'status': 'OK',
        'message': 'FEP+ Mapping Service is running'
    })

@app.route('/plan-fep-map', methods=['POST'])
def plan_fep_map():
    """
    Endpoint to process an SDF file and generate an FEP+ map.
    
    Expects a multipart/form-data POST request with:
    - 'file': SDF file containing molecules
    
    Optional parameters:
    - threed: boolean (default: True) - use 3D information for mapping
    - max3d: float (default: 1.0) - maximum 3D distance for atom mapping
    - element_change: boolean (default: False) - allow changes in atom elements
    """
    # Check if file was provided in request
    if 'file' not in request.files:
        return jsonify({
            'status': 'error',
            'message': 'No file part in the request'
        }), 400
    
    file = request.files['file']
    
    # Check if file was selected
    if file.filename == '':
        return jsonify({
            'status': 'error',
            'message': 'No file selected'
        }), 400
    
    # Check if file is allowed
    if not allowed_file(file.filename):
        return jsonify({
            'status': 'error',
            'message': f'File type not supported. Allowed types: {", ".join(ALLOWED_EXTENSIONS)}'
        }), 400
    
    # Securely save the file
    filename = secure_filename(file.filename)
    filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
    file.save(filepath)
    
    try:
        # Extract parameters from request
        threed = request.form.get('threed', 'true').lower() == 'true'
        max3d = float(request.form.get('max3d', 1.0))
        element_change = request.form.get('element_change', 'false').lower() == 'true'
        
        print(f"Processing file with parameters: threed={threed}, max3d={max3d}, element_change={element_change}")
        
        # Save SDF content for molecule rendering
        try:
            # Try to read as text first
            with open(filepath, 'r') as f:
                sdf_content = f.read()
        except UnicodeDecodeError:
            # If that fails, try reading as binary and decode
            with open(filepath, 'rb') as f:
                binary_content = f.read()
                sdf_content = binary_content.decode('utf-8', errors='replace')
        
        print(f"Read SDF file with {sdf_content.count('$$$$')} molecule blocks")
        
        # Generate a unique ID for this SDF file
        import uuid
        sdf_id = str(uuid.uuid4())
        
        # Store in the cache
        sdf_cache[sdf_id] = sdf_content
        print(f"Stored SDF in cache with ID: {sdf_id}")
        
        # Process file and generate FEP+ map
        result = process_sdf_file(filepath, threed, max3d, element_change)
        
        # Add the SDF ID to the result
        result['sdf_id'] = sdf_id
        
        # Clean up
        os.remove(filepath)
        
        return jsonify(result)
    
    except Exception as e:
        # Get detailed error information including traceback
        import traceback
        error_traceback = traceback.format_exc()
        print(f"Error processing file: {str(e)}")
        print(f"Traceback: {error_traceback}")
        
        # Clean up in case of error
        if os.path.exists(filepath):
            os.remove(filepath)
        
        return jsonify({
            'status': 'error',
            'message': f'Error processing file: {str(e)}',
            'traceback': error_traceback
        }), 500

@app.route('/get-sdf/<sdf_id>', methods=['GET'])
def get_sdf(sdf_id):
    """Return the SDF content for a given ID."""
    if sdf_id not in sdf_cache:
        return jsonify({
            'status': 'error',
            'message': 'SDF file not found. It may have expired.'
        }), 404
        
    return jsonify({
        'status': 'success',
        'sdf_content': sdf_cache[sdf_id]
    })

@app.route('/molecule-svg/<sdf_id>/<int:mol_index>', methods=['GET'])
def molecule_svg(sdf_id, mol_index):
    """Generate and return an SVG image for a specific molecule."""
    if sdf_id not in sdf_cache:
        print(f"SDF ID {sdf_id} not found in cache. Available IDs: {list(sdf_cache.keys())}")
        return jsonify({
            'status': 'error',
            'message': 'SDF file not found'
        }), 404
    
    try:
        sdf_content = sdf_cache[sdf_id]
        print(f"Retrieved SDF content for ID {sdf_id}, length: {len(sdf_content)}")
        
        # Split SDF content into individual molecules
        molecules = []
        current_mol = ""
        
        for line in sdf_content.splitlines():
            current_mol += line + "\n"
            if line.strip() == "$$$$":
                if current_mol.strip():  # Only add non-empty molecules
                    molecules.append(current_mol)
                current_mol = ""
                
        # Add the last molecule if it doesn't end with $$$$
        if current_mol.strip() and "$$$$" not in current_mol:
            current_mol += "$$$$\n"
            molecules.append(current_mol)
        
        print(f"Found {len(molecules)} molecules in SDF")
        
        # Check if the requested molecule index is valid
        if mol_index < 0 or mol_index >= len(molecules):
            print(f"Molecule index {mol_index} out of range (0-{len(molecules)-1})")
            return jsonify({
                'status': 'error',
                'message': f'Molecule index {mol_index} out of range (0-{len(molecules)-1})'
            }), 404
            
        mol_block = molecules[mol_index]
        print(f"Molecule block length: {len(mol_block)}")
        
        # Convert MolBlock to RDKit molecule
        mol = Chem.MolFromMolBlock(mol_block)
        if mol is None:
            print(f"Failed to parse molecule block: {mol_block[:100]}...")
            return jsonify({
                'status': 'error',
                'message': 'Failed to parse molecule'
            }), 500
        
        # Generate SVG
        from rdkit.Chem import Draw
        from rdkit.Chem.Draw import rdMolDraw2D
        drawer = rdMolDraw2D.MolDraw2DSVG(300, 300)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        
        return svg, 200, {'Content-Type': 'image/svg+xml'}
        
    except Exception as e:
        import traceback
        error_trace = traceback.format_exc()
        print(f"Error generating SVG: {str(e)}")
        print(f"Traceback: {error_trace}")
        return jsonify({
            'status': 'error',
            'message': f'Error generating SVG: {str(e)}',
            'traceback': error_trace
        }), 500

@app.route('/molecule-placeholder/<name>/<formula>', methods=['GET'])
def molecule_placeholder(name, formula):
    """Generate a placeholder SVG when molecule rendering fails."""
    svg = f'''<?xml version="1.0" encoding="UTF-8"?>
    <svg xmlns="http://www.w3.org/2000/svg" width="300" height="300" viewBox="0 0 300 300">
        <rect width="300" height="300" fill="#f9f9f9" stroke="#cccccc" stroke-width="1"/>
        <text x="150" y="100" font-family="Arial" font-size="16" text-anchor="middle" fill="#333333">
            {name}
        </text>
        <text x="150" y="130" font-family="Arial" font-size="14" text-anchor="middle" fill="#666666">
            {formula}
        </text>
        <text x="150" y="170" font-family="Arial" font-size="12" text-anchor="middle" fill="#999999">
            Structure rendering not available
        </text>
    </svg>'''
    return svg, 200, {'Content-Type': 'image/svg+xml'}

def process_sdf_file(filepath, threed=True, max3d=1.0, element_change=False):
    """
    Process an SDF file and generate an FEP+ map using Lomap atom mapper.
    
    Args:
        filepath: Path to the SDF file
        threed: Use 3D positions for mapping
        max3d: Maximum distance between atoms for mapping
        element_change: Allow element changes in mapping
    
    Returns:
        Dictionary containing the FEP+ mapping results
    """
    # Load molecules from SDF file
    ligands = []
    rdkit_mols = []  # Store the RDKit molecules
    
    for mol in Chem.SDMolSupplier(filepath, removeHs=False):
        if mol is not None:
            # Create SmallMoleculeComponent for each molecule
            ligand = openfe.SmallMoleculeComponent(mol)
            ligands.append(ligand)
            rdkit_mols.append(mol)
    
    if not ligands:
        raise ValueError("No valid molecules found in the SDF file")
    
    # Create the atom mapper with specified parameters
    mapper = LomapAtomMapper(
        threed=threed,
        max3d=max3d,
        element_change=element_change
    )
    
    # Debug print the mapper
    print(f"Mapper created: {mapper}")
    
    try:
        # Generate the network using minimal spanning tree approach
        network = generate_minimal_spanning_network(
            ligands=ligands,
            scorer=openfe.lomap_scorers.default_lomap_score,
            mappers=[mapper]
        )
    except Exception as e:
        print(f"Error generating network: {str(e)}")
        import traceback
        print(traceback.format_exc())
        raise ValueError(f"Failed to generate network: {str(e)}")

    print(f"Network created with {len(network.edges)} edges")
    
    # Extract edges and mappings from the network
    edges = []
    for i, edge in enumerate(network.edges):
        print(edge)
        mol_a = edge.componentA.name
        mol_b = edge.componentB.name
        
        
        #print(f"Processing edge {i}: {mol_a} -> {mol_b}")
        #print(f"Edge type: {type(edge).__name__}")
        #print(f"Edge attributes: {dir(edge)}")
        
        # For older versions of OpenFE (before 0.15.0), the edge itself may be the mapping
        try:
            # Try various attribute paths to find the mapping
            if hasattr(edge, 'componentA_to_componentB'):
                # Try to access the mapping directly from the edge
                print("Using edge.componentA_to_componentB")
                mapping = edge.componentA_to_componentB
            elif hasattr(edge, 'mapping') and hasattr(edge.mapping, 'componentA_to_componentB'):
                print("Using edge.mapping.componentA_to_componentB")
                mapping = edge.mapping.componentA_to_componentB
            elif hasattr(edge, 'atom_mapping') and hasattr(edge.atom_mapping, 'componentA_to_componentB'):
                print("Using edge.atom_mapping.componentA_to_componentB")
                mapping = edge.atom_mapping.componentA_to_componentB
            elif hasattr(edge, 'transformation') and hasattr(edge.transformation, 'componentA_to_componentB'):
                print("Using edge.transformation.componentA_to_componentB")
                mapping = edge.transformation.componentA_to_componentB
            else:
                # Default to an empty dict if we can't find the mapping
                print(f"Could not find mapping path. Available edge attributes: {dir(edge)}")
                if hasattr(edge, 'atom_mapping'):
                    print(f"atom_mapping attributes: {dir(edge.atom_mapping)}")
                mapping = {}
        except Exception as e:
            print(f"Error accessing mapping: {str(e)}")
            mapping = {}
        
        # Convert mapping (dict with int keys) to string keys for JSON serialization
        serializable_mapping = {str(k): v for k, v in mapping.items()}
        
        # Try to calculate the score
        try:
            # Check for score in annotations dictionary
            if hasattr(edge, 'annotations') and 'score' in edge.annotations:
                print("Using edge.annotations['score']")
                score = edge.annotations['score']
            # For newer versions of OpenFE, score may be a method or property on the edge
            elif hasattr(edge, 'score'):
                print("Using edge.score")
                if callable(edge.score):
                    score = edge.score()
                else:
                    score = edge.score
            # Otherwise try to calculate it using the lomap scorer
            elif hasattr(edge, 'mapping'):
                print("Calculating score with edge.mapping")
                score = openfe.lomap_scorers.default_lomap_score(edge.mapping)
            elif hasattr(edge, 'atom_mapping'):
                print("Calculating score with edge.atom_mapping")
                score = openfe.lomap_scorers.default_lomap_score(edge.atom_mapping)
            else:
                # Default to a placeholder score
                print("Using default score of 0.5")
                score = 0.5
        except Exception as e:
            print(f"Error calculating score: {str(e)}")
            score = 0.5
        
        edges.append({
            'molecule_a': mol_a,
            'molecule_b': mol_b,
            'mapping': serializable_mapping,
            'score': score
        })
    
    # Extract nodes (molecules) from the network
    nodes = []
    for i, ligand in enumerate(ligands):
        # Calculate molecular formula from atom counts
        rdmol = ligand.to_rdkit()
        atom_dict = {}
        for atom in rdmol.GetAtoms():
            symbol = atom.GetSymbol()
            atom_dict[symbol] = atom_dict.get(symbol, 0) + 1
        
        # Format the molecular formula
        formula = ''
        for symbol in sorted(atom_dict.keys()):
            count = atom_dict[symbol]
            if count == 1:
                formula += symbol
            else:
                formula += f"{symbol}{count}"
        
        # Generate SMILES string
        try:
            smiles = Chem.MolToSmiles(rdkit_mols[i])
        except:
            # If original molecule isn't available, use regenerated one
            smiles = Chem.MolToSmiles(rdmol)
                
        nodes.append({
            'name': ligand.name,
            'num_atoms': rdmol.GetNumAtoms(),
            'formula': formula,
            'smiles': smiles
        })
    
    # Return the network as a JSON-serializable dictionary
    return {
        'status': 'success',
        'network': {
            'nodes': nodes,
            'edges': edges
        }
    }

@app.route('/', methods=['GET'])
def index():
    """Serve the main web interface for users to upload files."""
    return render_template('index.html')

@app.route('/molecule-svg-from-smiles', methods=['POST'])
def molecule_svg_from_smiles():
    """Generate and return an SVG image for a molecule from a SMILES string."""
    try:
        data = request.get_json()
        if not data or 'smiles' not in data:
            return jsonify({
                'status': 'error',
                'message': 'No SMILES provided in request'
            }), 400
        
        smiles = data['smiles']
        print(f"Generating SVG for SMILES: {smiles}")
        
        # Optional parameters
        width = data.get('width', 300)
        height = data.get('height', 300)
        
        # Convert SMILES to RDKit molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({
                'status': 'error',
                'message': 'Failed to parse SMILES string'
            }), 400
        
        # Generate 2D coordinates
        mol = Chem.AddHs(mol)
        AllChem.Compute2DCoords(mol)
        
        # Generate SVG
        drawer = Draw.rdMolDraw2D.MolDraw2DSVG(width, height)
        drawer.SetFontSize(0.8)  # Slightly smaller font for better fit
        
        # Set draw options
        opts = drawer.drawOptions()
        opts.addAtomIndices = False
        opts.addStereoAnnotation = True
        opts.atomHighlightsAreCircles = True
        opts.additionalAtomLabelPadding = 0.15  # Add some padding around atom labels
        
        # Draw the molecule
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        
        return svg, 200, {'Content-Type': 'image/svg+xml'}
        
    except Exception as e:
        import traceback
        error_trace = traceback.format_exc()
        print(f"Error generating SVG from SMILES: {str(e)}")
        print(f"Traceback: {error_trace}")
        return jsonify({
            'status': 'error',
            'message': f'Error generating SVG: {str(e)}',
            'traceback': error_trace
        }), 500

@app.route('/molecule-svg-from-smiles/<smiles>', methods=['GET'])
def molecule_svg_from_smiles_get(smiles):
    """Generate and return an SVG image for a molecule from a SMILES string using GET method."""
    try:
        print(f"Generating SVG for SMILES: {smiles}")
        
        # Parse parameters from query string
        width = request.args.get('width', default=300, type=int)
        height = request.args.get('height', default=300, type=int)
        
        # Convert SMILES to RDKit molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({
                'status': 'error',
                'message': 'Failed to parse SMILES string'
            }), 400
        
        # Generate 2D coordinates
        mol = Chem.AddHs(mol)
        AllChem.Compute2DCoords(mol)
        
        # Generate SVG
        drawer = Draw.rdMolDraw2D.MolDraw2DSVG(width, height)
        drawer.SetFontSize(0.8)  # Slightly smaller font for better fit
        
        # Set draw options
        opts = drawer.drawOptions()
        opts.addAtomIndices = False
        opts.addStereoAnnotation = True
        opts.atomHighlightsAreCircles = True
        opts.additionalAtomLabelPadding = 0.15  # Add some padding around atom labels
        
        # Draw the molecule
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        
        return svg, 200, {'Content-Type': 'image/svg+xml'}
        
    except Exception as e:
        import traceback
        error_trace = traceback.format_exc()
        print(f"Error generating SVG from SMILES: {str(e)}")
        print(f"Traceback: {error_trace}")
        return jsonify({
            'status': 'error',
            'message': f'Error generating SVG: {str(e)}',
            'traceback': error_trace
        }), 500

if __name__ == '__main__':
    print("Starting FEP+ Mapping Service...")
    app.run(host='0.0.0.0', port=5001, debug=False)
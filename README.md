# GCH1 Multi-Tumor Network Visualization

## Overview
This project creates an interactive network visualization with GCH1 as the central node, connected to multiple tumor types and their associated genes.

## Features
- Interactive network visualization with GCH1 as central node
- Multiple tumor types displayed as second-level nodes
- Tumor-specific genes connected to their respective tumor types
- Cross-tumor genes (appearing in multiple tumors) shown with special diamond shape
- Node size and edge width scaled by correlation strength
- Interactive controls for clustering, expanding, and resetting layout
- Position freezing functionality to fix node positions
- PNG download options (regular and high-resolution)
- Responsive design with hover tooltips for additional information

## Requirements
- Python 3.6+
- Required Python packages:
  - pandas
  - numpy
  - pyvis
  - pathlib

## Usage
1. Place your tumor correlation data CSV files in the `data` folder
2. Run the visualization script:
   ```
   python multi_tumor_network.py
   ```
3. Open the generated HTML file in a web browser to explore the network

## Data Format
Input CSV files should contain the following columns:
- Gene Symbol: The gene identifier
- PCC: Pearson Correlation Coefficient with GCH1

## Output
The script generates an interactive HTML visualization in the `output` directory.

## Customization
You can modify various aspects of the visualization by editing parameters in the `multi_tumor_network.py` file:
- Color schemes for tumor types
- Correlation thresholds
- Physics parameters for layout
- Maximum genes per tumor type

## License
MIT 
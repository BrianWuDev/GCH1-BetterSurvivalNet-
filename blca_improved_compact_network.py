#!/usr/bin/env python
"""
BLCA Improved Compact Network Visualization with Enhanced Download

This script creates an optimized network visualization with even more compact node layout
and enhanced PNG download options (including high-resolution export) for BLCA tumors.
"""

import pandas as pd
import numpy as np
import logging
import os
import re
from pathlib import Path
import json
import time
from typing import Dict, List, Tuple, Optional, Union
from pyvis.network import Network
from PIL import Image

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger('BLCA_Improved_Compact_Network')


class ImprovedCompactNetworkVisualizer:
    """Class for creating optimized BLCA gene network with ultra-compact layout and enhanced PNG download."""
    
    def __init__(self):
        """Initialize the visualizer with optimized settings."""
        # Settings
        self.tumor_types = ['BLCA Tumor']
        self.min_correlation = 0.4
        self.max_genes_per_tumor = 25  # Further reduced to make visualization more compact
        
        # Node colors - professional color palette with improved contrast
        self.central_node_color = '#FF4136'  # Bright red for GCH1
        self.tumor_colors = {
            'BLCA Tumor': '#3D9970'  # Green
        }
        
        # Visual parameters - optimized for ultra-compact layout
        self.node_size_range = (5, 12)  # Smaller nodes for more compact display
        self.edge_width_range = (0.3, 2)  # Thinner edges
        
        # Output settings
        self.output_dir = 'output'
        self.output_filename = 'blca_improved_compact_network'
        self.title = 'BLCA Tumor with Decreased Survival Rates due to GCH1 Upregulation'
        
        # Create output directory
        self.output_dir = Path(self.output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Initialize statistics tracking
        self.gene_counts = {tumor: 0 for tumor in self.tumor_types}
        
        # Physics optimization parameters
        self.physics_params = {
            "gravity": -100000,  # Stronger gravity to pull nodes closer
            "central_gravity": 1.0,  # Much stronger central gravity
            "spring_length": 80,  # Very short springs
            "spring_strength": 0.15,  # Stronger springs
            "damping": 0.3,  # Moderate damping
            "avoid_overlap": 0.2  # Minimal overlap avoidance
        }

    def load_data(self, file_path: str = '../Network0501/data/tumor.csv') -> pd.DataFrame:
        """Load and filter tumor data from CSV file.
        
        Args:
            file_path: Path to the tumor data CSV
            
        Returns:
            DataFrame with filtered data
        """
        try:
            df = pd.read_csv(file_path)
            logger.info(f"Successfully loaded {len(df)} rows from {file_path}")
            
            # Filter for specified tumor types
            filtered_df = df[df['Tumor'].isin(self.tumor_types)]
            logger.info(f"Filtered to {len(filtered_df)} rows for selected tumor types")
            
            # Apply correlation threshold
            high_corr_df = filtered_df[filtered_df['PCC'] >= self.min_correlation]
            logger.info(f"Applied correlation threshold ≥ {self.min_correlation}, resulting in {len(high_corr_df)} rows")
            
            # Print gene counts by tumor type
            logger.info(f"Genes with correlation ≥ {self.min_correlation} by tumor type:")
            for tumor in self.tumor_types:
                count = len(high_corr_df[high_corr_df['Tumor'] == tumor])
                logger.info(f"  {tumor}: {count} genes")
                print(f"  {tumor}: {count} genes")
            
            self.high_corr_df = high_corr_df
            return high_corr_df
            
        except Exception as e:
            logger.error(f"Error loading data: {e}")
            raise RuntimeError(f"Failed to load data: {e}")

    def create_network(self) -> Network:
        """Create the ultra-compact PyVis network with nodes and edges.
        
        Returns:
            PyVis Network object
        """
        # Create network with appropriate configuration
        net = Network(
            height='800px', 
            width='1000px', 
            bgcolor='#ffffff', 
            font_color='black', 
            directed=False,
            notebook=False,
            heading='',
            cdn_resources='remote'  # Use remote resources to avoid encoding issues
        )
        
        # Configure network display options - ultra-compact layout
        net.toggle_hide_edges_on_drag(False)
        net.toggle_physics(True)
        net.toggle_drag_nodes(True)
        net.toggle_stabilization(True)
        net.barnes_hut(
            gravity=self.physics_params["gravity"],
            central_gravity=self.physics_params["central_gravity"],
            spring_length=self.physics_params["spring_length"],
            spring_strength=self.physics_params["spring_strength"],
            damping=self.physics_params["damping"],
        )
        
        # Add GCH1 as center node
        net.add_node(
            'GCH1', 
            label='GCH1', 
            size=20,  # Smaller than original
            color=self.central_node_color, 
            title='GCH1 (Central Gene)', 
            shape='dot', 
            borderWidth=2, 
            font={'size': 18, 'face': 'Arial', 'color': 'white', 'strokeWidth': 3, 'strokeColor': '#000000'}
        )
        
        # Add BLCA tumor node radially
        tumor = 'BLCA Tumor'
        # Add tumor node
        net.add_node(
            tumor, 
            label=tumor, 
            size=15,  # Smaller tumor nodes
            color=self.tumor_colors[tumor],
            title=f"{tumor}", 
            shape='dot',
            borderWidth=2,
            font={'size': 14, 'face': 'Arial', 'color': 'white', 'strokeWidth': 2, 'strokeColor': '#000000'},
            x=0,
            y=300
        )
        
        # Connect tumor node to GCH1
        net.add_edge('GCH1', tumor, width=2, color='rgba(150,150,150,0.8)')
        
        # Get genes for this tumor type, sorted by correlation
        tumor_genes = self.high_corr_df[self.high_corr_df['Tumor'] == tumor]
        tumor_genes = tumor_genes.sort_values('PCC', ascending=False)
        
        # Add top genes - position them closer to their tumor node
        # using a spiral pattern around each tumor node
        spiral_factor = 0.2
        
        for gene_idx, (_, gene_row) in enumerate(tumor_genes.head(self.max_genes_per_tumor).iterrows()):
            gene_id = gene_row['Gene Symbol']
            pcc = gene_row['PCC']
            
            # Calculate node size based on correlation
            size_min, size_max = self.node_size_range
            size = size_min + (pcc - self.min_correlation) * (size_max - size_min) / (1 - self.min_correlation)
            
            # Calculate edge width based on correlation
            width_min, width_max = self.edge_width_range
            width = width_min + (pcc - self.min_correlation) * (width_max - width_min)
            
            # Use tumor color with transparency
            color_rgb = self.tumor_colors[tumor]
            
            # Calculate spiral position around the BLCA tumor
            spiral_angle = 2 * np.pi * gene_idx / 8  # 8 genes per revolution
            spiral_distance = 50 + spiral_factor * gene_idx * 5  # Increasing distance
            gene_x = 0 + spiral_distance * np.cos(spiral_angle)
            gene_y = 300 + spiral_distance * np.sin(spiral_angle)
            
            # Add tooltip with details
            tooltip = f"{gene_id}<br>PCC: {pcc:.3f}<br>Tumor: {tumor}"
            
            net.add_node(
                gene_id, 
                label=gene_id, 
                title=tooltip,
                size=size, 
                color={'background': color_rgb, 'border': color_rgb},
                shape='dot',
                borderWidth=1,
                font={'size': 9},  # Smaller font
                x=gene_x,
                y=gene_y
            )
            
            # Connect gene to its tumor type
            net.add_edge(
                tumor, 
                gene_id, 
                width=width, 
                title=f"PCC: {pcc:.3f}",
                color={'color': color_rgb, 'opacity': 0.6}
            )
            
            self.gene_counts[tumor] += 1
        
        # Log statistics
        logger.info("Genes added to network:")
        print("\nGenes added to network:")
        for tumor, count in self.gene_counts.items():
            logger.info(f"  {tumor}: {count} genes")
            print(f"  {tumor}: {count} genes")
        
        # Configure physics for ultra-compact layout
        physics_options = {
            "enabled": True,
            "solver": "forceAtlas2Based",
            "forceAtlas2Based": {
                "gravitationalConstant": -150,  # Stronger gravitational pull
                "centralGravity": 0.15,  # More central gravity
                "springLength": 80,  # Shorter springs
                "springConstant": 0.25,  # Stronger springs
                "damping": 0.3,
                "avoidOverlap": 0.2  # Reduced overlap avoidance
            },
            "stabilization": {
                "enabled": True,
                "iterations": 1500,  # More iterations for better stabilization
                "updateInterval": 25,
                "fit": True
            }
        }
        
        # Configure options for ultra-compact layout
        options = {
            "nodes": {
                "borderWidth": 1,
                "borderWidthSelected": 2,
                "opacity": 0.9,
                "font": {
                    "face": "Arial",
                    "size": 11,  # Smaller font
                    "strokeWidth": 1,
                    "strokeColor": "#ffffff"
                }
            },
            "edges": {
                "color": {
                    "inherit": False,
                    "opacity": 0.7
                },
                "smooth": {
                    "enabled": True,
                    "type": "continuous",
                    "roundness": 0.5
                },
                "arrows": {
                    "to": {"enabled": False},
                    "from": {"enabled": False}
                }
            },
            "physics": physics_options,
            "interaction": {
                "dragNodes": True,
                "dragView": True,
                "zoomView": True,
                "navigationButtons": False,
                "keyboard": {
                    "enabled": True,
                    "speed": {"x": 10, "y": 10, "zoom": 0.1},
                    "bindToWindow": True
                }
            },
            "configure": {
                "enabled": False
            }
        }
        
        net.set_options(json.dumps(options))
        return net

    def create_html_with_title_and_legend(self, html_path: str, net: Network) -> str:
        """Create HTML with title, legend, and enhanced download buttons.
        
        Args:
            html_path: Path to save the HTML file
            net: PyVis Network object
            
        Returns:
            Path to the modified HTML file
        """
        # Save basic network HTML
        net.save_graph(html_path)
        
        # Read the HTML content with explicit utf-8 encoding
        with open(html_path, 'r', encoding='utf-8') as file:
            html_content = file.read()
        
        # Remove any horizontal rules or dividers in the original HTML
        html_content = html_content.replace('<hr>', '')
        html_content = html_content.replace('<hr/>', '')
        html_content = html_content.replace('<hr />', '')
        
        # Add our CSS to the head
        css_styles = '''
        <style>
            /* Force no borders and dividers */
            body, html, div, canvas, section, article, header, footer, aside, main, nav {
                border: none !important;
                outline: none !important;
            }
            
            /* Remove any border between elements */
            div + div, h1 + div, p + div, div + canvas {
                border-top: none !important;
                margin-top: 0 !important;
                padding-top: 0 !important;
            }
            
            /* Container styles */
            .network-container {
                max-width: 1200px;
                width: 100%;
                margin: 0 auto;
                padding: 0;
                display: flex;
                flex-direction: column;
                align-items: center;
                border: none !important;
                background-color: #f9f9f9;
            }
            
            /* Title styles */
            .header-container {
                text-align: center;
                margin-bottom: 20px;
                width: 100%;
                border-bottom: none !important;
                padding-bottom: 0 !important;
                box-shadow: none !important;
            }
            
            .main-title {
                color: #333;
                font-size: 24px;
                margin-bottom: 10px;
            }
            
            .subtitle {
                color: #666;
                font-size: 16px;
            }
            
            /* Network styles */
            #mynetwork {
                width: 1000px !important;
                height: 800px !important;
                margin: 0 auto !important;
                border: none !important;
                outline: none !important;
                box-shadow: none !important;
            }
            
            /* Legend styles */
            .legend-container {
                position: absolute;
                top: 20px;
                left: 20px;
                background-color: white;
                padding: 15px;
                border-radius: 5px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                z-index: 1000;
                max-width: 220px;
                font-size: 14px;
            }
            
            .legend-title {
                text-align: center;
                margin-bottom: 10px;
                font-weight: bold;
            }
            
            .legend-item {
                display: flex;
                align-items: center;
                margin-bottom: 5px;
            }
            
            .legend-color {
                width: 15px;
                height: 15px;
                border-radius: 50%;
                margin-right: 10px;
            }
            
            .legend-section {
                margin-top: 10px;
                font-weight: bold;
            }
            
            .legend-info {
                margin-top: 10px;
                font-size: 12px;
            }
            
            .legend-controls {
                margin-top: 10px;
                border-top: 1px solid #eee;
                padding-top: 10px;
                font-size: 12px;
            }
            
            .controls-title {
                font-weight: bold;
                margin-bottom: 5px;
            }
            
            /* Control buttons */
            .control-buttons {
                position: absolute;
                bottom: 20px;
                left: 20px;
                z-index: 1000;
                display: flex;
            }
            
            .control-buttons button {
                padding: 8px 12px;
                color: white;
                border: none;
                border-radius: 5px;
                cursor: pointer;
                margin-right: 10px;
                font-weight: bold;
                transition: all 0.2s ease;
            }
            
            #cluster-btn {
                background-color: #9C27B0;
            }
            
            #expand-btn {
                background-color: #FF9800;
            }
            
            #reset-btn {
                background-color: #607D8B;
            }
            
            /* Download buttons */
            .download-buttons {
                position: absolute;
                top: 20px;
                right: 20px;
                z-index: 1000;
                display: flex;
            }
            
            #download-btn, #download-hires-btn {
                padding: 8px 12px;
                color: white;
                border: none;
                border-radius: 5px;
                cursor: pointer;
                font-weight: bold;
                transition: all 0.2s ease;
            }
            
            #download-btn {
                background-color: #4CAF50;
                margin-right: 10px;
            }
            
            #download-hires-btn {
                background-color: #2196F3;
            }
            
            #download-btn:hover, #download-hires-btn:hover,
            .control-buttons button:hover {
                filter: brightness(1.1);
                box-shadow: 0 2px 5px rgba(0,0,0,0.2);
            }
            
            /* Hide divider elements and controls */
            hr, .vis-separator, .separator, .divider, 
            .vis-navigation, .vis-button, .vis-toolbar {
                display: none !important;
                height: 0 !important;
                margin: 0 !important;
                padding: 0 !important;
                border: none !important;
            }
        </style>
        '''
        
        # Create title HTML
        title_html = f'''
        <div class="header-container">
            <h1 class="main-title">{self.title}</h1>
            <p class="subtitle">Visualization of GCH1 gene correlations with BLCA tumor and associated genes</p>
        </div>
        '''
        
        # Create legend HTML
        legend_html = f'''
        <div class="legend-container">
            <div class="legend-title">Legend</div>
            <div class="legend-item">
                <div class="legend-color" style="background-color: {self.central_node_color};"></div>
                <span>GCH1 (Central)</span>
            </div>
            
            <div class="legend-section">Cancer Type:</div>
        '''
        
        # Add tumor type to the legend
        tumor = 'BLCA Tumor'
        legend_html += f'''
        <div class="legend-item">
            <div class="legend-color" style="background-color: {self.tumor_colors[tumor]};"></div>
            <span>{tumor}</span>
        </div>
        '''
        
        # Add correlation info to legend
        legend_html += f'''
            <div class="legend-info">
                <div>Node size: PCC correlation strength</div>
                <div>Edge width: Connection strength</div>
                <div>Minimum correlation: {self.min_correlation}</div>
            </div>
            
            <div class="legend-controls">
                <div class="controls-title">Interactive Controls:</div>
                <div>- Drag nodes to rearrange</div>
                <div>- Scroll to zoom in/out</div>
                <div>- Hover over nodes for details</div>
                <div>- Use download buttons for PNG</div>
            </div>
        </div>
        '''
        
        # Add control buttons
        control_buttons = '''
        <div class="control-buttons">
            <button id="cluster-btn">Cluster Nodes</button>
            <button id="expand-btn">Expand Nodes</button>
            <button id="reset-btn">Reset Layout</button>
        </div>
        '''
        
        # Add JavaScript for control buttons
        control_script = '''
        <script>
        // Get the network instance
        let network = null;
        
        // Wait for the network to be created
        document.addEventListener('DOMContentLoaded', function() {
            // The network object is stored in the global window object by vis.js
            // Wait a bit for the network to initialize
            setTimeout(function() {
                try {
                    if (window.network) {
                        network = window.network;
                        console.log("Network object found");
                    }
                } catch (e) {
                    console.error("Error accessing network:", e);
                }
            }, 1000);
        });
        
        // Cluster nodes
        document.getElementById('cluster-btn').addEventListener('click', function() {
            if (network) {
                // Increase gravity to pull nodes together
                network.physics.options.forceAtlas2Based.gravitationalConstant = -500;
                network.physics.options.forceAtlas2Based.centralGravity = 0.4;
                network.physics.options.forceAtlas2Based.springLength = 50;
                network.physics.options.forceAtlas2Based.springConstant = 0.3;
                
                // Update physics and stabilize
                network.setOptions({physics: network.physics.options});
                network.stabilize(100);
            }
        });
        
        // Expand nodes
        document.getElementById('expand-btn').addEventListener('click', function() {
            if (network) {
                // Decrease gravity to push nodes apart
                network.physics.options.forceAtlas2Based.gravitationalConstant = -50;
                network.physics.options.forceAtlas2Based.centralGravity = 0.05;
                network.physics.options.forceAtlas2Based.springLength = 200;
                network.physics.options.forceAtlas2Based.springConstant = 0.05;
                
                // Update physics and stabilize
                network.setOptions({physics: network.physics.options});
                network.stabilize(100);
            }
        });
        
        // Reset layout
        document.getElementById('reset-btn').addEventListener('click', function() {
            if (network) {
                // Reset to original physics settings
                network.physics.options.forceAtlas2Based.gravitationalConstant = -150;
                network.physics.options.forceAtlas2Based.centralGravity = 0.15;
                network.physics.options.forceAtlas2Based.springLength = 80;
                network.physics.options.forceAtlas2Based.springConstant = 0.25;
                network.physics.options.forceAtlas2Based.damping = 0.3;
                network.physics.options.forceAtlas2Based.avoidOverlap = 0.2;
                
                // Update physics and stabilize
                network.setOptions({physics: network.physics.options});
                network.stabilize(300);
                
                // Reset zoom
                network.fit();
            }
        });
        </script>
        '''
        
        # Add html2canvas library and download buttons with script
        html2canvas_lib = '''
        <script src="https://html2canvas.hertzen.com/dist/html2canvas.min.js"></script>
        '''
        
        # Create download buttons
        download_buttons = '''
        <div class="download-buttons">
            <button id="download-btn">Download PNG</button>
            <button id="download-hires-btn">Download Hi-Res PNG</button>
        </div>
        '''
        
        # Add JavaScript for downloading PNG with different resolutions
        download_script = '''
        <script>
        // Regular PNG download
        document.getElementById('download-btn').addEventListener('click', function() {
            downloadNetworkImage();
        });
        
        // High-resolution PNG download
        document.getElementById('download-hires-btn').addEventListener('click', function() {
            downloadNetworkImage(3.0); // 3x higher resolution
        });
        
        function downloadNetworkImage(scaleFactor = 1.0) {
            // Get the entire document for capturing, including legend and title
            const fullPage = document.querySelector('body');
            
            if (fullPage) {
                // Hide the download buttons during capture
                const downloadBtns = document.querySelector('.download-buttons');
                const controlBtns = document.querySelector('.control-buttons');
                
                if (downloadBtns) downloadBtns.style.display = 'none';
                if (controlBtns) controlBtns.style.display = 'none';
                
                // Use html2canvas to capture the visualization including legend and title
                html2canvas(fullPage, {
                    scale: scaleFactor,
                    backgroundColor: '#f9f9f9',
                    logging: false,
                    allowTaint: true,
                    useCORS: true
                }).then(canvas => {
                    // Create a temporary link element
                    const link = document.createElement('a');
                    link.href = canvas.toDataURL('image/png');
                    
                    // Set filename based on resolution
                    const filename = scaleFactor > 1.0 ? 'blca_network_hires.png' : 'blca_network.png';
                    link.download = filename;
                    
                    // Append, click, and remove the link
                    document.body.appendChild(link);
                    link.click();
                    document.body.removeChild(link);
                    
                    // Restore download button visibility
                    if (downloadBtns) downloadBtns.style.display = 'flex';
                    if (controlBtns) controlBtns.style.display = 'flex';
                });
            } else {
                alert('Could not find the visualization. Please wait until it is fully loaded.');
            }
        }
        </script>
        '''
        
        # Find the position of </head> to inject CSS
        head_end = html_content.find('</head>')
        if head_end != -1:
            html_content = html_content[:head_end] + css_styles + html2canvas_lib + html_content[head_end:]
        
        # Find the body tag to wrap network in our custom container
        body_start = html_content.find('<body')
        if body_start != -1:
            body_tag_end = html_content.find('>', body_start)
            if body_tag_end != -1:
                # Insert the start of container and title after opening body tag
                html_content = html_content[:body_tag_end + 1] + '\n<div class="network-container">\n' + title_html + html_content[body_tag_end + 1:]
        
        # Find where to insert the closing container tag
        body_end = html_content.find('</body>')
        if body_end != -1:
            # Add close of container div before body end
            html_content = html_content[:body_end] + '\n</div>\n' + legend_html + control_buttons + download_buttons + control_script + download_script + html_content[body_end:]
        
        # Write modified HTML back to file
        with open(html_path, 'w', encoding='utf-8', errors='ignore') as file:
            file.write(html_content)
        
        return html_path

    def run(self, data_path: str = '../Network0501/data/tumor.csv') -> str:
        """Run the visualization process.
        
        Args:
            data_path: Path to the tumor data CSV
            
        Returns:
            Path to the generated HTML file
        """
        start_time = time.time()
        
        # Load and filter data
        self.load_data(data_path)
        
        # Create ultra-compact network
        net = self.create_network()
        
        # Create HTML file with title, legend, and enhanced download buttons
        html_path = os.path.join(self.output_dir, f'{self.output_filename}.html')
        self.create_html_with_title_and_legend(html_path, net)
        
        elapsed_time = time.time() - start_time
        logger.info(f"Visualization created in {elapsed_time:.2f} seconds")
        logger.info(f"Interactive visualization with enhanced features created at: {html_path}")
        
        print(f"\nInteractive visualization created at: {html_path}")
        print("Enhanced Features:")
        print("- Ultra-compact node layout")
        print("- Regular and high-resolution PNG download options")
        print("- Node clustering controls (bottom-left)")
        print("- Interactive exploration tools")
        print(f"- Visualization includes {sum(self.gene_counts.values())} genes across {len(self.tumor_types)} tumor types")
        
        return html_path


if __name__ == "__main__":
    visualizer = ImprovedCompactNetworkVisualizer()
    html_file = visualizer.run()
    
    # Open the file in the default browser
    import webbrowser
    try:
        webbrowser.open('file://' + os.path.abspath(html_file))
        print("\nOpening visualization in web browser...")
    except Exception as e:
        print(f"\nCould not open browser automatically: {e}")
        print(f"Please open the file manually: {html_file}") 
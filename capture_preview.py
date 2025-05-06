#!/usr/bin/env python
"""
Preview Image Capture Script for GCH1 Network Visualization

This script automates the capturing of a preview image from the network visualization
for use in the GitHub README and other documentation.
"""

import os
import time
import logging
from pathlib import Path
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger('Preview_Capture')

def capture_preview(html_file: str = 'output/multi_tumor_network.html', output_path: str = 'preview.png'):
    """Capture a preview image of the network visualization.
    
    Args:
        html_file: Path to the HTML file containing the network visualization
        output_path: Path where to save the preview image
    """
    logger.info(f"Starting capture of preview from {html_file}")
    
    # Set up Chrome options
    chrome_options = Options()
    chrome_options.add_argument("--headless")  # Run in headless mode
    chrome_options.add_argument("--window-size=1600,1200")  # Set window size
    chrome_options.add_argument("--disable-gpu")  # Disable GPU acceleration
    chrome_options.add_argument("--no-sandbox")  # Required for some environments
    
    try:
        # Check if HTML file exists
        if not os.path.exists(html_file):
            logger.error(f"HTML file not found: {html_file}")
            return False
        
        # Get absolute path to HTML file
        abs_path = os.path.abspath(html_file)
        file_url = f"file://{abs_path}"
        
        # Initialize Chrome driver
        logger.info("Initializing Chrome driver")
        driver = webdriver.Chrome(options=chrome_options)
        
        # Load the HTML file
        logger.info(f"Loading URL: {file_url}")
        driver.get(file_url)
        
        # Wait for network visualization to load and stabilize
        logger.info("Waiting for visualization to stabilize")
        time.sleep(5)  # Give time for the visualization to render
        
        # Take screenshot
        logger.info(f"Capturing screenshot to {output_path}")
        driver.save_screenshot(output_path)
        
        # Close browser
        driver.quit()
        logger.info("Preview capture completed successfully")
        return True
        
    except Exception as e:
        logger.error(f"Error capturing preview: {e}")
        return False

def main():
    # Make sure output directory exists
    Path('output').mkdir(exist_ok=True)
    
    # Capture preview image
    success = capture_preview()
    
    if success:
        print(f"\nPreview image captured successfully and saved to 'preview.png'")
        print("You can now use this image in your GitHub README.")
    else:
        print("\nFailed to capture preview image. Check logs for details.")

if __name__ == "__main__":
    main() 
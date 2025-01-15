# Cell State Transition Simulation

This MATLAB script simulates the transition of cell states in a biological system, focusing on the progression from normal cells to cancer cells, complex cells, and necrotic cells.

## Features

- Reads cell data from CSV files
- Calculates cell areas
- Simulates cell state transitions based on probabilistic rules
- Visualizes the cell state changes over time
- Provides statistical output for each time step

## Prerequisites

- MATLAB (version R2019b or later recommended)
- CSV files containing cell data:
  - `body_elems.csv`: Cell connectivity information
  - `body_point.csv`: Coordinates of cell vertices
  - `body_esuel.csv`: Neighbors of each cell

## Usage

1. Ensure all required CSV files are in the same directory as the script.
2. Open the script in MATLAB.
3. Run the script.

## Simulation Parameters

- `p1`: Probability for a normal cell to become cancerous
- `p2`: Probability for a cancer cell to become complex
- `p3`: Probability for a complex cell to become necrotic
- `numIterations`: Number of time steps to simulate
- `cancerCell`: Index of the first cancer cell

## Cell States

- 0: Normal cell (Green)
- 1: Cancer cell (Red)
- 2: Complex cell (Blue)
- 3: Necrotic cell (Gray)

## Visualization

The script generates a 3D surface plot that updates at each time step, showing the progression of cell states. The color of each cell represents its current state.

## Output

For each time step, the script prints the count of cells in each state:


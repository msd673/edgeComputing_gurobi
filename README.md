## Overview
This project explores offloading SPARQL query processing to edge servers to overcome latency and bandwidth limitations associated with cloud-based RDF data management solutions. By integrating edge computing techniques, RDF graph data processing is moved closer to the edge, significantly enhancing query performance.

## Description
This project implements an optimization model for query assignment and computational resource allocation in an edge computing system. The system involves multiple terminal devices (EUs) and edge servers (ESs) connected through high-bandwidth networks. The primary goal is to minimize the query response time by optimally assigning queries to edge servers and allocating computational resources based on the available bandwidth and computational capabilities.

The optimization problem is solved using the branch and bound method and the Gurobi solver is used for efficient solving.

## Features
- Query Assignment Optimization: Assigns queries to edge servers based on bandwidth and computational resources.
- Branch-and-Bound Algorithm: Explores all possible allocations to find the optimal solution.
- Gurobi Solver: Solves the optimization problem using Mixed-Integer Linear Programming (MILP).

## Requirements
- CMake: To configure and build the project.
- GUROBI: Required for optimization tasks in the project.
- C++ Compiler: Make sure your system has a working C++ compiler.

## Installation
1. Download GUROBI
   To use the optimization model, you need to install GUROBI. Follow the steps below to download and set it up:
   - Go to the GUROBI website and download the appropriate version for your operating system.
   - Install GUROBI following the installation instructions provided on the website.
   - Set the GUROBI_HOME environment variable to the GUROBI installation path:
2. Build the Project
   Once GUROBI is set up, you can proceed with building the project. Follow these steps:
   ```bash
   git clone <url>
   
   cd edgeComputing_gurobi
   
   make
   
3. Running the Application
   After building the project, run the compiled binary:
   ```bash
   ./gurobi_EC

# Contact
For any questions or inquiries, please contact us at msd673@hnu.edu.cn. 

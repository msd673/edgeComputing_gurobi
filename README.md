## Overview
This project explores offloading SPARQL query processing to edge servers to overcome latency and bandwidth limitations associated with cloud-based RDF data management solutions. By integrating edge computing techniques, RDF graph data processing is moved closer to the edge, significantly enhancing query performance.

## Description
This project implements an optimization model for query assignment and computational resource allocation in an edge computing system. The system involves multiple terminal devices (EUs) and edge servers (ESs) connected through high-bandwidth networks. The primary goal is to minimize the query response time by optimally assigning queries to edge servers and allocating computational resources based on the available bandwidth and computational capabilities.

The optimization problem is solved using the branch and bound method and the Gurobi solver is used for efficient solving.

## Requirements
- CMake: To configure and build the project.
- GUROBI: Required for optimization tasks in the project.
- C++ Compiler: Make sure your system has a working C++ compiler.

## Installation
### 1. Download GUROBI
To use the optimization model, you need to install GUROBI. Follow the steps below to download and set it up:
- Go to the GUROBI website and download the appropriate version for your operating system.
- Install GUROBI following the installation instructions provided on the website.
- Set the GUROBI_HOME environment variable to the GUROBI installation path. 
### 2. Build the Project
Once GUROBI is set up, you can proceed with building the project. Follow these steps:
   ```bash
   git clone <url>
   
   cd edgeComputing_gurobi
   
   make
   ```
### 3. Input Parameters
When executing the binary, the user will be prompted to provide the following parameters in sequence:
#### System Scale Parameters
- Number of EUs (`n`): the total number of end users.
- Number of ESs (`k`): the total number of edge servers.
#### Network Configuration
Servers should be configured with the `iperf` tool to measure the available bandwidth between the servers and end users.
- Cloud server IP: IP address of the cloud server. In our experimental setup, this can be an AWS cloud instance.
- Edge server IPs: IP addresses of all edge servers, provided sequentially.
#### Input Data Files
- Query executable vector (`e`): an `n × k` matrix file.
- Computation amount vector of queries (`c`): a vector of size `n`.
- Result size vector of queries (`w`): a vector of size `n`.
- Computational capability vector of edge servers (`F`): a vector of size `k`.
### 4. Running the Application
After building the project, run the compiled binary:
   ```bash
   ./gurobi_EC
   ```
# Contact
For any questions or inquiries, please contact us at msd673@hnu.edu.cn. 

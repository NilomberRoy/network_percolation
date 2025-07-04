# MDA_percolation

## Overview

This program simulates **network percolation** using the **Mediated Degree Attachment (MDA) model**. It builds a network by preferentially attaching new nodes via mediators and then runs a percolation process using the Achlioptas product rule. The code computes various physical quantities (order parameter, entropy, specific heat, susceptibility, Binder cumulant, etc.) and outputs them for further analysis.

## Features

- **Network Construction:**  
  Builds a connected seed network, then iteratively grows the network by adding nodes using the MDA process.

- **Percolation Process:**  
  Simulates bond percolation using the Achlioptas process (product rule), tracking the formation and growth of clusters.
  
# ER_percolation
This program simulates a percolation process on a random network of N nodes, where at each step, M random pairs of nodes are considered and the pair whose merging minimizes the product of their cluster sizes is connected (the product rule, related to explosive percolation). The simulation tracks the growth of clusters, specifically monitoring the size of the largest cluster, entropy, and other statistical properties as edges are added, thereby modeling the critical behavior and phase transition in network connectivity. Results are output for further analysis of quantities such as entropy, susceptibility, and largest cluster size during the percolation process.

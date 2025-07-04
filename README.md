# MDA_percolation

## Overview

This program simulates **network percolation** using the **Mediated Degree Attachment (MDA) model**. It builds a network by preferentially attaching new nodes via mediators and then runs a percolation process using the Achlioptas product rule. The code computes various physical quantities (order parameter, entropy, specific heat, susceptibility, Binder cumulant, etc.) and outputs them for further analysis.

## Features

- **Network Construction:**  
  Builds a connected seed network, then iteratively grows the network by adding nodes using the MDA process.

- **Percolation Process:**  
  Simulates bond percolation using the Achlioptas process (product rule), tracking the formation and growth of clusters.

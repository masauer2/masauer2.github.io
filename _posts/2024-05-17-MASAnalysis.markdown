---
layout: post
title:  "An Introduction to MASAnalysis"
date:   2024-05-17 17:23:56 -0700
categories: jekyll update
---
This is a guide to using the MASAnalysis package found at https://github.com/masauer2/MASContactCalculator. 

## Setup

```
cd MASContactCalculator
pip install -e .
```
Currently only requiring numpy for array usage and pytest for unittests.

# Example Usage

Note: These instructions are made available in `scripts/fromMAS.py`, which will read the dcd file `data/R1.dcd` and compute the distance matrix for all atoms over the 500 frame trajectory.

## Step 1 - Setting up the DCD Reader
Instantiate DCDReader object with the trajectory filename. Will read in header frame.<br/> 
```
dcd = DCDReader("trajectory.dcd")
```

## Step 2 - Using the DCD Reader to read from file

Create an array to hold empty Frame objects. Each Frame will store the atomic coordinates of the system at a specific timestep. For each frame in the trajectory, we will use the read_DCD_Frame() function to read the data for the next timestep. <br/>
```
frames = np.empty(len(dcd), dtype=object)
for frameNum in range(options["nRead"]):
  frames[frameNum] = dcd.read_DCD_Frame()
```

## Step 3 - Calculate Molecular Properties

After reading in the trajectory, we can calculate the molecular properties of each Frame. <br/>

For example, we can compute the contact distance matrix. The Frame class includes static methods to calculate the distances between all atoms in a given Frame. <br/>

```
for frameNum in range(len(dcd)):
  distance_matrix = Frame.compute_distance_matrix(frames[frameNum], frames[frameNum])
  Frame.output_distance_matrix(distance_matrix, "matrix.out")
```

The distance matrix can be calculated for a subset of atomic coordinates in the system. Selections can be made using the frame.get_selection(arr_slice) function. <br/>

The following code computes the distance matrix between atoms 0-19 and atoms 50-249 for all timesteps in the trajectory:

```
sel1 = np.arange(0,20,1)
sel2 = np.arange(50,250,1)
for frameNum in range(len(dcd)):
  distance_matrix = Frame.compute_distance_matrix(frames[frameNum].get_selection(sel1), frames[frameNum].get_selection(sel2))
  Frame.output_distance_matrix(distance_matrix, "matrix.out")
```


# To-Do List
- Data types -- the user should have to make selections with a predetermined format - not array slices. (implement Selection class)
- One class per type of calculation i.e distance calculations should be stored in a `Distance` class w/ static methods to compute properties on a frame.
- Topology Reader
- Other file formats? Trr? XYZ/GRO/PDB?

[jekyll-docs]: https://jekyllrb.com/docs/home
[jekyll-gh]:   https://github.com/jekyll/jekyll
[jekyll-talk]: https://talk.jekyllrb.com/

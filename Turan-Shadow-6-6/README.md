## getting started
To run this code, you need one more piece
clone the repo `https://github.com/snap-stanford/snap` under `SNAP_directory`    
The two files to add/change are provided in this repo under `SNAP_directory/snap/examples/node2vec`    
Then `cd SNAP_directory/snap/examples/node2vec` and run `make`   

## sample run
check `sample_run.jl` or `sample_run.ipynb` for a few datasets ready to run
For access to data, contact me.

## The main function of interest is:
`visualize_all(filename,A,displayedges;labels=[],TSfunction=x->x,fromk=3,tok=50)`
- `filename` is the prefix of the filename for this dataset
- `A` is the adjacency matrix
- `displayedges` is true if you want to show edges, false if you want to show just a scatter plot

Note that the last two three methods are dependent on compiling a shared object. Instructions on the Node2Vec code is provided here. If you need to run the ECC method, contact me to see if I can share the C code (it's not my code). Otherwise comment out the last two methods (lines 64-81).

A dropbox folder to save all figures: https://www.dropbox.com/sh/ltlsl4pf5ee7toz/AAAYALVzHy9Ww4cZMYXOWiyra?dl=0


Base.run(`
    for school in "American75.mat" "Amherst41.mat" ; do \
            /p/mnt/software/julia-1.2.0/bin/julia run_fb100data.jl $$school ; \
    done`)
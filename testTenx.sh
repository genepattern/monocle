rm -rf job_1/*

docker run  -v ${PWD}:$PWD -w $PWD/job_1 genepattern/monocle Rscript /monocle/run_monocle3.R --input.10x.file $PWD/GpUnit/data/pbmc3k_filtered_gene_bc_matrices.tar.gz --output.file test1 --max_dim 3 --resolution 2.0 --reduction yes

rm -rf job_1/*

docker run  -v ${PWD}:$PWD -w $PWD/job_1 genepattern/monocle Rscript /monocle/run_monocle3.R --input.seurat.rds.file $PWD/GpUnit/data/mat_regressed.Robj --output.file test1 --root.cells spheroid_ACTGTGAAGACCAAAT  --plot.color.by cluster,orig.ident,Phase --plot.post.color.by orig.ident,Phase 

~                                                                
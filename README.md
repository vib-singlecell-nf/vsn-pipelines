# SCENIC

## Running the pipeline

### Generate the config file

- Single SCENIC run

*Note*: The `qsub` profile if you are not running the pipeline on a cluster.

```{bash}
nextflow config \
   -profile scenic,qsub,singularity vib-singlecell-nf/vib-singlecell-nf \
   > nextflow.config
```

- Multi-runs SCENIC

*Note*: The `qsub` profile if you are not running the pipeline on a cluster.

```{bash}
nextflow config \
   -profile scenic_multiruns,qsub,singularity vib-singlecell-nf/vib-singlecell-nf \
   > nextflow.config
```

### Update the config file

Make sure the following parameters are correctly set:
- `params.global.project_name`
- `params.global.qsubaccount` if running on a cluster (SGE cluster)
- `params.sc.scenic.filteredLoom`
- `params.sc.scenic.grn.TFs`
- `params.sc.scenic.cistarget.mtfDB`
- `params.sc.scenic.cistarget.mtfANN`
- `params.sc.scenic.cistarget.trkDB` if commented, track-based cisTarget won't run
- `params.sc.scenic.cistarget.trkDB` if commented, track-based cisTarget won't run
- `params.sc.scenic.numRuns` if running SCENIC in multi-runs mode
- `singularity.runOptions` Specify the paths to mount
- `params.sc.scope.tree`

Additionally, you can update the other paraemeters for the different steps.

### Run 

```{bash}
nextflow -C nextflow.config run \
   vib-singlecell-nf/vib-singlecell-nf \
      -entry scenic \
      -with-report report.html \
      -with-trace
```

### Miscellaneous

Here is the DAG summarizing the multi-runs SCENIC workflow:

![Multi-Runs Motif and Track based SCENIC](assets/multi_runs_motif_track_scenic.svg)

## Testing the pipeline

```{bash}
nextflow -C conf/test.config,conf/test_multi_runs.config run main.nf --test
```
# nf-core/eager: Configuration for other clusters

It is entirely possible to run this pipeline on other clusters, though you will need to set up your own config file so that the pipeline knows how to work with your cluster.

> If you think that there are other people using the pipeline who would benefit from your configuration (eg. other common cluster setups), please let us know. We can add a new configuration and profile which can used by specifying `-profile <name>` when running the pipeline.

If you are the only person to be running this pipeline, you can create your config file as `~/.nextflow/config` and it will be applied every time you run Nextflow. Alternatively, save the file anywhere and reference it when running the pipeline with `-c path/to/config` (see the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more).

A basic configuration comes with the pipeline, which runs by default (the `standard` config profile - see [`conf/base.config`](../conf/base.config)). This means that you only need to configure the specifics for your system and overwrite any defaults that you want to change.

## Cluster Environment
By default, pipeline uses the `local` Nextflow executor - in other words, all jobs are run in the login session. If you're using a simple server, this may be fine. If you're using a compute cluster, this is bad as all jobs will run on the head node.

To specify your cluster environment, add the following line to your config file:

```nextflow
process.executor = 'YOUR_SYSTEM_TYPE'
```

Many different cluster types are supported by Nextflow. For more information, please see the [Nextflow documentation](https://www.nextflow.io/docs/latest/executor.html).

Note that you may need to specify cluster options, such as a project or queue. To do so, use the `clusterOptions` config option:

```nextflow
process {
  executor = 'SLURM'
  clusterOptions = '-A myproject'
}
```

## Software Requirements
To run the pipeline, several software packages are required. How you satisfy these requirements is essentially up to you and depends on your system. If possible, we _highly_ recommend using either Docker or Singularity.
Please see the [`installation documentation`](../installation.md) for how to run using the below as a one-off. These instructions are about configuring a config file for repeated use.

### Docker
Docker is a great way to run nf-core/eager, as it manages all software installations and allows the pipeline to be run in an identical software environment across a range of systems.

Nextflow has [excellent integration](https://www.nextflow.io/docs/latest/docker.html) with Docker, and beyond installing the two tools, not much else is required - nextflow will automatically fetch the [nfcore/eager](https://hub.docker.com/r/nfcore/eager/) image that we have created and is hosted at dockerhub at run time.

To add docker support to your own config file, add the following:

```nextflow
docker.enabled = true
process.container = "nfcore/eager"
```

Note that the dockerhub organisation name annoyingly can't have a hyphen, so is `nfcore` and not `nf-core`.


### Singularity image
Many HPC environments are not able to run Docker due to security issues.
[Singularity](http://singularity.lbl.gov/) is a tool designed to run on such HPC systems which is very similar to Docker.

To specify singularity usage in your pipeline config file, add the following:

```nextflow
singularity.enabled = true
process.container = "shub://nf-core/eager"
```
If you intend to run the pipeline offline, nextflow will not be able to automatically download the singularity image for you.
Instead, you'll have to do this yourself manually first, transfer the image file and then point to that.

First, pull the image file where you have an internet connection:

```bash
singularity pull --name nf-core-eager.simg shub://nf-core/eager
```

Then transfer this file and point the config file to the image:

```nextflow
singularity.enabled = true
process.container = "/path/to/nf-core-eager.simg"
```

By default nextflow will store a singularity image in the working directory of a job. You can alternatively further specify a 'central' singularity cache to keep all singularity contains for a(ll) user(s). This can be
done by either setting a central environmental variable `NXF_SINGULARITY_CACHEDIR` or specifying the location in a nextflow config file with `singularity.cacheDir`.

### Conda
If you're not able to use Docker or Singularity, you can instead use conda to manage the software requirements.
To use conda in your own config file, add the following:

```nextflow
process.conda = "$baseDir/environment.yml"
```

## Job Resources
#### Automatic resubmission
Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

#### Custom resource requests
Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files in [`conf`](../conf) for examples.

### AWS Batch specific parameters
Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use the `-awsbatch` profile and then specify all of the following parameters.
#### `--awsqueue`
The JobQueue that you intend to use on AWS Batch.
#### `--awsregion`
The AWS region to run your job in. Default is set to `eu-west-1` but can be adjusted to your needs.

Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.

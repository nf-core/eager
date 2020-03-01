# Frequently Asked Questions

## I specified a module and it didn't produce the expected output

Possible options:

1. Check there wasn't a typo in the parameter name. Nextflow _does not_ check for this
2. Check that an upstream module was turned on (if a module requires the output of a previous module, it will not be activated unless it receives the output)

## The pipeline cashes almost immediately with an early pipeline step (`fastqc`, `outputdocumentation` etc.)

If you're running singularity, it could be that nextflow cannot access your singularity image properly - often due to missing bind paths.

See [here](https://nf-co.re/usage/troubleshooting#cannot-find-input-files-when-using-singularity) for more information.
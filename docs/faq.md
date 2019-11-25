# Frequently Asked Questions

## I specified a module and it didn't produce the expected output?

Possible options:

1. Check there wasn't a typo in the parameter name. Nextflow _does not_ check for this
2. Check that an upstream module was turned on (if a module requires the output of a previous module, it will not be activated unless it recieves the output) 
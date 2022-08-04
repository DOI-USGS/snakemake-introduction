# Issue 2: Reading Snakefile Parameters into Your Script

Pipeline steps: No new pipeline step

Concepts learned: Snakefile python module to read in variables from snakefile

## Defining Named Variables in a Snakefile
In the last issue, we wrote a Snakefile rule to fetch a zipped folder containing temperature prediction data:
```
rule get_sb_data:
    output:
        "1_fetch/tmp/pgdl_predictions_04_N45.50-48.00_W92.00-93.00.zip"
    script:
        "1_fetch/sb_get.py"
```

There were no file inputs to this step (if there were, we would have defined them as a Snakefile `input`). However, if you look at the `1_fetch/sb_get.py` script, you will notice that there are two pieces of information that could be considered "inputs" to this pipeline step - the ScienceBase item that we are fetching data from, and the name of the file we want to download. These pieces of information have been hardcoded into the script as variables:
```
if __name__ == "__main__":
    sb_item = "5e5d0bb9e4b01d50924f2b36"
    sb_file = "1_fetch/tmp/pgdl_predictions_04_N45.50-48.00_W92.00-93.00.zip"
    main(sb_item, sb_file)
"
```

It can be helpful to specify these variables in the Snakefile instead. That way, you get a more comprehensive view of the inputs and outputs for your pipeline step. It also allows you to more easily update your pipeline if, for example, you were to switch to a different (but identically structured) data source from ScienceBase.

Let's define these variables in the Snakefile now. In just a bit, we'll see the properties of a Snakefile's rule - `input`, `params`, `output`, etc. - can be accessed from within the Python script called by that rule. We have already defined the string stored in the `sb_file` variable in the output of our current rule. Let's give this string a name so that we can reference it by name in our script:
```
rule get_sb_data:
    output:
        sb_file = "1_fetch/tmp/pgdl_predictions_04_N45.50-48.00_W92.00-93.00.zip"
    script:
        "1_fetch/sb_get.py"
```

We have not defined the string stored in the `sb_item` variable of our Python script anywhere in our Snakefile yet. Since it is not an input or output file for the script called in this step, we will define this variable as a parameter for this rule by using the `params` directive:
```
rule get_sb_data:
    params:
        sb_item = "5e5d0bb9e4b01d50924f2b36"
    output:
        sb_file = "1_fetch/tmp/pgdl_predictions_04_N45.50-48.00_W92.00-93.00.zip"
    script:
        "1_fetch/sb_get.py"
```

## Reading Named Variables from a Snakefile into a Python Script
Now that we have defined all the arguments for the `main()` function of our `1_fetch/sb_get.py` script in our Snakefile, let's replace the hardcoded strings in our script. We can access variables defined in a Snakefile rule using the `snakemake` Python object. The `snakemake` object is a global object that Snakemake provides in any script that is called via the `script` directive. We call the variables by the property of the rule (input, output, params), and since we have given them names, we can reference them by name:
```
if __name__ == "__main__":
    sb_item = snakemake.params['sb_item']
    sb_file = snakemake.output['sb_file']
    main(sb_item, sb_file)
```

Save your Snakefile and Python script with these updates, and go ahead and delete the zip file we downloaded earlier (in `1_fetch/tmp/`). Execute the rule again now that its script uses the `snakemake` object, and see if the file downloads correctly for you again:
```
snakemake --cores 1 1_fetch/tmp/pgdl_predictions_04_N45.50-48.00_W92.00-93.00.zip
```

Hopefully the snakemake pipeline executed without any errors and you were able to re-download the output file for our current pipeline: `1_fetch/tmp/pgdl_predictions_04_N45.50-48.00_W92.00-93.00.zip`. If you had any trouble, look back over your work, and contact your instructor if you can't find the problem on your own.

In the next issue, we'll learn some best practices for unzipping files using snakemake.

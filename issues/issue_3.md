# Issue 3: Unzip data files

Pipeline steps: Unzip the temperature predictions

Concepts learned:
- The `input` directive
- Linking rules together via inputs and outputs
- Multiple output files

## Adding a second rule

The next step in our pipeline is to unzip the file with temperature predictions.
That's the file we downloaded in the first step: `pgdl_predictions_04_N45.50-48.00_W92.00-93.00.zip`.
We already have a script that unzips the file for us: `1_fetch/unzip_file.py`.
Like we did with the `get_sb_data` rule, all we need to do is:
1. Write a Snakemake rule to call the script
2. Generalize the script to use filenames provided by Snakemake

First, let's examine the existing script.
Two strings are hardcoded in the script:
1. The name of the file to unzip
2. The directory where we'd like to unzip files to

Instead of hardcoding them in the script, we should specify those paths in the Snakefile.
That will keep the script general and the Snakefile specific.
Let's add the new Snakemake rule to the Snakefile.
```
rule unzip_sb_data:
```

## Using the `input` directive

We want this rule to run after the zip file is downloaded.
We can ensure that the zip file is downloaded first by naming the zip file as an input to this rule.
Just as outputs to a rule are defined using the `output` directive, inputs are defined using the `input` directive.
Let's add the input to our new rule:
```
rule unzip_sb_data:
    input:
        zip_file_path = "1_fetch/tmp/pgdl_predictions_04_N45.50-48.00_W92.00-93.00.zip"
```

Now, when the rule `unzip_sb_data` is called, Snakemake will look for the file `1_fetch/tmp/pgdl_predictions_04_N45.50-48.00_W92.00-93.00.zip`.
If the file does not exist or if it is out-of-date, then Snakemake will look for a way to create it - a rule that has that file as its output.
In this case, that rule is `get_sb_data`.
So, Snakemake will execute `get_sb_data` first, check that `1_fetch/tmp/pgdl_predictions_04_N45.50-48.00_W92.00-93.00.zip` has been downloaded, and then execute `unzip_sb_data`.

> NOTE: We named the input file `zip_file_path`.
> Naming inputs, outputs, and parameters isn't actually required, but it's a recommended practice because it makes code more readable, especially in Python scripts called by Snakemake.

## Multiple outputs

Snakemake works best when it tracks all files in the pipeline properly.
For this unzipping step, the easiest way to let Snakemake track all the files is to explicitly list the name of each file in the zip archive.
These files are outputs of this step, so we use the `output` directive as follows:
```
rule unzip_sb_data:
    input:
        zip_file_path = "1_fetch/tmp/pgdl_predictions_04_N45.50-48.00_W92.00-93.00.zip"
    output:
        "1_fetch/out/tmp/pgdl_nhdhr_{FC091A8F-FC45-46C0-91F9-18379CF0EAAE}_temperatures.csv",
        "1_fetch/out/tmp/pgdl_nhdhr_69545713_temperatures.csv",
        "1_fetch/out/tmp/pgdl_nhdhr_80006805_temperatures.csv",
        "1_fetch/out/tmp/pgdl_nhdhr_86444267_temperatures.csv",
        "1_fetch/out/tmp/pgdl_nhdhr_86445115_temperatures.csv",
        "1_fetch/out/tmp/pgdl_nhdhr_105954753_temperatures.csv",
        "1_fetch/out/tmp/pgdl_nhdhr_107071276_temperatures.csv",
        "1_fetch/out/tmp/pgdl_nhdhr_107071492_temperatures.csv",
        "1_fetch/out/tmp/pgdl_nhdhr_107072210_temperatures.csv",
        "1_fetch/out/tmp/pgdl_nhdhr_111726865_temperatures.csv",
        "1_fetch/out/tmp/pgdl_nhdhr_120018402_temperatures.csv",
        "1_fetch/out/tmp/pgdl_nhdhr_120018788_temperatures.csv",
        "1_fetch/out/tmp/pgdl_nhdhr_120018790_temperatures.csv",
        "1_fetch/out/tmp/pgdl_nhdhr_120020150_temperatures.csv",
        "1_fetch/out/tmp/pgdl_nhdhr_120020163_temperatures.csv",
        "1_fetch/out/tmp/pgdl_nhdhr_120020166_temperatures.csv",
        "1_fetch/out/tmp/pgdl_nhdhr_120020167_temperatures.csv",
        "1_fetch/out/tmp/pgdl_nhdhr_120020444_temperatures.csv",
        "1_fetch/out/tmp/pgdl_nhdhr_120020465_temperatures.csv",
        "1_fetch/out/tmp/pgdl_nhdhr_120020466_temperatures.csv",
        "1_fetch/out/tmp/pgdl_nhdhr_120020478_temperatures.csv",
        "1_fetch/out/tmp/pgdl_nhdhr_120020480_temperatures.csv",
        "1_fetch/out/tmp/pgdl_nhdhr_120020497_temperatures.csv",
        "1_fetch/out/tmp/pgdl_nhdhr_120020979_temperatures.csv"
```
Specifying multiple outputs is as simple as listing many strings separated by commas.
You can specify multiple inputs or parameters in exactly the same way.
Add the code above to the Snakefile.

Snakemake stores everything specified in the `input`, `output`, and `params` directives as Python list-like objects.
That means that you can treat the object `snakemake.output` as a regular Python list, for example. You are able to call items in the Snakemake directive by index (not just by name, like we did before).
So, in a Python script that is called by the Snakefile, you could use `snakemake.output[-1]` to access the last file listed under the output directive in the Snakefile. You could also pass `snakemake.input` as an argument to a function that expects a list-like object without a problem. Something to keep in mind for the future!

> NOTE: We've included the relative path to each output file.
> Snakemake requires all paths to files to be specified relative to the Snakefile.


## Completing the rule

Remember, there are two strings hardcoded in the script file: the zip file path, and the destination directory.
Our rule now hardcodes the zip file path, but not the destination directory.
So, let's hardcode the destination directory in the rule!
The question is which directive to use.
The destination directory obviously isn't an input file.
It isn't an output, because outputs should (almost) always be files, not directories.
So, let's use the `params` directive.

Go ahead and add the output directory as a parameter named `out_dir`.
We can do this just like in the previous issue.
If you're not sure how, you can go back and review the previous issue to refresh your memory.

We'll also need to specify the script to run: `1_fetch/unzip_file.py`.
As with the previous rule, we can use the `script` directive.
Edit the Snakefile and specify the script to run for this rule.

## Modifying the script

Our Snakemake rule is in good shape, but now we need to generalize the `unzip_file.py` script.
Open the script and replace the hardcoded strings with references to our `input` and `params`, using the `snakemake` Python module.
Again, if you're not sure how to do this, refer back to the previous issue.

## Test the new rule

Once you've got everything ready, you can test your new rule by executing it.
There are two ways to do this.
One way is to call snakemake using the name of any of the unzipped files, like this:
```
snakemake --cores 1 1_fetch/out/tmp/pgdl_nhdhr_120020979_temperatures.csv
```
Snakemake will look for and execute the rule with this file as an output - the `unzip_sb_data` rule.
When that rule is executed, all of that rule's outputs will be created, not just the one you specified.

The other way to execute the rule is to call it by name in the snakemake command:
```
snakemake --cores 1 unzip_sb_data
```
Use either method to test the rule.
Hopefully you'll find all the unzipped files in the directory `1_fetch/out/tmp/`.
If snakemake shows an error message, read it closely and see if you can fix the error.
If you get stuck, ask your instructor for help.

Let's make sure that the pipeline runs from beginning to end correctly.
We'll delete all output files and then create them again.
Try this command:
```
snakemake --cores 1 unzip_sb_data --delete-all-output
```
All the unzipped files and the downloaded archive should be gone.
Now execute the `unzip_sb_data` rule again.
```
snakemake --cores 1 unzip_sb_data
```
Did Snakemake automatically execute the download step first, as we expect it to, and then unzip the archive?

As a one-to-many operation, unzipping is a tricky operation to include in a pipeline.
If you've successfully unzipped the archive using Snakemake, congratulations!
If you have instructor support, then submit a new PR for their review on a branch named "issue-3".

In the next issue we'll add a rule to calculate mean temperatures by day of year for any lake using Snakemake wildcards. When you're all finished here, head over to [Issue 4](issue_4.md).

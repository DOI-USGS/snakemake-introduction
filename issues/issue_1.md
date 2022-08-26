# Issue 1: Introduction to Snakefiles

Pipeline steps: Download the temperature predictions

Concepts learned: 
- Snakefile structure
- The `script` directive
- Basic `snakemake` commands 

## Tutorial Overview

Welcome to the USGS Data Science Practitioner's Snakemake tutorial!
This tutorial consists of several issues.
Each issue will guide you in building a piece of a Snakemake data pipeline.
This pipeline will download, aggregate, and visualize predictions of lake temperatures at multiple depths.
At the same time, each issue will introduce you to new Snakemake concepts.

Most of the Python code for the pipeline is already written.
It's organized into three phases.
- `1_fetch` downloads the lake temperature predictions from this data release on ScienceBase: https://www.sciencebase.gov/catalog/item/5e5d0bb9e4b01d50924f2b36 and unzips the files. This download includes daily temperature predictions at 0.5 meter depth intervals for 24 lakes. The daily predictions span multiple years.
- `2_process` computes the mean predicted temperature on each day of the year (at all depths) for each lake we select and then concatenates the data for each lake into a single .csv file.
- `3_plot` plots the mean temperatures on each day of the year for each lake at a selection of depths.

The code for each of the phases is located in the correspondingly named folder of this repository.
Within each folder are all the Python scripts needed to execute the tasks of that phase. 
If we wanted to carry out these tasks using the scripts as they are, we could execute each script in order.
Alternatively, we can use a Snakemake pipeline.
Snakemake allows us to:

- Describe the order in which scripts should run, and what the expected inputs/outputs of each step are
- Automate the execution of the scripts
- Track which outputs have already been made
- Only run the parts of the pipeline that need to be run
- Rerun parts of the pipeline, or the entire pipeline at once
- Easily expand the pipeline to accommodate more data.

In this issue we'll create a Snakefile and write our first rule.

## What's a Snakefile?

A Snakefile is a file that describes every step a pipeline will take, which files that step depends upon, and which files that step produces.
The Snakefile can also hold or point to configuration settings.
It's like a recipe for how to create outputs, or a blueprint of the pipeline.
A Snakefile is broken down into a sequential set of steps, each of which is called a rule.

## What's a rule?

A rule describes one piece of a pipeline - one step that it can take.
A rule specifies:
- the files that are inputs to that step
- the files that are outputs of that step
- how to execute the step
- parameters and settings that pertain to that step

## Running the script without the pipeline

The rule that we will create in this issue will download the temperature predictions from ScienceBase that we will process and visualize in later steps.
As mentioned earlier, the Python code for this step is already written.
You can find it in the script `sb_get.py`.
Let's test out the Python code before we start to build our pipeline. Navigate to your local clone of this repository in either Anaconda Navigator or the command line. Be sure that you're in the Conda environment we created for this tutorial - if not, activate the right environment first:
```
conda activate snakemake-tutorial
```

Now, try running the script. I'll wait!
```
python 1_fetch/sb_get.py
```

The script should download a single zip file to `1_fetch/tmp/pgdl_predictions_04_N45.50-48.00_W92.00-93.00.zip`.
If the file wasn't downloaded, let your instructor know and we'll get it sorted out.

If the file is there, great!
Now, let's write a Snakemake rule to download that file and build the first step of our pipeline.

## Writing a rule

In the root directory of the repository, create a new empty file called `Snakefile` (yes, with the S capitalized).
Open that file in a text editor, and type
```
rule get_sb_data:
```
That's how you define a Snakemake rule.
> NOTE: It's a lot like the way functions get defined in Python, except instead of `def` we type `rule`.
> In fact, Snakefiles are basically Python files with some added functionality.
> You can include Python code directly in a Snakefile - define variables, or functions, or anything you want - and that code will work the same as if it were in a `.py` file.

Next, we'd usually list the files that are inputs to this rule.
However, since this is the first rule of our pipeline, we don't have any inputs.
So, let's skip straight to outputs.
We'll use the output directive to define the outputs of this rule.
In this case, there's only one output: the zip file.
Easy enough!
Let's add this to our rule:
```
rule get_sb_data:
    output:
        "1_fetch/tmp/pgdl_predictions_04_N45.50-48.00_W92.00-93.00.zip"
```

Finally, we need to tell Snakemake how to create this zip file - what code to run.
There are a few ways to do this, but we'll start with our favorite way: the `script` directive.
The `script` directive tells Snakemake to run a script in order to make the output file(s) for this rule.
In this case, we know that all we need to do is to run `1_fetch/sb_get.py` to download the file.
So, let's specify that Python file as our script.
```
rule get_sb_data:
    output:
        "1_fetch/tmp/pgdl_predictions_04_N45.50-48.00_W92.00-93.00.zip"
    script:
        "1_fetch/sb_get.py"
```
Like with any Python file, be sure to pay attention to indentation!
Since the script has the extension `.py`, Snakemake knows to use `python` to run it.
Unless we specify otherwise, the version of Python that's used to run the script will be the version that's called when you type `python` into the command line.
So, when you want to run the pipeline, be sure that you're in the correct Conda environment.

Okay, we've got our Snakefile written with its first rule!
Let's try it out.

## Executing the pipeline

First, go ahead and delete the zip file we downloaded earlier - the one that's in `1_fetch/tmp/`.
Next, execute our new rule using the following command:
```
snakemake --cores 1 1_fetch/tmp/pgdl_predictions_04_N45.50-48.00_W92.00-93.00.zip
```
This is the main way to run the pipeline.
We call the `snakemake` executable, passing in the file or files we want it to create as arguments.
Notice the `--cores 1` argument.
That specifies how many cores we want to use, and it is a required argument that you must pass each time you run snakemake.
Snakemake is great for running many tasks in parallel.
In this case, though, there's only one task, so parallelization won't help.
So, we're only using one core.

> NOTE: The flags --cores, -c, --jobs, and -j all mean the same thing for local execution.

How did it go?
If there was an error, don't worry.
Look back over the Snakefile and make sure there's no errors in its text, or changes in indentation.
Also, make sure that you're running the command from the root directory of the repository, and that the Snakefile is in the root directory also.
If the `snakemake` command can't be found, then you might not have the Conda environment activated, or you may need to revisit the installation instructions.
If there is a message "Nothing to be done", then make sure you've deleted the zip file and run the `snakemake` command again.
If after all that, there's still a problem, let your instructor know.

If you've been able to download the file using Snakemake, then congratulations!
You've written your first Snakemake rule!

If you're working with an instructor, then create a pull request (PR) for them to review.
Create a new branch in your fork entitled "issue-1", commit your changes to that branch, and push them to your fork.
```
git checkout -b issue-1
git add Snakefile
git commit
git push -u origin issue-1
```
Then, you can submit a PR through GitHub. Make sure that your PR aims to merge into the main branch of YOUR FORK and not the main branch of USGS-R/snakemake-introduction. Request a review from your instructor through GitHub. See the instructions in [this repository's README](../README.md) for more details about the process of submitting PRs.

In the next issue, we'll generalize our workflow so that the ScienceBase item and filename aren't hardcoded into the script, but rather specified in the Snakefile. Go to [Issue 2](issue_2.md) once you're ready.

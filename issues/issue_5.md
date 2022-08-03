# Issue 5: Combining wildcards

Pipeline steps:
- Combine site-specific mean temperature files
- Plot mean temperatures by day of year
- Update rule all 
- Add in third lake

Concepts learned:
- The split-apply-combine strategy
- Rebuilding parts of the pipeline

## The split-apply-combine strategy

The split-apply-combine strategy is a powerful and widespread data analysis approach in which a dataset is split apart into manageable pieces, analysis is performed on each piece, and the results of the analyses are re-combined.
In Snakemake pipelines, this strategy makes it easy to execute the analysis step for each of the pieces in parallel.
That parallelization can save hours or days if the analyses are computationally costly.

Let's consider our pipeline thus far in the context of the split-apply-combine strategy.
We have three rules so far that accomplish the following three tasks:

1. Download a dataset as a zip file (all lakes)
2. Unzip the file (obtain lake-specific files)
3. Calculate mean temperatures by day-of-year (one lake at a time)

So far, we've SPLIT our dataset apart by lake (actually, it came that way - all we had to do was unzip a file).
We've APPLIED the computation of means by day-of-year to each lake we're interested in.
Next, we'll COMBINE those means back together to make it easier to compare and visualize temperatures in different lakes at different depths.

## Combining the site files

Our next rule will combine the processed site files into a single file.
Let's call the rule `combine_site_files`.
Once again, most of the code is already written.
This time, the relevant code can be found in `2_process/combine_site_files.py`.
Our job is to Snakemake-ify it!

Let's start by writing the rule.
The inputs to this combining step should be the day-of-year averaged .csv files that are outputs of `calc_doy_means`.
Notice that these are the same files we used earlier to create rule `all`.
We can copy that code and paste it as the inputs for our new rule.
```
rule combine_site_files:
    input:
        "2_process/out/doy_120020150.csv",
        "2_process/out/doy_107072210.csv"
```
The wildcard in rule `calc_doy_means` will be resolved in the same way that they are when rule `all` requests these files as inputs.

There is only one output file for this step: a single .csv file with all the averaged temperatures for every lake.
In `2_process/combine_site_files.py` you can see that the file is saved to `2_process/out/combined_doy.csv`.
So, we'll include that and the path to the Python script to complete our new rule.
```
rule combine_site_files:
    input:
        "2_process/out/doy_120020150.csv",
        "2_process/out/doy_107072210.csv"
    output:
        out_file = "2_process/out/combined_doy.csv"
    script:
        "2_process/combine_site_files.py"
```

Now it's time to modify the Python code.
Open `2_process/combine_site_files.py` and change the code to use the rule's `input` and `output` values in place of hardcoded file paths, using the `snakemake` module like before.
Refer to the previous issues for a refresher on how to do this.

Once done, use the `snakemake` command to test the new rule and create the file `2_process/out/combined_doy.csv`.
Check that the file has been created, and that it contains the contents of both `2_process/out/doy_107072210.csv` and `2_process/out/doy_120020150.csv`.
If everything looks good, let's move forward!
Reach out to your instructor if you're stuck.

## Plotting the mean temperatures

One more rule to add!
This final rule will plot the average temperatures at multiple depths for each of the sites we specified as inputs to `combine_site_files`.
We'll name it `ploy_doy_mean`.
For this rule we'll need one input file and one output file.
The input file is the output of the previous step: `2_process/out/combined_doy.csv`.
The output file is the plot we'll make as a `.png` file: `3_plot/out/doy_plot.png`.
The Python script we'll adapt to use for this step is `3_plot/plot_doy_mean.py`.

If you read the Python script, you'll see that there's one other hardcoded bit of information: the lake depths at which to plot temperatures.
Let's include that in our rule as well using the `params` directive.

Try creating this rule yourself.
Then, edit the Python script as needed to use the values you defined in the Snakefile instead of hardcoded values.
Run `snakemake` to create the plot and test your new addition to the pipeline.
Does the plot look like it contains the information you expect it to?


## Updating rule `all`

Now that we've added to the end of the pipeline, our rule `all` is out of date.
Let's update it to create the final output of the pipeline: the plot image.
Change the input to rule all to be that single file.
Now, let's clear all our outputs and test rule `all`.
Remember how to delete all output files?
```
snakemake --cores 1 3_plot/out/doy_plot.png --delete-all-output
```
Next, we'll try a dry run without specifying any output files.
First, what do you expect to see?
Now carry out the dry run and see if the output matches your expectations.
If it all looks good, then execute the actual run.
```
snakemake --cores 1
```
Did the pipeline run as expected?
If not, track down any problems, or ask your instructor for assistance.
If so, congratulations! The pipeline is complete.

## Changing and rebuilding pipeline targets

Now, with a complete pipeline, we can play around a bit!
Let's run some experiments to test Snakemake's limits and better understand when Snakemake does or does not rebuild a part of a pipeline.
First things first: let's delete all output files and rebuild the pipeline from scratch.
```
snakemake -c1 3_plot/out/doy_plot.png --delete-all-output
snakemake -c1 3_plot/out/doy_plot.png
```
> NOTE: -c1 is short for --cores 1

Now let's see what happens when we add another lake to the plot.
From among all the site IDs of the csvs that are outputs of `unzip_sb_data`, choose one site and add it to the inputs of `combine_site_files`.
For instance, you could add "2_process/out/doy_86444267.csv" as a third output of `combine_site_files`.
Now try a dry run to see what happens when we build the pipeline again.
```
snakemake -n
```
> NOTE: -n is short for --dry-run

You'll probably see the following message:
> Building DAG of jobs...
> The input used to generate one or several output files has changed:
>     To inspect which output files have changes, run 'snakemake --list-input-changes'.
>     To trigger a re-run, use 'snakemake -R $(snakemake --list-input-changes)'.
> Nothing to be done (all requested files are present and up to date).

We know that Snakemake won't rebuild any pipeline targets because it prints the message, "Nothing to be done".
This is important to notice - when the Snakefile is changed, Snakemake doesn't rebuild everything affected by that change by default.
However, Snakemake did detect the change to inputs that we made and offers a way to re-run the pipeline with the changes.
First, let's inspect which output files have changes, as it suggests.
```
snakemake --list-input-changes
```
One file appears: `2_process/out/combined_doy.csv`. 
That makes sense - the inputs for that file changed.
Let's try re-running the pipeline in the way the Snakemake warning message suggests.
First, let's do a dry run.
```
snakemake -n -R $(snakemake --list-input-changes)
```
> NOTE: The `-R` flag is short for --forcerun. The flag forces a file to be re-created, even if Snakemake would otherwise not remake it. All dependent (downstream) files are also re-created, but the files that the target depends on (upstream) are not remade.

You can see that, if we execute this command without the `-n` flag, we'll create the averaged temperatures for the new lake in the file `2_process/out/doy_86444267.csv` and re-create `2_process/out/combined_doy.csv` and `3_plot/out/doy_plot.png`.
This entails re-running four rules: `calc_doy_means`, `combine_site_files`, `plot_doy_mean`, and `all` (rule `all` doesn't actually create any files).
Notice that only the new lake's averaged file gets created; `doy_120020150.csv` and `doy_107072210.csv` remain untouched because adding a new lake has no effect on those target files.
However, `2_process/out/combined_doy.csv` concatenates averaged temperatures for all specified lakes, so it gets re-created along with `3_plot/out/doy_plot.png`, which depends on it.

Let's try re-running the pipeline using the command Snakemake suggested.
We'll need to make one change to the command as written: we need to include the `-c` flag to tell Snakemake how many cores to use, or it won't run.
```
snakemake -c1 -R $(snakemake --list-input-changes)
```
Now open `3_plot/out/doy_plot.png`.
There's our newly added lake, plotted with the others!

You can similarly change the depths shown in the plots by changing the value of `depths` in the params of the `plot_doy_mean` rule.
Try that out and see if Snakemake rebuilds the file by default after your change to the params, or if it only warns you of a change.
Feel free to continue making changes and exploring under what conditions does Snakemake recreates files or provides warning messages.


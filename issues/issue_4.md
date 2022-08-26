# Issue 4: Generalizing Rules with Wildcards

Pipeline steps: Compute the mean predicted temperature on each day of the year (at all depths) for each lake of interest

Concepts learned:
- wildcards
- rule all
- `rule all`
- dry-run
- the `dry-run` argument

## Update `calc_doy_means.py` script to use properties from Snakemake rule
Our current workflow is able to fetch and unzip raw lake temperature prediction data from ScienceBase. In the `out` file of our `1_fetch` step, we have a folder called `tmp` containing a set of csv files, each of which contains daily temperature predictions for a single lake (at multiple depths).

The daily dataset spans nearly 40 years, and we would like to process it by calculating the day-of-year mean temperature predictions (at all depths) for each lake. Let's add a new rule to do this.

For this tutorial, we will calculate the day-of-year mean temperatures for two lakes. Those lakes have the IDs `120020150` and `107072210`.

Let's process the data based on our current knowledge of Snakemake pipelining. We can write a rule that uses the `2_process/calc_doy_means.py` script. We will need to indicate the input and output files for the two lakes we have chosen in our rule as well:
```
rule calc_doy_means:
    input:
        "1_fetch/out/tmp/pgdl_nhdhr_120020150_temperatures.csv",
        "1_fetch/out/tmp/pgdl_nhdhr_107072210_temperatures.csv"
    output:
        "2_process/out/doy_120020150.csv",
        "2_process/out/doy_107072210.csv"
    script:
        "2_process/calc_doy_means.py"
```

This rule is perfectly acceptable and should work if you test it out by running `snakemake --cores 1 calc_doy_means` in your command line/Anaconda prompt.

> Note: If you do test out the rule, be sure to delete the outputs before moving forward (`snakemake --cores 1 calc_doy_means --delete-all-output`). If you don't delete the outputs, the Snakemake pipeline will not execute next time because the outputs have already been generated!

Next, let's replace the hardcoded file names and lake IDs in the `2_process/calc_doy_means.py` script with our Snakemake rule properties.

In our Snakefile, we have multiple files as inputs (and multiple files as outputs too), none of which are given variable names in the Snakefile. When we call `snakemake.input` in our Python script, a list of the inputs from the Snakefile will be returned (`["1_fetch/out/tmp/pgdl_nhdhr_120020150_temperatures.csv", "1_fetch/out/tmp/pgdl_nhdhr_107072210_temperatures.csv"]`). Likewise, a list of output file names will be returned when we call `snakemake.output`.

Since the Snakemake rule inputs and outputs already have the lake IDs included, we can remove the list of `lake_ids` from our script. However, since we were using the lake IDs to indicate how many iterations our for loop needs to go through, we will now need to change that. We can choose the number of iterations by using the length of our input file list: `for i in range(len(snakemake.input))`. We can call each input and output file name from the list by index as we go through our for loop.

Let's put that all together:
```
if __name__ == '__main__':
    for i in range(len(snakemake.input))
        out_file = snakemake.output[i]
        in_file = snakemake.input[i]
```

Lastly, let's extract the `lake_id` from our `out_file` name in each iteration of the for loop. Here's the full code block to update your `calc_doy_means.py` script with:
```
if __name__ == '__main__':
    for i in range(len(snakemake.input)):
        out_file = snakemake.output[i]
        in_file = snakemake.input[i]
        lake_id = os.path.splitext(out_file)[0].split('_')[-1]
        main(out_file, in_file, lake_id)
```

Try running this step again with `snakemake --cores 1 calc_doy_means` and make sure it works for you. If it does, delete the outputs before we move on with `snakemake --cores 1 calc_doy_means --delete-all-output`.

## Generalize your rule with wildcards
As you might imagine, if we wanted to process data for a large number of lakes, listing every input and every output could become very tedious and result in a very long Snakefile. Because the input filepaths are identical aside from just one varying substring (the lake ID), we can simplify our Snakefile down to just one line to represent all the input files. The same is true for the output files. This generalization will be done using "wildcards".

Wildcards are indicated by using `{}` around the wildcard name (`lake_id` in this case). Here's our updated rule that uses wildcards:
```
rule calc_doy_means:
    input:
        in_file = "1_fetch/out/tmp/pgdl_nhdhr_{lake_id}_temperatures.csv"
    output:
        out_file = "2_process/out/doy_{lake_id}.csv"
    script:
        "2_process/calc_doy_means.py"
```
Don't forget to name the variables in the rule (`in_file` and `out_file`)! They aren't required, but it's good practice to specify names for inputs, outputs, and params. It will make the code in your script easier to read and understand. Of course, if there's a long list of input or output files it may not make sense to assign a variable name to every file in the list, but it's good to name variables whenever practical.

## Resolve wildcards and introduction of `rule all`
We've set up our Snakefile with some wildcards, but we've lost the information about which `lake_ids` specifically we want to process. We would typically resolve our wildcards in the next pipeline rule as inputs, like this:
```
rule whatever_we_want_to_do_next:
    input:
        "2_process/out/doy_120020150.csv",
        "2_process/out/doy_107072210.csv"
```

This would allow Snakemake to intuit that the `lake_id` wildcard should have values `120020150` and `107072210`.

We aren't ready to build the next step of our pipeline yet, but we can still create another rule to resolve our wildcards. We could create a new rule with any name just like above and run our pipeline by calling that named rule (in the example above, that would be `snakemake --cores 1 whatever_we_want_to_do_next`).

However, it is generally considered a good Snakemake practice to collect your final outputs in a `rule all`. The naming of this rule does not make any difference, but this is a standard convention that is recommended by the source documentation.

The `rule all` is typically defined at the top of your Snakefile. The reason for placing it at the top of your Snakefile is that if no target rule is specified in the command line, Snakemake will use the first rule of the Snakefile as the target. This means that we could just run `snakemake --cores 1`, and Snakemake would use the first rule as your target. 

Let's add a `rule all` to the top of our Snakefile that will resolve our `lake_id` wildcard:
```
rule all:
    input:
        "2_process/out/doy_120020150.csv",
        "2_process/out/doy_107072210.csv"
```

## Test your Snakefile structure with a dry-run
Let's check to make sure we have built our wildcards correctly. So far, we've been testing our pipeline by running it and actually generating the outputs as we go. This can be a tedious way to test the pipeline because some steps may be slow to run.

We may want to test out our Snakefile to see if we have resolved all of our wildcards correctly without actually running the pipeline. We can do this by doing a `dry-run`. The `dry-run` command line argument will allow you to see how the wildcard value is propagated to your input and output filenames in your command line. Try it out!
```
snakemake --dry-run 
```

## Update `calc_doy_means.py` script to use the new rule properties
Now that we know our wildcards can be correctly resolved from our Snakefile, we just have one more update to make.

Our `calc_doy_means.py` is currently set up to loop through our list of input files within the script. However, we have replaced our list of input files with wildcards. We want to update our `calc_doy_means.py` to expect a single input file - Snakemake will now call the script many times, once for each `lake_id`. The input and output properties of our rule will be read into our python script with any wildcards already resolved, so we don't need to worry about formatting them with the `lake_id` in the script. We can also make use of the wildcards property of our rule to determine the `lake_id` instead of pulling it out of a filename. Here's the code we want to update our `calc_doy_means.py` script with:
```
if __name__ == '__main__':
    out_file = snakemake.output['out_file']
    in_file = snakemake.input['in_file']
    lake_id = snakemake.wildcards['lake_id']
    main(out_file, in_file, lake_id)
```

We're all done - time to run our pipeline! Now that we have used wildcards in our Snakefile, that rule can actually be parallelized. Let's tell snakemake to use two cores so that each lake's data can be processed in parallel. Try it out:
```
snakemake --cores 2
```

If you're working with an instructor, submit a new PR on a branch named "issue-4".

We've successfully learned how to implement wildcards in our pipeline and introduced the concepts of the `rule all` and how to test your pipeline with a `dry-run`.

In the final issue, we'll learn how to update your Snakefile with new wildcards, as well as how to combine these wildcards back together to produce a single file and plot for all of the lakes. Head to [Issue 5](issue_5.md) once you're ready.

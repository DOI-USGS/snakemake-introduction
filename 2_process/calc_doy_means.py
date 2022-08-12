import pandas as pd
import os

def read_pred_csv(csv_path):
    df = pd.read_csv(csv_path,
                     parse_dates = ['date'],
                     infer_datetime_format = True)
    return df

def calc_doy_means(df, site_id, out_file=None):
    doy_means = df.groupby(df.date.dt.day_of_year).mean()
    doy_means.index.name = "doy"
    doy_means["site_id"] = site_id
    if out_file:
        doy_means.to_csv(out_file)
    return doy_means

def main(out_file, in_file, lake_id):
    out_dir = os.path.dirname(out_file)
    if not os.path.exists(out_dir):
	       os.makedirs(out_dir)
    df = read_pred_csv(in_file)
    calc_doy_means(df, lake_id, out_file)

if __name__ == '__main__':
    out_file = str(snakemake.output)
    in_file = str(snakemake.input)
    lake_id = str(snakemake.wildcards)
    main(out_file, in_file, lake_id)

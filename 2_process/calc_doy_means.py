import pandas as pd
import os

def read_pred_csv(csv_path):
    df = pd.read_csv(csv_path,
                     parse_dates = ['date'],
                     infer_datetime_format = True)
    return df

def calc_doy_means(df, site_id, out_file=None):
    doy_means = df.groupby([df.date.dt.month, df.date.dt.day]).mean()
    doy_means = doy_means.reset_index(drop=True)
    doy_means.index.name = "doy"
    doy_means["site_id"] = site_id
    if out_file:
        doy_means.to_csv(out_file)
    return doy_means

def main(out_dir, in_dir, lake_ids):
    if not os.path.exists(out_dir):
	       os.makedirs(out_dir)
    for i in lake_ids:
        pred_csv_file = os.path.join(in_dir,
                                 f"pgdl_nhdhr_{i}_temperatures.csv")
        df = read_pred_csv(pred_csv_file)
        calc_doy_means(df, i, f"2_process/out/doy_{i}.csv")

if __name__ == '__main__':
    out_dir = "2_process/out"
    in_dir = "1_fetch/out/tmp/"
    lake_ids = ["120020150", "107072210"]
    main(out_dir, in_dir, lake_ids)

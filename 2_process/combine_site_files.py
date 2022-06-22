import pandas as pd
import os

def combine_site_files(site_files, out_file=None):
    df_list = [pd.read_csv(f) for f in site_files]
    combined_df = pd.concat(df_list)
    if out_file:
        combined_df.to_csv(out_file)
    return combined_df

def main(out_dir, lake_ids):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    site_files = [os.path.join(out_dir, f"doy_{i}.csv") for i in lake_ids]
    out_file = os.path.join(out_dir, "combined_doy.csv")
    combine_site_files(site_files, out_file)

if __name__ == '__main__':
    out_dir = "2_process/out"
    lake_ids = ["120020150", "107072210"]
    main(out_dir, lake_ids)

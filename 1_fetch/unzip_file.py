import zipfile
import os

def unzip_file(zip_file_path, out_dir):
    with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
        zip_ref.extractall(out_dir)

def main(zip_file_path, out_dir):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    unzip_file(zip_file_path, out_dir)

if __name__ == '__main__':
    zip_file_path = "1_fetch/tmp/pgdl_predictions_04_N45.50-48.00_W92.00-93.00.zip"
    out_dir = "1_fetch/out/"
    main(zip_file_path, out_dir)

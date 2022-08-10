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
    zip_file_path = snakemake.input['zip_file_path']
    out_dir = snakemake.params['out_dir']
    main(zip_file_path, out_dir)

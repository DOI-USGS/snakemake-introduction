import os
from sciencebasepy import SbSession

def sb_get(item_id, sb_data_file=None, destination_dir='.'):
    """
    Download an item from ScienceBase by ID.
    Written by @AndyMcAliley

    :param item_id: ScienceBase ID of item to download
    :param sb_data_file: Name of file to download. If None, download all files. (Default value = None)
    :param destination_dir: directory to save to (Default value = '.')
    :returns: ScienceBase JSON response

    """
    sb_session = SbSession()
    if not (sb_session.ping()['result'] == 'OK'):
        raise ConnectionError('ScienceBase ping unsuccessful')
    # If data are public, no need to log in
    # if not sb_session.is_logged_in():
        # sb_session.login()

    item_json = sb_session.get_item(item_id)
    if sb_data_file is None:
        response = sb_session.get_item_files(item_json, destination=destination_dir)
    else:
        all_file_info = sb_session.get_item_file_info(item_json)
        file_found = False
        for file_info in all_file_info:
            if file_info['name'] == sb_data_file:
                file_found = True
                response = sb_session.download_file(file_info['url'], file_info['name'], destination_dir)
        if not file_found:
            raise FileNotFoundError(f'{sb_data_file} not found on ScienceBase')

    return response

def main(sb_item, sb_file):
    out_dir = os.path.dirname(sb_file)
    sb_data_file = os.path.basename(sb_file)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    sb_get(sb_item, sb_data_file, out_dir)

if __name__ == "__main__":
    sb_item = "5e5d0bb9e4b01d50924f2b36"
    sb_file = "1_fetch/out/pgdl_predictions_04_N45.50-48.00_W92.00-93.00.zip"
    main(sb_item, sb_file, out_dir)

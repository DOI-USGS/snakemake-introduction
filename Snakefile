rule get_sb_data:
    params:
        sb_item = '5e5d0bb9e4b01d50924f2b36'
    output:
        sb_file = "1_fetch/tmp/pgdl_predictions_04_N45.50-48.00_W92.00-93.00.zip"
    script:
        "1_fetch/sb_get.py"
    



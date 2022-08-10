rule get_sb_data:
    output:
        "1_fetch/tmp/pgdl_predictions_04_N45.50-48.00_W92.00-93.00.zip"
    script:
        "1_fetch/sb_get.py"



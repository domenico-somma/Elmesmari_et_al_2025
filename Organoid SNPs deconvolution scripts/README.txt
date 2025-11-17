# Souporcell Adaptation for BD Rhapsody Data
# The standard souporcell pipeline expects barcodes in a different format.
# Specifically, BD Rhapsody barcodes are numerical and need to be converted
# to their corresponding nucleotide sequences for souporcell compatibility.
#
# Steps:
# 1. Convert the numerical whitelist barcodes to nucleotide sequences.
python3 convert_BD_barcode_enh.py
#
# 2. Use the converted barcodes to extract the correct read sequences from the original BAM file.
nohup ./renamer_BD.sh &
#
# 3. Run the main souporcell pipeline using the newly created, compatible FASTQ files.
nohup ./souporcell_BD_pipeline.sh &
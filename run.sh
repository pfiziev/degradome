echo "running 010determine_read_type.py ${@: -1}"
python 010determine_read_type.py ${@: -1}

echo "running 015map_hits_to_regions.py ${@: -1}"
python 015map_hits_to_regions.py ${@: -1}

echo "running 020hit_statistics.py ${@: -1}"
python 020hit_statistics.py ${@: -1}

echo "running 050check_peak_regions_for_known_snoRNAs.py ${@: -1}"
python 050check_peak_regions_for_known_snoRNAs.py ${@: -1}

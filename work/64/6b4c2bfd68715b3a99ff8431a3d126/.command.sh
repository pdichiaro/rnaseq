#!/usr/bin/env bash -C -e -u -o pipefail
echo "Files received:" > collected_files.txt
for f in file1_process2.txt file2_process2.txt file1_process3.txt file1_process1.txt file2_process1.txt; do
    echo "  - $f" >> collected_files.txt
done

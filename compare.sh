# Which test to run from 0..2
test="$1"

# Test sequence pairs
file1[0]="data/LITMUS28i.gbk.txt"
file2[0]="data/LITMUS28i-mal.gbk.txt"
file1[1]="data/M13KE.gbk.txt"
file2[1]="data/M13KO7.gbk.txt"
file1[2]="data/Adenovirus-C.txt"
file2[2]="data/Enterobacteria-phage-T7.txt"

# Make executables
cd smwt-org
make all
cd ../smwt-paral 
make all
cd ../

# Run executables 
./smwt-org/smwt 0.33 1.33 ${file1[test]} ${file2[test]} > org.txt
./smwt-paral/smwt-paral 0.33 1.33 ${file1[test]} ${file2[test]} > paral.txt

# Compare outputs
diff org.txt paral.txt

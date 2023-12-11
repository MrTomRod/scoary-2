# Running Scoary and Scoary2 on the same data

Dataset: 100 randomly picked and binarized traits from the Scoary2 dataset.

**1) Run Scoary**

```bash
echo "s1 start: $(date +"%T")" >> runtime.txt

podman run --user 0:0 --rm -it -v ./:/data:Z biocontainers/scoary:v1.6.16-1-deb_cv1 \
scoary -t 100_traits.csv -g N0_count.csv -s 2 -o s1_out --permute 1000 --correction I -p 0.1

echo "s1 end: $(date +"%T")" >> runtime.txt
```

**2) Run Scoary2**

```bash
echo "s2 start: $(date +"%T")" >> runtime.txt

podman run --rm -v ./:/data:Z troder/scoary-2 \
scoary2 \
--genes N0_count.csv \
--gene-data-type 'gene-count:,' \
--traits 100_traits.csv \
--trait-data-type 'binary:,' \
--multiple_testing native:0.1 --n-permut 1000 \
--n-cpus 8 \
--random-state 42 \
--outdir s2_out \
--trait_wise_correction

echo "s2 end: $(date +"%T")" >> runtime.txt
```

## Results

```bash
$ cat runtime.txt
s2 start: 14:58:12
s2 end: 14:58:35
s1 start: 14:58:58
s1 end: 15:21:31
```

- Scoary2 took 23 seconds
- Scoary took 22 minutes and 33 seconds or 1353 seconds
- Scoary2 is 1353 / 23 = **59 times** faster than Scoary on this dataset

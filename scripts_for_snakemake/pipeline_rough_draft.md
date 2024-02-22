### step 1: run scansnv
https://github.com/parklab/scan-snv
example script ./run_scan_snv.sh
### step 2: run bicseq2-norm
https://www.math.pku.edu.cn/teachers/xirb/downloads/software/BICseq2/BICseq2.html
example script ./run_bicseq2_norm.sh
### step 3: run scanner
```python
import hiscanner
# define your json file path
path = "your/dir/to/json"
# preprocess
hiscanner.pp.preprocess(path)
# segment
hiscanner.tl.segment(path)
# infer copy number
hiscanner.tl.infer_copy_number(path)
# visualize 
hiscanner.pl.plot_whole_genome_track(path)
```
## Questions
1. snakemake (scansnv) inside snakemake?
2. need intermediat steps (create config file for bicseq2 and scanner)
3. how to best monitor memory/time usage
4. how to best monitor progress and errors
5. environment management
# Channel Model with Memory for DNA Data Storage with Nanopore Sequencing

Simulator to simulate the entire DNA Data Storage process (from synthesis to basecalling)

The Channel Model with Memory for DNA Data Storage with Nanopore Sequencing is published in ISTC2021 conference. If you use this Simulator in your research, please cite this manuscript:
<br>
> B.  Hamoum,  E.  Dupraz,  L.  Conde-Canencia,  and  D.  Lavenier, “Channel Model with Memory for DNA Data Storage with Nanopore Sequencing,” in ISTC 2021 - 11th International Symposium on Topics in Coding.  Montreal, Canada: IEEE, Aug. 2021, pp. 1–5. (https://hal.archives-ouvertes.fr/hal-03337117/).
<br>

This is the first version of the Channel Model with Memory for DNA Data Storage with Nanopore Sequencing based on Julia (https://julialang.org/downloads/).
Julia was used with 1.5.3 version and recent ones should also be compatible (tested until v1.6.2).

After setting Julia correctly, simulations can be launched using next command:
```sh
> ./DNA_data_storage_channel.sh -i example/ref.txt -n 100 -o example/sim.fastq -k 6
```	
**Parameters**:

 * **-i**:  Path to the input sequence (should be on fasta format) to simulate. An example of such sequence is available on "example" folder with the name 'ref.txt' and length=1086. [required]
 * **-n**:  Number of reads to simulate. [required]
 * **-o**:  Simulated sequences output path. Will be presented in a fastq format (without included scores). [required]
 * **-k**:  Channel memory length.  fixed to k=6 by default [recommended]

**Error Profile**:

Error profiles were computed using DNA data storage data which went through multiple steps including:
* **_Synthesis_**: chemical technique, with oligonucleotide assembling  called "GeneArt Strings DNA Fragments" and done by *Thermofisher*.
* **_Sequencing_**: MinIon sequencer with "R9.4_180mv_450bps" nanopore model  
* **_Basecalling_**: Guppy v5.0.7 using "super-accurate" mode. 


More documentation will be added as soon as possible. Meanwhile feel free to reach me by e-mail (**belaid.hamoum@gmail.com**) for more details.

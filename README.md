**A NOVEL TOOL FOR THE DETECTION OF COPY NUMBER VARIATIONS**
The scope of the thesis span on the basic techniques required to develop a tool for the detection of copy number variations.Using a selection of Illumina sequencing datasets, this manual walks users through the fundamentals of reproducing result as generated in the thesis. The approach incorporate a set of R, Python packages and command-line toolbox for finding copy number variations throughout the genome using high-throughput sequencing.

**Abbreviations**

-HTS: High-Throughput Sequencing

-NGS: Next Generation Sequencing

-BAM: Binary Alignment Map

-BED: Browser Extensible Data

-CNVs: Copy Number Variations


**Python dependencies**
As a first step, if user haven't already satisfied these dependencies on your system, install these Python packages via pip or conda:

    Biopython
    matplotlib
    NumPy
    SciPy
    Pandas
    pyfaidx
    pysam

**Step by step Process**
In order for the method to function properly, it first converts files that are in the BAM format into read count matrices or genomic ranges objects. These are the objects that comprise the input. A file in the BED format is used to construct the aforementioned genomic ranges.


# mRNA display sequencing analysis

A set of scripts written to analyze illumina data for an mRNA display selection. 
Scripts: 
(1) process the data from the initial FASTQ file format 
(2) analyze the peptide sequences obtained and can be used for motif analysis, copy number counts, sequence clustering and family generation (among other things)  

As a note: I do not advocate for keeping only sequences that match exactly the theoretical library length, but instead to trim by the known N- and C-terminal regions. Things much shorter can be and often are real...the majority of my Ph.D was spent studying a sequence that was 29 AAs in length, not the 30 of the library.   

Also, this code is quite old and specific to the library i used, would only use it as a guide / help suggest analyses to run. 

# _kite_: kallisto index tag extractor

This package offers a few utilities to help with demultiplexing samples using kallisto bus. 



kallisto bus, introduced in v0.45, allows us to quickly demultiplex single cell samples by just having a list of acceptable tags, regardless of where the tag appears.

The notebook [kite_citeseq_SRR8281307](https://github.com/pachterlab/kite/blob/master/docs/kite_citeseq_SRR8281307.ipynb) in the `docs` folder implements the complete kite pipeline. Note that you must have [kallisto](https://pachterlab.github.io/kallisto/download) and [BUStools](https://github.com/BUStools/bustools/releases) installed to run the entire notebook (or if you can't install just call call the executable file for kallisto and BUStools)

A poster giving an overview of the ideas behind kite is availabke at https://drive.google.com/file/d/1TKyjlFvHwwxfF7UqPn_-LMl7TBavp9Xm/view?usp=drivesdk

## _kite_ Utilities

#### `make_mismatch_map(filename)`

This function returns all sample tags and and their single base mismatches (hamming distance 1). 

It returns a dictionary 4 * tag_length * number_of_tags the size of the original tag.

e.g. tags AA and TT  result in 4*2*2 = 16 tags.


#### `save_mismatch_map(tagmap, tagmap_file_path)`

Saves the generated tagmap to file in fasta format for building the kallisto index.



## Retrieving tags for sample demultiplexing with kallisto bus: Cell hashing as an example

The notebook [kite_citeseq_SRR8281307](https://github.com/pachterlab/kite/blob/master/docs/kite_citeseq_SRR8281307.ipynb) in the `docs` folder implements the complete kite pipeline

The examplation below uses as an example the tag structure from the article **Publication: Cell hashing enable sample multiplexing, multiplet identification and super-loading on droplet-based single cell RNA-sequencing platforms**, puublished on Genome Biology: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1603-1

kallisto bus, introduced in v0.45, allows us to quickly demultiplex single cell samples by just having a list of acceptable tags, regardless of where the tag appears.

To do so, instead of looking for transcript fragments of length `k` (the `k-mer` length) in the reads, which is what kallisto usually does, we look for tag fragments of lenght `k`. Thus, instead of building a transcriptome index and using kallisto bus as normal, we build an index from the list of tags and use kallisto bus as usual. The result will be a BUS file: a list of the cell barcode, UMI, tag equivalence class, and multiplicity (how many times the same entry was found). Because tags are designed to ber far appart from each other in the sequence space, each k-mer should map uniquely to one tag, thus the tag equivalence class should have a one-to-one correspondence to a hashtag (A-H) 

In the cell hashing study, the hashtag sequences have length 12, followed by `C,T,G`. Because kallisto needs to use an odd value for `k`, we have two options: 
We could use a smaller k-mer (11 or 9) or we could extend the barcode sequence by appending one extra base in our list of expected barcode, using `k = 13`, which is what we do here. If the sequence flanking the barcodes is constant, we just include them and that's it. In the cell hashing tag structure they are folowed by a `B`, which is `C,T,G` but not `A`, so for each hashtag we need to make 3 variants. So for example, for hashtag A we use the following sequences:

```
hashtagA-C: AGGACCATCCAAC
hashtagA-T: AGGACCATCCAAT
hashtagA-G: AGGACCATCCAAG
```

Thus the 8 original tags become 24. We could just use those 24 extended hashtags to build our index and things would work fine, but kallisto does only perform exact matches. So a read that has a single letter mismatch in the tag would be thrown away, and give that's about 1-2% of the reads, we can do a little better. The solution is to use a list of all the correct tags plus each tag sequence with one mismatch (hamming distance 1). Since each tag is length 12, and each site can mutate to the 3 other bases, that becomes `3*12 = 36` new sequences, plus the correct sequence, so 37. 

So because we're extending each tag sequence in 3 different ways, in total we're building the index with `8*37*3 = 888` sequences

Because the mapping will be one-to-one, collapsing them to retrieve the 8 original sequences is very straightforward with pandas. 

After doing the UMI collapsing, we will have a count matrix of cell barcodes by hashtag, which we can use to decide which hashtag to ascribe to each barcode, and then demultiplex the data.

# PRguideline
Marjan Farahbod  
February 2, 2016  
## Where to look for the paper

There are several journals you can check, below is a list of some journals that I have come across, based on my field of research. There are plenty more.  
 
* Plos Computational Biology 
* Bioinformatics
* Nature Biotechnology
* Nature Methods
* BMC Biomedical 
* Genome Research
* Molecular Systems Biology
* Genome Biology

You can also look at conference proceedings. Here are some relevant conferences:

* ISMB
* Recomb
* Pacific Symposium on Biocomputing 

It is a good idea to check for the impact factor of the journal, that would give you a sense of journal's quality in your field. 

## How to pick it

Well, you probably already have some interests in the field. A particular disease or tissue, for example. Look for that in the title of the papers. Papers with some sort of **expression data** (i.e: RNA-seq, microarray) are a good choice, as you can compare their analytical methods with ones you learn in the course. Some journals are specific to a field, like "Molecular Psychiatry" or "The journal of Allergy and Clinical Immunology". You can use Google Scholar to find papers about the subject of your interest as well. 

Paul and I browsed the some recent issues of Nature and Science and were able to spot some articles (quite easily) which have some sort of _expression data_ analysis and therefore are related to the course material. 

#### Don't panic if the Biology is new
Different subjects in Biology might have many new terms and phrases for someone who is not a Biologist or Bioinformatician working in that particular field. Look up the words and try to get to the main purpose of the paper. You don't need to understand every single Biological aspect of the paper. 

## OK, I got my paper, where do I start? 

### Read the title and abstract, see if you can answer these questions:
#### What is the contribution of the paper?
* Is it a specific Biological finding? Like: A gene being associated to a disease. 
* Is it introducing a new analytical method (which is supposedly better than the previously known methods)? Like: a new clustering method for gene expression data. 
* Is it introducing a new form of molecular experiment, which could potentially result in new findings? Like: a new method for Single Cell experiments. 

#### What are the questions (if any) that the authors were trying to answer in their research? 
Like: A study on the genetics of Autism. _The authors were trying to find (or claim to have found) a group of genes associated with ASD._

At this point you should have a brief idea of what is the paper about. You probably could even guess what kind of data is used in the study (or it is mentioned in the abstract), or what are the key analytical processes.  

### Through the main body and the methods section (or even supplementary), learn about the data they used  
* Do they use any expression datasets (i.e. RNA-seq or microarray) in their study? Do they have an expression matrix for genes or transcripts? What is the size of the data matrix? How many factors do they have? Look for some of the QC steps in the supplementary or the main body of the paper. How would you compare that to those you learned from the course? 

* Look for other kinds of data in the paper. Are there Genomics, Proteomics or Methylation data? Are they using previously published datasets or they have generated the data themselves? Do they use a specific DataBase to complement their analysis and results? What is that database, how is the data collected there? 

### What kind of analysis do they do? 

It is common that you observe many analytical methods in a paper. Multiple methods could be applied (and compared with each other) on one dataset. Also, different methods might be used for different types of data in one paper. Look for the _key analysis/method used to answer the main questions in the paper_.

* Check out the figures or methods section for familiar terms like: p-value, clustering, PCA, etc. 
* What was the statistical test from which the p-value was obtained? What is the p-value exactly showing? 
* What is their clustering method? How is it different/similar to the clustering methods we see in the course? 
* What kind of model are they using for analyzing the data? 
* Did they use any Bootstrap/resampling? 
* Do you see any gene function/ network analysis? How could you comment/explain/criticize that? 

# NIREUS DNA Synthesis - Novel Iterative Restriction Enzyme digest Utilizing Sticky-end alignment DNA Synthesis

## Background

### Origins

NIREUS DNA Synthesis originally started as a continuation of my 9th-grade science fair project in 2016. Since then I have continued to work on the software, reworking and rewriting most of its elements. Now compiled in Python, this project successfully illustrates an iterative restriction enzyme and sticky-end alignment process for synthesizing DNA. 

### Rational

Stemming from the need to decrease the cost of synthetic synthesis, the rationale behind this research begins with the introduction of a hypothesized catabolic synthesis method (referred to as NIREUS DNA Synthesis): a process that iteratively digests DNA into fragments using many endonucleases to synthesize a gene. 

This process is an extension of the BioBricks™ standard proposed by Dr. Thomas Knight. In their publication, BioBricks™ is defined as a process that “employs iterative restriction enzyme digestion and ligation reactions to assemble small basic parts into larger composite parts.” These larger components are built using fragments flanked by two restriction enzymes, Xbal and SpeI. One key difference between restriction synthesis and BioBricks™ is that the latter restricts digestion to those two enzymes: XbaI and SpeI. Conversely, NIREUS DNA Synthesis relies on a dataset of more than 1000 enzymes all with unique restriction cutting sites. In flanking subsequences with a greater number of available restriction enzymes, restriction synthesis was originally thought to build many more composite subsequences, which are homologous to those referenced by the BioBricks™ standard.

Since the project’s inception, the estimated cost of synthesis has shown to decrease rapidly as the simulation grew more representative of the proposed lab procedure. With costs now as low as $0.019/bp, this simulation is ready for its next stage of research, developmental testing.

## Methods

Given a query sequence and a larger DNA sequence (also known as the reference sequence), the algorithm simulates the process by first searching for fragments in the reference that sequentially appear in the query. These matches are considered subsequence matches, which is one of the three conditions required for a fragment to be used for the synthesis process. The length of these fragments varies with each search depending on what subsequences are available in the reference. It is more likely, for instance, to find a subsequence match of length eight compared to a subsequence match of length sixteen. Taking this into account, each subsequence search begins with parsing for subsequences of length t, and then continues by decrease this value until a subsequence match is found or t reaches a baseline value like s (t is set to 16 and s is set to 4 by default). Through this process, each iteration for a subsequence match ensures that the longest possible subsequence matches are found. 

An instance match, unlike a subsequence match, is based on whether two enzymes exist that can flank and digest the subsequence match that appears in the reference. This condition is critical because while there exist thousands of subsequence matches, only a select few may be applicable for the synthesis process-—in order for a subsequence to be included from the set, it needs to be removed by flanked enzymatic digest. Using Wikipedia's restriction enzyme database, the algorithm has access to over 1000 unique restriction enzymes and their respective restriction sites. Despite this relatively low number of enzymes compared to the number of sequences possibly available, many have versatile recognition sequences, cutting a variety of DNA sites. Access to these versatile enzymes is one of the reasons why this process can successfully complete a variety of genes. 

In addition to a subsequence match and an instance match, a third match is required: a ligation match. Some restriction enzymes cut DNA in a Z-like manner, where the product of the enzymatic digest contains nucleotide overhangs, commonly referred to as sticky-ends. Because many of the enzymes create these sticky-ends on digested fragments, a ligation match is needed to ensure that sequential fragments are compatible. For instance, a subsequence that ends in an overhang of AGTA on the top strand 5’ end needs to be followed with an overhang that begins with TCAT on the bottom 3’ end. Also, fragments cut with enzymes that produce blunt ends, or ends without any overhangs, need to be followed with fragments that have these produced blunt ends as well. This final end-to-end compatibility indicates whether a fragment is a complete match and whether the search process continues. 

## Next Steps

At this point, NIREUS DNA Synthesis is able to simulate this proposed DNA synthesis method. The following list outlines the next steps for NIREUS DNA Synthesis.

- The last cost analysis was completed almost four years ago, so the price of the method should be heavily inspected and reported.
- Design and perform tests related to the query-reference length ratio. 
- Design a final restriction map plot of the reference sequence that indicates what enzymes should be used where in order to get the desired fragments in lab
- Create a method for determining the order of applying enzymes for the digestion of the reference
- Implement a pre-built machine learning model for subsequence classification
- Build and implement a machine learning model for subsequence alignment prediction


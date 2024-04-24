### (1) sta141c.ipynb 
The main aim of this project is to predict whether an asteroid should be classified as hazardous based on explanatory variables. During the Exploratory analysis portion, one can observe that the explanatory varibles are highly correlated. VIF selection was used to select variables which were not highly correlated. Ridge regression was implemented from scratch to classify asteroids. Other machine learning methods like Logistic Regression, Boosting methods and tree based methods were also used.

A secondary aim of this project was a compare the performance of different matrix decomposition methods such as QR, LU and Cholesky decomposition when solving the normal equation $X^TX\beta = X^TY$

### (2) dna2prot.py

Command line python tool to convert DNA sequence to expected translated protein.
#### Running script.py
This script will generate an amino acid sequence of a protein from a given DNA sequence. It expects a .txt file containing a DNA sequence ordered in a 5' to 3' direction. This file should contain only the letters A, T, G, and C, and a sample input file (`input.txt`) is provided in the `input/` directory.

This will save the final protein sequence to an output .txt file in the `output/` directory. A sample `output.txt` file is provided in the `output/` directory. By default, the script writes to `output.txt`, but an optional flag can be specified by the user to give a custom output file name. Additionally, the detected genes are also outputted to a file (`output_genes.txt`) in the `output/` directory.

Usage:
```
python3 script.py <input_file>
python3 script.py -o <output_file> <input_file>
```

Additionally, the script contains options to perform mutation simulations on the DNA sequence. The -m or --mutation flag finds the effect of changing the base at a specific position on the resulting protein structure. The -p or --probability flag finds the probability that a mutation at a position impacts the protein product. These results are outputted to the terminal

### (3) yelp_scrape.ipynb

A python notebook contains the scapring, cleaning and visualization of data scraped from yelp.com. The difference metrics like the price, rating and health inspection violations (from sacremento county website) across different cuisines were analyzed.

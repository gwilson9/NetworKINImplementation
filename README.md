# NetworKINImplementation
Edit of NetworKIN3.0_release

###### Written by GMW - 20181005 ######

NetworKIN source code implementation for mouse and human phoshopeptide data sets. 

Source code was downloaded from http://networkin.info/download.shtml#changelog (NetworKIN3.0_release.zip)

Source code contains a file containing STRING interactions with Ensembl74 Human ENSP identifiers. It is necessary that the protein identifiers from your data be converted to corresponding Ensemb74 identifiers NetworKIN to return a kinase prediction. The script “UniprotIDsToEnsembl74.py” will perform this function. Input file must contain column "Majority Protein ID" with UniprotIDs, "Gene name" with human HGNC or Mouse MGI, "Amino Acid" with S, T, or Y site of phosphorylation, and "Position within Protein" filled with the residue number of modification. This will write out the "PhosphoSites.txt" which will be one of the input files to NetworKIN. 

Usage: ./UniprotIDsToEnsembl74 -o ("10090" or "9606") -r (results file.csv)

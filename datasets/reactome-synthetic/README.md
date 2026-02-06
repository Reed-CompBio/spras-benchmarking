# Reactome Pathways (Signal Transduction)

This dataset contains 33 pathways listed under the "Signal Transduction" category of Reactome (downloaded February 4, 2026). Each pathway in the [list of Signal Transduction pathways](https://reactome.org/PathwayBrowser/#/R-HSA-162582) contains a unique standard identifier. 

TODO: Check why there are 33 and not 17 (like the list indicates).

## Get BioPAX-Formatted Pathways

The script `scripts/extract-reactome-pathways.py` uses Reactome's RESTful API to get BioPAX Level 3 files for the 33 pathways, and then we use [PaxTools](https://www.biopax.org/Paxtools/) to convert BioPAX `.owl` files into extended SIF formats. A final script will eventually convert the SIF-formatted files into processed files for SPRAS inputs.

You must be in the `scripts/` directory.

```
python extract-reactome-pathways.py
```

## PaxTools to Convert to SIF-Formatted Pathways

In the `scripts/` directory, [Download the "fat" PaxTools JAR](https://sourceforge.net/projects/biopax/files/latest/download). Version 6.0.0.

```
java -jar paxtools-6.0.0.jar
```

If needed install Java (tested with `openjdx@17` installed via Homebrew). A Paxtools call to convert one pathway to SIF is:

```
java -jar paxtools-6.0.0.jar toSIF ../raw/R-HSA-157118__Signaling_by_NOTCH.owl ../raw/R-HSA-157118__Signaling_by_NOTCH.sif -extended seqDb=uniprot,hgnc,refseq
```

This specifies UniProt identifiers, then HGNC identifiers if UniProt isn't found, then RefSeq identifers if neither UniProt nor HGNC identifiers are found.

To convert all `*.owl` files into `*.sif` files, 

```
ls ../raw/*.owl | sed 's/.owl//g' | awk '{print "java -jar paxtools-6.0.0.jar toSIF "$1".owl "$1".sif -extended seqDb=uniprot,hgnc,refseq"}' | bash
```

TODO: Need to add molecules to a blacklist. There's a recommended one somewhere on PathwayCommons.

## If needed, install java via Homebrew

These are Anna's notes, will be removed if we decide to merge.

```
brew install openjdk@17
```

then symlink it

```
sudo ln -sfn /opt/homebrew/opt/openjdk@17/libexec/openjdk.jdk \
  /Library/Java/JavaVirtualMachines/openjdk-17.jdk
```

set `JAVA_HOME`

In ~.zshrc file:
```
export JAVA_HOME=$(/usr/libexec/java_home -v 17)
export PATH="$JAVA_HOME/bin:$PATH"
```

```
java -jar paxtools-6.0.0.jar toSIF ../raw/R-HSA-157118__Signaling_by_NOTCH.owl ../raw/R-HSA-157118__Signaling_by_NOTCH.sif -extended seqDb=uniprot,hgnc,refseq
```
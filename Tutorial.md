# Tutorial accompanying the chapter "Managing historical linguistic data for computational phylogenetics and computer-assisted language comparison"
 
# Tresoldi, Tiago; Rzymski, Christoph; Forkel, Robert; Greenhill, Simon J.; List, Johann-Mattis

This tutorial supplements the aforementioned chapter, intentionally following a similar progression of analysis. The aim is to facilitate understanding of the more theoretical topics of the chapter while allowing readers to see step-by-step examples of the input, intermediate, and output files of a common phylogenetic analysis.

The workflow uses the Python programming language and a number of libraries developed at the Max-Planck-Institute for the Science of Human History (MPI-SHH), in particular the `lingpy` library. All these libraries are free and open source. For more information on LingPy, we recommend readers see the tutorial published as [supplementary material](https://github.com/lingpy/lingpy-tutorial) to "[Sequence comparison in computational historical linguistics](https://academic.oup.com/jole/article/3/2/130/5050100)" by List et al. (2018). 

## 1 Getting Started

We assume readers are familiar with the essentials of command-line operation on either Unix-like systems (such as Linux and MacOS) or Windows systems. Likewise, we assume that the system is updated and that [Python](https://www.python.org/) is installed, along with its package manager `pip` and the `wheel` package. The latter might be necessary for installing the supporting libraries and, in those situations, it must be available *before* installing them.

We strongly recommend using virtual environments for reproducing the steps outlined here. Virtual environments are a practical solution for creating independent configurations for testing and experimenting, with no interference on the system-wide installation. More information on virtual environments can be found on the [Python documentation](https://docs.python.org/3/tutorial/venv.html). In most systems, a virtual environment can be created by calling `python -m venv env`. Once created, the environment can be activated with `source env/bin/activate` and deactivated with the `deactivate` command.

After all the dependencies have been installed, a virtual environment with the libraries required by this tutorial can be installed with:
 
```
python -m venv env
source env/bin/activate
pip install --upgrade wheel
pip install -r pip-requirements.txt
```

All the code snippets below assume that the necessary libraries have already been imported in Python. The full list of imports, which can be executed upon initialization and for testing that all libraries were properly installed, is:

```python
import lingpy
from lingpy.convert.strings import matrix2dst, write_nexus
from lingpy.compare.sanity import average_coverage, mutual_coverage_check, synonymy
from lingpy.evaluate.acd import bcubes
from tabulate import tabulate
```
## 2 Illustration of the Phylogenetic Data Life-Cycle

### 2.1 Data Handling (Stage 1)

As this tutorial is intended to illustrate the management of data for phylogenetic analyses, we will not discuss common operations of data collection, normalization, etc. A dataset for Kho-Bwa languages derived from Lieberherr and Bodt (2017) is distributed along with this document, and mirrors the final stage of data preparation that is expected. While `lingpy` can load directly from single-table wordlists, such as those produced with a spreadsheet program, we recommend converting data to a common format, [CLDF](https://cldf.clld.org/), before adapting this tutorial to other needs. More detailed explanations are given in `lingpy` and `CLDF` documentation, and the single-table wordlist that served as a source the CLDF data, which can be compared to the CLDF dataset, can be found as a supplementary material to Lieberherr and Bodt (2017).

#### 2.1.1 Creating and Curating Data in CLDF

While the details and advantages of CLDF have been discussed and presented in Forkel et al. (2018), we would like to briefly highlight the benefits of using CLDF as a basis for datasets in historical linguistics and, by extension, phylogenetic analyses. While, at a first glance, CLDF and the outer appearance of a dataset (e.g. information split across multiple files), might seem daunting or even plainly unnecessary, it gives researchers a powerful toolset via the expressiveness of relational data structures in well-supported serialization formats, namely CSV and JSON. The format allows for the creation of truly linked data (see for example the embedded information concerning languoids and concepts in the Kho-Bwa sample dataset), which at the same time exists independently of any other software stack, given that all the information are stored in plain (but structured) text files. The fact that multiple files are 'glued' together by a JSON file, which specifies how files are related as well as metadata information for every file and column, allows for unique identifiers for every record, and guarantees maximum flexibility with regards to research designs and experimental setups.

We support the creation as well as the conversion of CLDF-compatible data from a [variety of source materials](https://github.com/cldf/cookbook/tree/master/recipes) and CLDF can even be easily created from [within spreadsheet software like Microsoft Excel](https://github.com/cldf/cookbook/tree/master/recipes/excel). For the more systematic creation of CLDF datasets from different sources, we have created [`pylexibank`](https://github.com/lexibank/pylexibank/), a Python toolchain employed in the [Lexibank project](https://github.com/lexibank/). `pylexibank` makes it easy to programmatically create and curate CLDF datasets.

A first impression regarding the status of the data can be gained with the help of the Python `pycldf` package, which also makes a `cldf` command-line utility available. To validate the data, and to operate with CLDF datasets in general, we only need to refer to its JSON metadata file, as in: 
 
```
$ cldf validate cldf-khobwa/Wordlist-metadata.json
```

Note that `cldf validate` will generate no output if no issues are found.

Fundamental datasets statistics can be obtained with the `cldf stats` command:
 
```
$ cldf stats cldf-khobwa/Wordlist-metadata.json
```

Which for this dataset will output:
 
```
<cldf:v1.0:Wordlist at cldf-khobwa>
key            value
-------------  --------------------------------------------
dc:conformsTo  http://cldf.clld.org/v1.0/terms.rdf#Wordlist

Path            Type              Rows
--------------  --------------  ------
forms.csv       FormTable         1948
languages.csv   LanguageTable       20
parameters.csv  ParameterTable     100
cognates.csv    CognateTable      1948
```
 
We can quickly verify that the dataset contains 1948 forms, 100 concepts ("parameters") and 20 languages. Two of the languages included by Lieberherr and Bodt (2017) from other sources were excluded, as suggested by the authors.

#### 2.1.2 Loading Data in LingPy

Due to the CLDF integration, it is possible to load CLDF datasets as a single, in-memory wordlist with `lingpy`. Here we will store the data in the `wl` variable:

```python
# load the wordlist from CLDF
wl = lingpy.Wordlist.from_cldf("cldf-khobwa/Wordlist-metadata.json",
    columns=('parameter_id',
             'concept_name',
             'language_id',
             'language_name',
             'value',
             'form',
             'segments',
             'language_glottocode',
             'concept_concepticon_id',
             'cognacy',
             'cogid_cognateset_id'))
```

There are various ways to iterate over the data in the wordlist, such as directly indexing rows or by using the built-in `Wordlist.iter_rows()` method while using the column names listed in `wl.columns`. A quick check for consistency can be performed by tabulating the first rows of data:

```python
print(tabulate([wl[idx] for idx in range(1, 5)]))
```

```
---  ---  -------  -------    ---------  ----------------          ---
1sg  1SG  duhumbi  duhumbi    ga         g a                       101
2sg  2SG  duhumbi  duhumbi    naŋ        n a ŋ                     301
3sg  3SG  duhumbi  duhumbi    wɔj        w ɔ j                     501
ant  ant  duhumbi  duhumbi    kʰin-ʨʰɔk  kʰ i n + tɕʰ ɔ k          701
---  ---  -------  -------    ---------  ----------------          ---
```

The statistics reported by the command-line `cldf` utility can be confirmed by inspecting the in-memory wordlist. Besides confirming the number of languages, concepts, and forms, we can extract more advanced statistics such as minimal and mutual average coverage of concepts across the languages in the dataset [@List2018f]. These statistics provide useful information that can help us ascertain if we have an acceptable amount of comparable material for a phylogenetic analysis.

Basic statistics can be obtained by querying directly the wordlist object:

```python
# count number of languages, number of rows, number of concepts
print("Wordlist has {0} languages and {1} concepts across {2} rows.".format(wl.width, wl.height, len(wl)))
```

With `lingpy` providing functions to easily compute coverage statistics:

```python
# Check mutual coverage
for i in range(wl.height, 0, -1):
    if mutual_coverage_check(wl, i):
        print("Minimal mutual coverage is at {0} concept pairs.".format(i))
        break


# check avarage coverage
print('Avarage coverage is {0:.2f}.'.format(average_coverage(wl)))
```

Which, for this dataset, will report a reasonable minimal coverage of 87 concepts and 94% of average coverage. An additional important statistic about a dataset is the number of potential synonyms, which should be as low as possible. Once more, `lingpy` offers methods and functions to easily collect this information:

```python
# check for synonyms
synonyms = synonymy(wl)
num_synonyms = len(wl) - len(synonyms)
syn_ratio = 1 - (len(synonyms)/len(wl))
if num_synonyms == 0:
    print('Found {0} potential synonyms.'.format(num_synonyms))
else:
    print('Found {0} potential synonyms ({1:.2}%):'.format(num_synonyms, syn_ratio*100.0))
    for (language, concept), count in sorted(synonyms.items(), key=lambda x: x[1], reverse=True):
        if count > 1:
            print('{0:15}  {1:12}  {2}'.format(language, concept, count))
```

Which will inform that we have 13 potential synonyms, a ratio low enough not to impact our analyses.

```
Found 13 potential synonyms (0.67%):
rahung           3SG           2
rahung           to do/make    2
rahung           water         2
rahung           yesterday     2
khoitam          bone          2
khoitam          water         2
jerigaon         house         2
jerigaon         liver         2
dikhyang         bone          2
dikhyang         eye           2
dikhyang         flesh/meat    2
bichom           1SG           2
bichom           black         2
```

#### 2.1.3 Retrieving Existing Data 

In this tutorial, we are loading a CLDF dataset directly, but for most experiments, especially when exploring wide cross-linguistic questions, it is worth it to obtain the data as a resource built upon the `pylexibank` library. In this case, a dataset similar to the one here distributed (but with the two languages removed for our exploratory purposes) can be installed as a `lexibank_lieberherrkhobwa` package. Data packages are not released on PyPI, but they can be installed with `pip` from `git` repositories, such as with:

```
pip install -e git+https://github.com/lexibank/lieberherrkhobwa.git#egg=lexibank_lieberherrkhobwa
```

When installed with `pip`, a CLDF dataset can be loaded in different ways, the easiest of which is by loading from the metadata file:

```
from lexibank_lieberherrkhobwa import Dataset as DS

cldf_metadata = DS().cldf_dir.joinpath('cldf-metadata.json')
wl = lingpy.Wordlist.from_cldf(cldf_metadata,
    columns=('parameter_id',
             'concept_name',
             'language_id',
             'language_name',
             'value',
             'form',
             'segments',
             'language_glottocode',
             'concept_concepticon_id',
             'language_latitude',
             'language_longitude',
             'cognacy',
             'cogid_cognateset_id',),
)
```

Please remember that the dataset distributed along with this tutorial is derived and simplified from the original. The code for running analyses does not change, but the results will be different whether one loads from this teaching material or from the actual, complete original dataset.

### 2.2 Computer-Assisted Language Comparison (Stage 2)

#### 2.2.1 Sequence Alignment

Our dataset provides cognate sets coded by experts, but no information on sequence alignment. Sequence alignment is a necessary step for many methods of automatic cognate detection, and can be used to support other tasks of comparative linguistics such as identification of sound correspondences.

LingPy can be used to automatically align all forms belonging to same cognate set:

```python
alms = lingpy.Alignments(wl, transcription='form', ref="cogid", segments="segments")
alms.align(method='library')
```

We can explore the aligned forms with LingPy itself, such as for forms belonging to cognateset #5701, for `FAT (ORGANIC SUBSTANCE)`, which is shared by all languages in the dataset:

```python
# Print the alignment for the first cogid (index #1)
print(lingpy.SCA(alms.msa['cogid'][5701]))
```

Which will return a nicely formatted table:

```
e       -       j       u       -
a       ʒ       -       ɔː      -
a       z       -       ua      -
ə       -       j       ou      -
ɔ       -       j       ɔ       w
a       -       j       ɔː      -
e       -       j       o       -
ɔ       -       j       u       -
a       -       j       ɔː      -
a       -       j       ɔː      -
a       z       j       aː      -
a       ʤ       -       ua      -
e       -       j       o       -
a       -       j       ɔ       w
a       z       j       aː      -
a       -       j       ɔː      -
a       z       -       ua      -
a       -       j       ɔː      -
e       -       j       u       -
e       -       j       u       -
```

The alignment object is a wordlist which carries information from its source, and can written to disk as a single-table or CLDF dataset. The alignment can therefore be easily exported into other tools, for example Edictor (mentioned in #3.4), which is designed for exploration and manipulation of this kind of data. A tab-separated output, for example, can be generated with:

```python
alms.output('tsv', filename='khobwa-aligned', ignore="all", prettify=False)
```

#### 2.2.2 Cognate Identification and Evaluation
 
Even though this dataset has expert cognate decisions already, we can also automatically infer cognates as well. Automatic cognate inference is useful in cases where there are no cognates, or for benchmarking the accuracy of the different available methods.

Here, we will use the pre-processed data from the file `khobwa-aligned.tsv` to create a `LexStat` scorer and perform detection using the three most common algorithms, `edit-dist`, `sca`, and `lexstat` (added as columns `editid`, `scaid`, and `lexstatid`, respectively):

```python
lex = lingpy.LexStat('khobwa-aligned.tsv', segments='tokens', check=True)
lex.get_scorer(runs=10000)
```

It is normal for this to take a while and to obtain a long output detailing the alignment and random correspondence. Once these steps are finished, we can run the methods for automatic cognate detection, each with its reference threshold and clustering method:

```
cognate_detect = [
    ('edit-dist', 'editid', 0.75, 'upgma'),
    ('sca', 'scaid', 0.45, 'upgma'),
    ('lexstat', 'lexstatid', 0.55, 'infomap')
]

for method, ref, threshold, clustering in cognate_detect:
    lex.cluster(method=method, threshold=threshold, ref=ref, cluster_method=clustering)
```

Once all the methods have run, we can evaluate their performance:


```python
# Compute overall scores in terms of precision, recall, and f-score, and present it
ret = []
for column in ['editid', 'scaid', 'lexstatid']:
    precision, recall, fscore = bcubes(lex, 'cogid', column, pprint=False)
    ret.append([
        column, '%.3f' % precision, '%.3f' % recall, '%.3f' % fscore
    ])

print(tabulate(ret, headers=['Precision', 'Recall', 'F-score'], tablefmt='pipe'))
```

Which will detail how, for this dataset, LexStat performs better than SCA, which performs better than Edit Distance (a performance that is replicated in most datasets):

```
|           |   Precision |   Recall |   F-score |
|:----------|------------:|---------:|----------:|
| editid    |       0.871 |    0.789 |     0.828 |
| scaid     |       0.92  |    0.766 |     0.836 |
| lexstatid |       0.916 |    0.812 |     0.861 |
```

The reported numbers might be slightly different due to the random permutations used in the detection method.

LingPy also allows for the detection of partial cognates, which can be used to compute cognate identifications based on shared or non-shared partial information. We will not explore this option here, as the resulting data does not change, but it is an option worth exploring for many datasets and better detailed in the aforementioned tutorial.

### 2.2.3 Manual Exploration and Manipulation 

Data can be explored and manipulated with the [EDICTOR tool](http://edictor.digling.org), a web-based interface for the creation and curation of etymological dictionaries. EDICTOR essentially supports to edit and correct cognate judgments manually, and also allows to inspect and correct phonetic alignments. In addition, one can use EDICTOR to inspect cognate set distributions, phoneme inventories, and morphological structures (full and partial colexifications), also allowing to export the data to Nexus format.

![Edictor](edictor.png)

### 2.3 Exploratory Data Analysis (Stage 3)

#### 2.3.1 Distance Measures

The cognate set of each form, here reported in the `cogid` column, is the main material for phylogenetic analyses. For initial data exploration base in language distance, we can use LingPy's built-in function to obtain the distance matrix used by most methods, saving it to a `khobwa.dst` file.

```python
wl.get_distances(ref='cogid')
wl.output('dst', filename='khobwa')
```

The distance matrix file is a pure-textual file with the number of languages in the first row, preceded by a space, followed by one language and its distances per line, separated by a space (in which case the language name should not be longer than ten characters, or by a tabstop, which is not reflecting the official standard, but accepted by many software packages). The file can be checked with standard shell commands or passed to other programs like [SplitsTree](http://splitstree.org):

```
$ head -n 3 khobwa.dst
 20
bichom  0.0000  0.5506  0.5618  0.0909  0.5730  0.6092  0.1461  0.5843  0.6250  0.6180  0.5169  0.5618  0.0787  0.6180  0.5227  0.6180  0.5281  0.6250    0.0449  0.0674
bulu    0.5506  0.0000  0.3300  0.5567  0.5600  0.5204  0.5730  0.5400  0.5354  0.5000  0.1700  0.3400  0.5618  0.5200  0.2062  0.5400  0.2700  0.5455    0.5618  0.5618
```


#### 2.3.2 Distance-Based Phylogenetic Trees

Distance matrices can be used for additional methods of data exploration, such as bulding distance-based trees with algorithms like Neighbor joining (NJ) and Unweighted Pair Group Method with Arithmetic-mean (UPGMA). Many software packages, including _SplitsTree_, can use such matrices, with LingPy already providing methods for easy computation and display of such trees:

```python
tree = lingpy.Tree(wl.get_tree(ref='cogid', tree_calc='upgma', force=True))
print(tree.asciiArt(compact=True))
```

Which will print a nice ASCII representation:

```
                    /edge.0-- /-duhumbi
          /edge.6--|          \-khispi
         |          \edge.5-- /-shergaon
         |                    \edge.4-- /-rupa
         |                              \edge.3-- /-jerigaon
         |                                        \edge.2-- /-khoina
-root----|                                                  \edge.1-- /-khoitam
         |                                                            \-rahung
         |                    /-kaspi
         |          /edge.11-|          /-namphri
         |         |          \edge.10-|          /edge.7-- /-bichom
         |         |                    \edge.9--|          \-singchung
          \edge.17-|                              \edge.8-- /-dikhyang
                   |                                        \-wangho
                   |          /edge.13- /-saria
                    \edge.16-|          \edge.12- /-chayangtajo
                             |                    \-lasumpatte
                              \edge.15- /-bulu
                                        \edge.14- /-kojorojo
                                                  \-rawa
```

The tree itself can be saved in standard Newick format for use in tree visualisation and manipulation software:

```python
# save tree to file:
tree.writeToFile("lieberherrkhobwa-upgma.trees")
```

### 2.4 Bayesian Phylogenetic Analysis (Stage 4)

Data can also be exported to perform Bayesian analysis with tools such as `BEAST2`, `MrBayes`, `RevBayes`, `APE`, etc., whose usage will not be discussed here. From the point of view of data-management, almost all software for phylogenetic analysis use data in some dialect of the NEXUS format, which can be generated from a wordlist with LingPy:

```python
# Save output for splitstree
nex = write_nexus(wl, mode="SPLITSTREE", filename="lieberherrkhobwa.nex")

# Save output for BEAST
nex = write_nexus(wl, mode="BEAST", filename="lieberherrkhobwa-beast.nex")
```

One important software of the phylogenetic toolbelt is `beastling`, which facilitates the interfacing of linguistic data with `BEAST2` and had full CLDF integration. As such, we can export our wordlist data in CLDF format and use it with the expanding range of software supporting CLDF:

```python
from lingpy.convert.cldf import to_cldf
to_cldf(wl, path='cldf-new')
```

### 2.5 Data Sharing (Stage 5)

Being able to export wordlists, with any additional information such as cognate sets coming from automatic cognate judgment, allows users to easily export, deploy, and share datasets. Both wordlists and CLDF datasets can be either provided as-is to users or, preferably, made available via on-line services such as GitHub and Zenodo. For CLICS, the Database of Cross-Linguistic Colexifcations, we make heavy use of this and automatically publish all curated datasets from GitHub to the [CLICS Zenodo Community](https://zenodo.org/communities/clics/).

An interesting possibility is to rely on the infrastucture of `lexibank`, to deploy datasets as Python packages that can then be conveniently installed from the command line, and [`clld`](https://clld.org), to deploy them as online services. This makes publishing and serving curated data straightforward, and also has the advantage of facilitating anonymous sharing of data for review purposes, deploying it as a Python package on GitHub and using it by means of on-line services such as the Open Science Framework.

## 3 Conclusion

Given the complexity of the topic, the current tutorial can only provide an incomplete introduction to the fascinating world of data managment in computational phylogenetics and computer-assisted language comparison. We hope, however, to have shown enough to ease the the entry for readers with little experience in data managment and phylogenetic analysis.  If you have questions on the management of linguistic data for phylogenetic analyses or CLDF, don't hesitate to contact us at <lexibank@shh.mpg.de>. Questions about LingPy can be posted either as an issue on its [GitHub repository](https://github.com/lingpy/lingpy) or sent to <info@lingpy.org>.

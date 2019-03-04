# The tissue analysis: `one_tissue_analysis.py`

This file contains several functions to be able to analyse and represent the same results as TheGeneAnalysis file, but looking at all the genes of a whole tissue file.

These are the all the possible arguments we can use:

| Argument | Type | Default | Description |
|:-|:-:|:-:|:-|
| `--data` | string | - | The *.csv.gz* or *.csv* file where the tissue data is saved. |
| `--minexp` | float | 0.8 | The minimum expression you want to be explained. Imagine you have a gene with 5 isoforms, and the proportion of each one is 0.18, 0.72, 0.05, 0.02 and 0.13 in one of the samples, respectively. The second isoform explains the 72% of the expression of the gene, which is very high, but maybe you need a minimum of 80% to be explained. In this case, this is what you should write: `--minexp 0.8`. |
| `--minsamps` | integer | 10 | The minimum samples where a gene has to be expressed to take it into account. This threshold is useful when you have a big number of samples. |
| `--out_tsdir` | string | Working directory | The name of the directory where the tissue summary file will be saved. This argument is useful in case you decide to perform various analysis with different thresholds. |
| `--out_tsfile` | string | *tissue* _ *minexp*_*minsamps*.csv | The name of the file where all the results will be saved (tissue summary file). |
| `--in_tsfile` | string | - | The name of the tissue summary file (saved before). |
| `--savefile` | boolean | `True` | In case of the tissue statistics function, do you want to create a file with the results? Then, write `--savefile` in your call and it will save them in a *.csv* file. If you do not write anything, it will show the results in the terminal. |
| `--out_statsdir` | string | Working directory | The name of the directory where the tissue summary file will be saved. |
| `--out_statsfile` | string | *tissuefile*_statistics.csv | The name of the file where all the statistics will be saved (tissue statistics file). |
| `--in_statsfile` | string | - | The name of the file where all the statistics has been saved (tissue statistics file). |
| `--genetype` | string, list of strings | All | The gene types you are interested in. |
| `--drop_tsfile` | boolean | `True` | Do you want to remove the tissue summary file? Then, write `--drop_tsfile` in your call. This is useful when you have a big number of tissue summary files, that could take up gigabytes of memory, and you only want the statistics, not the whole summaries. **NOTE**: In case of the tissue statistics function, this argument just works if you have saved the tissue statistics file (so, if you have used in the same call the argument `--savefile`. |
| `--seqnumsamps` | integer, list of integers | - | A sequence of minimum number of samples (`minsamps`) you want to compare. If you write `--seqnumsamps 1 2 3`, the function you call will analyse the data with a `minsamps` of 1, a `minsamps` of 2 and a `minsamps` of 3. |
| `--seqexp` | float | 0.05 | The difference between the `minexp` thresholds you want to compare. As the default is 0.05, the analysis will be done taking into account that `minexp` is 0, 0.05, 0.1, 0.15, ... until 0.95, so there will be 20 analysis for each `minsamps` you wrote in the `--seqnumsamps` argument; if you wrote `--seqexp 0.01`, there will be 100 analysis for each `minsamps` you wrote in the `--seqnumsamps` argument, and so on. |
| `--out_thresdir` | string | *tissue*_DiffThres_Summaries | The name of the directory where all the tissue summary files will be saved. |
| `--in_thresdir` | string | - | The name of the directory where all the different tissue summary files has been saved. Inside this folder, there cannot be anything else but the tissue summary files. |
| `--ncpus` | integer | 1 | The number of CPUs you will use. |
| `--num_cores` | integer | 1 | The number of cores you will use. |
| `--plotfile` | string | *statisticsfile*_*plottype*.png | The name of the file where the plot will be saved. It can be *.pdf* or *.png*. |
| `--expressed` | boolean | `True` | Do you want to draw just the results of the expressed genes? Then, write `--expressed` in your call.  |
| `--samplots` | integer, list of integers | All | For how many samples do you want to see the histograms? |

In some of the functions, some of these arguments are mandatory. We will see this for each function.

## Numeric functions

These are the functions that perform the analysis:

- Tissue summary (TSu)
- Tissue statistics (TSt)
- Tissue different thresholds summaries (TDSu)
- Tissue different thresholds statistics (TDSt)

#### Tissue summary function: `-TSu`

This function uses the gene classification function of the `gene_analysis.py` program to classifies each gene of a whole tissue into Monoform, Biform, Triform, Multiform, NotExpressed and FewSamples, and creates a *.csv* file with this information.

Mandatory arguments:

* `--data`

Optional arguments:

* `--out_tsdir`
* `--out_tsfile`
* `--minexp`
* `--minsamps`

Let's see how it works:

```bash
$ python one_tissue_analysis.py -TSu --data AllTissues_Initial/SMTS_Fallopian_Tube.csv.gz --minexp 0.9 --minsamps 1
```

The resulting file is like this:

| GeneId | GeneName | GeneType |MinimumExpression | MinimumSamples | ExpressedSamples | NumberOfExpressedSamples | GeneClassification | Transcripts | TranscriptTypes | Mean | CumulativeMean | NewMean | NewCumulativeMean |
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| ENSG00000000971.11 | CFH | protein_coding | 0.9 | 1 | ['GTEX-OHPK-2326-SM-3MJH2', 'GTEX-S32W-1326-SM-4AD5Q', 'GTEX-S341-0826-SM-4AD73', 'GTEX-SE5C-0926-SM-4BRUF', 'GTEX-T5JW-0326-SM-4DM6J', 'GTEX-T6MO-1026-SM-4DM72', 'GTEX-U3ZN-1126-SM-4DXUL'] | 7 | Triform | ['ENST00000367429.4', 'ENST00000439155.2', 'ENST00000359637.2'] | ['protein_coding', 'protein_coding', 'protein_coding'] | [0.5935667812018571, 0.27960656394042854, 0.07837390839278571] | [0.5935667812018571, 0.8731733451422856, 0.9515472535350713] | [0.6237911769454547, 0.29384411851504866, 0.08236470453949671] | [0.6237911769454547, 0.9176352954605034, 1.0] |
| ENSG00000004848.6 | ARX | protein_coding | 0.9 | 1 | ['GTEX-S341-0826-SM-4AD73'] | 1 | Monoform | ['ENST00000379044.4'] | ['protein_coding'] | [1.0] | [1.0] | [1.0] | [1.0] |
| ENSG00000000005.5 | TNMD | protein_coding | 0.9 | 1 | [] | 0 | NotExpressed | - | - | - | - | - | - |

In this example, the function has created the file *SMTS_Fallopian_Tube_0.900_0001.csv* inside the working directory.

#### Tissue statistics function: `-TSt`

> **ATTENTION**: It is mandatory to have done the **Tissue summary function** step to be able to do this one.

With this function, we perform an anaylisis of the results of the tissue summary function.

Mandatory arguments:

* `--in_tsfile`

Optional arguments:

* `--savefile`
* `--out_statsdir`
* `--out_statsfile`
* `--genetype`
* `--drop_tsfile`

Let's see an example of usage:

```bash
$ python one_tissue_analysis.py -TSt --in_tsfile SMTS_Fallopian_Tube_0.900_1.csv

GeneType     MinimumExpression MinimumSamples 3prime_overlapping_ncrna IG_C_gene IG_C_pseudogene IG_D_gene IG_J_gene IG_J_pseudogene IG_V_gene IG_V_pseudogene Mt_rRNA Mt_tRNA TR_C_gene TR_D_gene TR_J_gene TR_J_pseudogene TR_V_gene TR_V_pseudogene antisense lincRNA miRNA misc_RNA polymorphic_pseudogene processed_transcript protein_coding pseudogene rRNA sense_intronic sense_overlapping snRNA snoRNA  Total
Monoform                 0.900              1                       10        13               1         0         0               0       109               9       2       0         4         0         0               0        21               1      1460    1244    12      121                      5                  121           5218       1627    7            205                90    66     81  10427
Biform                   0.900              1                        2         0               0         0         0               0         1               1       0       0         0         0         0               0         0               0       222     174     0        0                      6                   77           3905        224    0              4                 5     0      0   4621
Triform                  0.900              1                        0         0               0         0         0               0         0               0       0       0         0         0         0               0         0               0        82      62     0        0                      2                   36           2516         70    0              1                 2     0      0   2771
Multiform                0.900              1                        0         0               0         0         0               0         0               0       0       0         0         0         0               0         0               0        47      34     0        0                      2                   38           3752         55    0              1                 1     0      0   3930
NotExpressed             0.900              1                        9         1               8        37        18               3        28             177       0      22         1         3        74               4        76              26      3465    5600  3043     1913                     30                  243           4954      11955  520            531               104  1850   1376  36071
FewSamples               0.900              1                        0         0               0         0         0               0         0               0       0       0         0         0         0               0         0               0         0       0     0        0                      0                    0              0          0    0              0                 0     0      0      0
```

Here we decided not to save the statistics file. This was to show which is the answer of this function and what will be the saved data in case you decide to save it.

Then, if we write...

```bash
$ python one_tissue_analysis.py -TSt --in_tsfile SMTS_Fallopian_Tube_0.900_1.csv --savefile
```

... a new file called *SMTS_Fallopian_Tube_0.900_1_statistics.csv* will be created inside your working directory, but we decided not to drop the tissue summary file yet.

Let's see another example:

```bash
$ python one_tissue_analysis.py -TSt --in_tsfile SMTS_Fallopian_Tube_0.900_1.csv --genetype protein_coding pseudogene


GeneType     MinimumExpression MinimumSamples protein_coding pseudogene  Total
Monoform                 0.900              1           5218       1627   6845
Biform                   0.900              1           3905        224   4129
Triform                  0.900              1           2516         70   2586
Multiform                0.900              1           3752         55   3807
NotExpressed             0.900              1           4954      11955  16909
FewSamples               0.900              1              0          0      0
```

Now, we just wanted to see (but not save, so ~~`--savefile`~~) the results of the protein coding genes and the pseudogenes (so `--genetype protein_coding pseudogene`), without losing the tissue summary file (so ~~`--drop_tsfile`~~).

If we want to save these results in a file called, for instance, *FallopianTube_protcod_pseudo.csv* in a new folder called *Fallopian_Tube_Statistics*, we write...

```bash
$ python one_tissue_analysis.py -TSt --in_tsfile SMTS_Fallopian_Tube_0.900_1.csv --genetype protein_coding pseudogene --out_statsdir Fallopian_Tube_Statistics --out_statsfile FallopianTube_protcod_pseudo.csv --savefile
```

#### Tissue different threshold summaries function: `-TDSu`

This function does exactly the same as the tissue summary one, but adding the possibility of obtaining the results for more than one threshold, for both the *minexp* argument and the *minsamps* one.

Mandatory arguments: 

* `--data`
* `--seqnumsamps`

Optional arguments:

* `--out_thresdir`
* `--seqexp`
* `--ncpus`
* `--num_cores`

Let's see its functioning:

```bash
$ python one_tissue_analysis.py -TDSu --data AllTissues_Initial/SMTS_Fallopian_Tube.csv.gz --seqnumsamps 1 3 5 7
```

This call has created a folder called *SMTS_Fallopian_Tube_DiffThres_Summaries* and saved inside `4 x 20 = 80` tissue summary files, one per each pair of `minexp` and `minsamps` values:

```bash
$ ls SMTS_Fallopian_Tube_DiffThres_Summaries/
SMTS_Fallopian_Tube_0.000_1.csv  SMTS_Fallopian_Tube_0.150_5.csv  SMTS_Fallopian_Tube_0.350_1.csv  SMTS_Fallopian_Tube_0.500_5.csv  SMTS_Fallopian_Tube_0.700_1.csv  SMTS_Fallopian_Tube_0.850_5.csv
SMTS_Fallopian_Tube_0.000_3.csv  SMTS_Fallopian_Tube_0.150_7.csv  SMTS_Fallopian_Tube_0.350_3.csv  SMTS_Fallopian_Tube_0.500_7.csv  SMTS_Fallopian_Tube_0.700_3.csv  SMTS_Fallopian_Tube_0.850_7.csv
SMTS_Fallopian_Tube_0.000_5.csv  SMTS_Fallopian_Tube_0.200_1.csv  SMTS_Fallopian_Tube_0.350_5.csv  SMTS_Fallopian_Tube_0.550_1.csv  SMTS_Fallopian_Tube_0.700_5.csv  SMTS_Fallopian_Tube_0.900_1.csv
SMTS_Fallopian_Tube_0.000_7.csv  SMTS_Fallopian_Tube_0.200_3.csv  SMTS_Fallopian_Tube_0.350_7.csv  SMTS_Fallopian_Tube_0.550_3.csv  SMTS_Fallopian_Tube_0.700_7.csv  SMTS_Fallopian_Tube_0.900_3.csv
SMTS_Fallopian_Tube_0.050_1.csv  SMTS_Fallopian_Tube_0.200_5.csv  SMTS_Fallopian_Tube_0.400_1.csv  SMTS_Fallopian_Tube_0.550_5.csv  SMTS_Fallopian_Tube_0.750_1.csv  SMTS_Fallopian_Tube_0.900_5.csv
SMTS_Fallopian_Tube_0.050_3.csv  SMTS_Fallopian_Tube_0.200_7.csv  SMTS_Fallopian_Tube_0.400_3.csv  SMTS_Fallopian_Tube_0.550_7.csv  SMTS_Fallopian_Tube_0.750_3.csv  SMTS_Fallopian_Tube_0.900_7.csv
SMTS_Fallopian_Tube_0.050_5.csv  SMTS_Fallopian_Tube_0.250_1.csv  SMTS_Fallopian_Tube_0.400_5.csv  SMTS_Fallopian_Tube_0.600_1.csv  SMTS_Fallopian_Tube_0.750_5.csv  SMTS_Fallopian_Tube_0.950_1.csv
SMTS_Fallopian_Tube_0.050_7.csv  SMTS_Fallopian_Tube_0.250_3.csv  SMTS_Fallopian_Tube_0.400_7.csv  SMTS_Fallopian_Tube_0.600_3.csv  SMTS_Fallopian_Tube_0.750_7.csv  SMTS_Fallopian_Tube_0.950_3.csv
SMTS_Fallopian_Tube_0.100_1.csv  SMTS_Fallopian_Tube_0.250_5.csv  SMTS_Fallopian_Tube_0.450_1.csv  SMTS_Fallopian_Tube_0.600_5.csv  SMTS_Fallopian_Tube_0.800_1.csv  SMTS_Fallopian_Tube_0.950_5.csv
SMTS_Fallopian_Tube_0.100_3.csv  SMTS_Fallopian_Tube_0.250_7.csv  SMTS_Fallopian_Tube_0.450_3.csv  SMTS_Fallopian_Tube_0.600_7.csv  SMTS_Fallopian_Tube_0.800_3.csv  SMTS_Fallopian_Tube_0.950_7.csv
SMTS_Fallopian_Tube_0.100_5.csv  SMTS_Fallopian_Tube_0.300_1.csv  SMTS_Fallopian_Tube_0.450_5.csv  SMTS_Fallopian_Tube_0.650_1.csv  SMTS_Fallopian_Tube_0.800_5.csv
SMTS_Fallopian_Tube_0.100_7.csv  SMTS_Fallopian_Tube_0.300_3.csv  SMTS_Fallopian_Tube_0.450_7.csv  SMTS_Fallopian_Tube_0.650_3.csv  SMTS_Fallopian_Tube_0.800_7.csv
SMTS_Fallopian_Tube_0.150_1.csv  SMTS_Fallopian_Tube_0.300_5.csv  SMTS_Fallopian_Tube_0.500_1.csv  SMTS_Fallopian_Tube_0.650_5.csv  SMTS_Fallopian_Tube_0.850_1.csv
SMTS_Fallopian_Tube_0.150_3.csv  SMTS_Fallopian_Tube_0.300_7.csv  SMTS_Fallopian_Tube_0.500_3.csv  SMTS_Fallopian_Tube_0.650_7.csv  SMTS_Fallopian_Tube_0.850_3.csv
```

You might add `--seqexp 0.02`, for instance, and take `4 x 50 = 200` files, and so on.

It will take a while to create all the summaries, depending on the memmory taked up by the original file (`--data`) and the `--seqnumsamps` and `--seqexp` arguments you wrote. You can use the `--ncpus` and/or `--numcores` arguments to reduce this time.

#### Tissue different threshold statistics function: `-TDSt`

> **ATTENTION**: It is mandatory to have done the **Tissue different threshold summaries function** step to be able to do this one.

This function does the same anaylisis as the tissue statistics function, but for all the tissue summary files in the folder you specify in the `--in_thresdir` argument.

Mandatory arguments: 

* `--in_thresdir`

Optional arguments:

* `--out_statsdir`
* `--out_statsfile`
* `--genetype`
* `--drop_tsfile`

Let's see its functioning:

```bash
$ python one_tissue_analysis.py -TDSt --in_thresdir SMTS_Fallopian_Tube_DiffThres_Summaries/
```

It has created a unique file called *SMTS_Fallopian_Tube_statistics.csv* in the working directory with the statistics results for each pair of thresholds:

| Classification | MinimumExpression | MinimumSamples | 3prime_overlapping_ncrna | IG_C_gene | ... | protein_coding | pseudogene | ... | Total |
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| Monoform | 0 | 1 | 12 | 13 | ... | 15391 | 1976 | ... | 21749 |
| Biform | 0 | 1 | 0 | 0 | ... | 0 | 0 | ... | 0 |
| Triform | 0 | 1 | 0 | 0 | ... | 0 | 0 | ... | 0 |
| Multiform | 0 | 1 | 0 | 0 | ... | 0 | 0 | ... | 0 |
| NotExpressed | 0 | 1 | 9 | 1 | ... | 4954 | 11955 | ... | 36071 |
| FewSamples | 0 | 1 | 0 | 0 | ... | 0 | 0 | ... | 0 |
|...|...|...|...|...|...|...|...|...|...|
| Monoform | 0.85 | 3 | 5 | 13 | ... | 4939 | 822 | ... | 7754 |
| Biform | 0.85 | 3 | 1 | 0 | ... | 3878 | 192 | ... | 4446 |
| Triform | 0.85 | 3 | 0 | 0 | ... | 2306 | 44 | ... | 2484 |
| Multiform | 0.85 | 3 | 0 | 0 | ... | 2721 | 40 | ... | 2846 |
| NotExpressed | 0.85 | 3 | 9 | 1 | ... | 4954 | 11955 | ... | 36071 |
| FewSamples | 0.85 | 3 | 6 | 0 | ... | 1547 | 878 | ... | 4219 |
|...|...|...|...|...|...|...|...|...|...|

As in case of the tissue statistics function, we can tell it which gene types we want as result, with the `--genetype` argument, change the name of the *.csv* file with the `--out_statsfile` argument, and so on:

```bash
$ python one_tissue_analysis.py -TDSt --in_thresdir SMTS_Fallopian_Tube_DiffThres_Summaries/ --out_statsfile FT_protcod_pseudo.csv --genetype protein_coding pseudogene
```

So as not to lose the first file we created, we add the argument `--out_statsfile` to create a new file called *FT_protcod_pseudo.csv*, and the resultant file is like this:

| Classification | MinimumExpression | MinimumSamples | protein_coding | pseudogene | Total |
|:-:|:-:|:-:|:-:|:-:|:-:|
| Monoform | 0 | 1 | 15391 | 1976 | 17367 |
| Biform | 0 | 1 | 0 | 0 | 0 |
| Triform | 0 | 1 | 0 | 0 | 0 |
| Multiform | 0 | 1 | 0 | 0 | 0 |
| NotExpressed | 0 | 1 | 4954 | 11955 | 16909 |
| FewSamples | 0 | 1 | 0 | 0 | 0 |
|...|...|...|...|...|...|...|...|...|...|
| Monoform | 0.85 | 3 | 4939 | 822 | 5761 |
| Biform | 0.85 | 3 | 3878 | 192 | 4070 |
| Triform | 0.85 | 3 | 2306 | 44 | 2350 |
| Multiform | 0.85 | 3 | 2721 | 40 | 2761 |
| NotExpressed | 0.85 | 3 | 4954 | 11955 | 16909 |
| FewSamples | 0.85 | 3 | 1547 | 878 | 2425 |
|...|...|...|...|...|...|...|...|...|...|

Note that, now, the totals are different than before. This is because, as we are just interested in two types of genes instead of all of theme, we recalculated the totals.

In addition, you may see that, now, we have not added any argument to tell the function that saves the statistics file. This is because there is much more information in this tables, with all the different thresholds.

## Plot functions

These functions return different plots to be able to visualise the results of the tissue statistics functions:

- Tissue genes barplot (TB)
- Tissue different thresholds barplots (TDB)

#### Tissue genes barplot function: `-TB`

> **ATTENTION**: It is mandatory to have done the **Tissue statistics function** step or the **Tissue different threshold statistics function** one to be able to do this one.

This function draws a stacked barplot with the statistics of the tissue statistics file, showing the results for the pair of `minexp` and `minsamps` you write, so that you can see the proportions of the gene classification of the tissue.

Mandatory arguments:

* `--in_statsfile` or `--in_thresfile`

Optional arguments:

* `--plotfile`
* `--expressed`
* `--minexp`
* `--minsamps`

You have 2 possible ways to call this, depending on how you have created the tissue statistics file.

If you have used the tissue statistics function, creating only the statistics for one unique pair of thresholds, you do not need to write which thresholds the function have to search for:

```bash
$ python one_tissue_analysis.py -TB --in_statsfile SMTS_Fallopian_Tube_0.900_1_statistics.csv
```

This call has created a *.png* file called *SMTS_Fallopian_Tube_0.9_1_AllGenesBarplot.png*:

![ ](/home/aserrano/PROJECTS/3.Mattia/1.- IsoformStats/plots/SMTS_Fallopian_Tube_0.9_1_AllGenesBarplot.png  "FT 0.9 1 AllGenes Barplot")

If you have used the tissue different thresholds statistics function, creating one unique file for the statistics of all pairs of thresholds, you should use the `--in_thresfile` instead of the `--in_statsfile` one; then, the functions takes the thresholds you write:

```bash
$ python one_tissue_analysis.py -TB --in_thresfile SMTS_Fallopian_Tube_statistics.csv --minexp 0.8 --minsamps 5
```

This call has created the file *SMTS_Fallopian_Tube_0.8_5_AllGenesBarplot.png* in the working directory:

![ ](/home/aserrano/PROJECTS/3.Mattia/1.- IsoformStats/plots/SMTS_Fallopian_Tube_0.8_5_AllGenesBarplot.png  "FT 0.8 5 AllGenes Barplot")

In this case, if you do not write the thresholds, the function use the defaults: `--minexp 0.8` and `--minsamps 10`.

In addition, you can draw this barplot only with the expressed genes (this is, the function would not take into account the genes classified neither as NotExpressed nor as FewSamples):

```bash
$ python one_tissue_analysis.py -TB --in_statsfile SMTS_Fallopian_Tube_0.900_1_statistics.csv --expressed
```

Now, the new *.png* file is called *SMTS_Fallopian_Tube_0.9_1_ExpressedGenesBarplot.png*:

![ ](/home/aserrano/PROJECTS/3.Mattia/1.- IsoformStats/plots/SMTS_Fallopian_Tube_0.9_1_ExpressedGenesBarplot.png  "FT 0.9 1 ExpressedGenes Barplot")

Or with the different threshold statistics file:

```bash
$ python one_tissue_analysis.py -TB --in_thresfile SMTS_Fallopian_Tube_statistics.csv --minexp 0.8 --minsamps 5 --expressed
```

![ ](/home/aserrano/PROJECTS/3.Mattia/1.- IsoformStats/plots/SMTS_Fallopian_Tube_0.8_5_ExpressedGenesBarplot.png  "FT 0.8 5 ExpressedGenes Barplot")

#### Tissue different thresholds barplots function: `-TDB`

> **ATTENTION**: It is mandatory to have done the **Tissue different threshold statistics function** step to be able to do this one.

This function draws a set of histograms to compare the results of the tissue different thresholds statistics file.

Mandatory arguments:

* `--in_thresfile`

Optional arguments:

* `--plotfile`
* `--expressed`
* `--genetype`
* `--samplots`

We can now use a new argument: `--samplots`. This argument refers to the `minsamps` one. For each `samplots` you write in the argument, a new plot will appear. If you do not write anything, it will show one plot per `minsamps` you have saved in the tissues different threshold statistics file before.

Let's see how the function works:

```bash
$ python one_tissue_analysis.py -TDB --in_thresfile SMTS_Fallopian_Tube_statistics.csv
```

This call has created a file called *SMTS_Fallopian_Tube_AllDiffThresHistograms.png* with the next plot:

Remember we created the *SMTS_Fallopian_Tube_statistics.csv* with the `--seqnumsamps 1 3 5 7`, so we have obtained 4 plots:

![ ](/home/aserrano/PROJECTS/3.Mattia/1.- IsoformStats/plots/SMTS_Fallopian_Tube_AllDiffThresHistograms.png  "DiffThres Histogram AllSamps")

If we add the `--samplots` argument with just 2 of this `seqnumsamps`, we just obtain 2 plots:

```bash
$ python one_tissue_analysis.py -TDB --in_thresfile SMTS_Fallopian_Tube_statistics.csv --samplots 3 5
```

![ ](/home/aserrano/PROJECTS/3.Mattia/1.- IsoformStats/plots/SMTS_Fallopian_Tube_3-5_AllDiffThresHistograms.png  "DiffThres Histogram 3-5 Samps")

Besides, we can add the `--expressed` argument, as before:

```bash
$ python one_tissue_analysis.py -TDB --in_thresfile SMTS_Fallopian_Tube_statistics.csv --expressed
```

![ ](/home/aserrano/PROJECTS/3.Mattia/1.- IsoformStats/plots/SMTS_Fallopian_Tube_ExpressedDiffThresHistograms.png  "DiffThress Histogram AllSamps ExpressedGenes")

In addition, we might add the `--genetype` argument in order to draw only the results of those gene types we want:

```bash
$ python one_tissue_analysis.py -TDB --in_thresfile SMTS_Fallopian_Tube_statistics.csv --genetype protein_coding pseudogene
```

This call has created the file *SMTS_Fallopian_Tube_protein_coding-pseudogene_AllDiffThresHistograms.png*, which have the next plot:

![ ](/home/aserrano/PROJECTS/3.Mattia/1.- IsoformStats/plots/SMTS_Fallopian_Tube_protein_coding-pseudogene_AllDiffThresHistograms.png  "DiffThres AllSamps PC&pseudo")
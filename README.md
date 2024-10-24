# sTPLS: identifying common and specific correlated patterns under multiple biological conditions

sTPLS could simultaneously identify common and specific correlated patterns between two types of features across multiple biological conditions. Given  paired data matrices $X^i \in R^{n_i\times p}$ and $Y^i \in R^{n_i\times q}\ (i=1,\cdots, M)$ from $M$ biological conditions, in which rows and columns respectively represent the same set of samples and two types of features, sTPLS integrates them to identify common comodules in which selected features demonstrate highly correlated patterns across multiple conditions and specific comodules in which selected features are correlated specifically in certain condition. 

### Packages

sTPLS is implemented in R. Examples could be found under the directory './sTPLS/examples' 

- Simulated data: refer to the [R Markdown file](https://htmlpreview.github.io/?https://github.com/Jinyu2019/sTPLS/blob/main/examples/GDSC_data/GDSC_application.html).
- Paired gene expression and drug response data from the GDSC dataset:  refer to the [R Markdown file](https://htmlpreview.github.io/?https://github.com/Jinyu2019/sTPLS/blob/main/examples/GDSC_data/GDSC_application.html).

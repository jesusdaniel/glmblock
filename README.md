# Simultaneous prediction and community detection for networks 

This repository contains R code that implements a block-structured regularization for regression problems with network-valued covariates to perform supervised community detection and regularized regression simultaneously (see the arxiv manuscript [Arroyo and Levina (2019)](https://arxiv.org/abs/2002.01645)).

## Methodology

Given a sample of N pairs of networks with <a href="https://www.codecogs.com/eqnedit.php?latex=n" target="_blank"><img src="https://latex.codecogs.com/gif.latex?n" title="n" /></a> aligned vertices (represented with their <a href="https://www.codecogs.com/eqnedit.php?latex=n\times&space;n" target="_blank"><img src="https://latex.codecogs.com/gif.latex?n\times&space;n" title="n\times n" /></a> adjacency matrices) and responses <a href="https://www.codecogs.com/eqnedit.php?latex=(A^{(1)},&space;Y_1),&space;\ldots,&space;(A^{(m)},&space;Y_n)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?(A^{(1)},&space;Y_1),&space;\ldots,&space;(A^{(m)},&space;Y_n)" title="(A^{(1)}, Y_1), \ldots, (A^{(m)}, Y_n)" /></a>, where <a href="https://www.codecogs.com/eqnedit.php?latex=Y_i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?Y_i" title="Y_i" /></a> is a real-valued variable, the method fits a regression model to predict <a href="https://www.codecogs.com/eqnedit.php?latex=Y_i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?Y_i" title="Y_i" /></a> using a linear combination of the entries of <a href="https://www.codecogs.com/eqnedit.php?latex=A^{(i))}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?A^{(i))}" title="A^{(i))}" /></a>. A matrix <a href="https://www.codecogs.com/eqnedit.php?latex=B\in\mathbb{R}^{n\times&space;n}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?B\in\mathbb{R}^{n\times&space;n}" title="B\in\mathbb{R}^{n\times n}" /></a> encodes the coefficients of the linear model to form a prediction rule of the form
<p style="text-align: center;"><a href="https://www.codecogs.com/eqnedit.php?latex=\widehat{Y}_i&space;=&space;f(\langle&space;A^{(i))},&space;B\rangle&space;&plus;&space;b)&space;=&space;f\left(\sum_{u=1}^n\sum_{v=1}^n&space;A^{(i)}_{uv}B_{uv}&space;&plus;&space;b&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\widehat{Y}_i&space;=&space;f(\langle&space;A^{(i))},&space;B\rangle&space;&plus;&space;b)&space;=&space;f\left(\sum_{u=1}^n\sum_{v=1}^n&space;A^{(i)}_{uv}B_{uv}&space;&plus;&space;b&space;\right&space;)" title="\widehat{Y}_i = f(\langle A^{(i))}, B\rangle + b) = f\left(\sum_{u=1}^n\sum_{v=1}^n A^{(i)}_{uv}B_{uv} + b \right )" /></a>,</p>
where <a href="https://www.codecogs.com/eqnedit.php?latex=f" target="_blank"><img src="https://latex.codecogs.com/gif.latex?f" title="f" /></a> is a link function and <a href="https://www.codecogs.com/eqnedit.php?latex=b\in\mathbb{R}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?b\in\mathbb{R}" title="b\in\mathbb{R}" /></a> is the intercept of the regression. Currently, only linear and logistic regression are implemented.

To enforce regularization and structure in the coefficients, the methods imposes a block-structured constraint by dividing the rows and columns of <a href="https://www.codecogs.com/eqnedit.php?latex=B" target="_blank"><img src="https://latex.codecogs.com/gif.latex?B" title="B" /></a> into <a href="https://www.codecogs.com/eqnedit.php?latex=K" target="_blank"><img src="https://latex.codecogs.com/gif.latex?K" title="K" /></a> different groups. This is analogous to the community detection problem, and by fitting a regularization of this form, the method can cluster the nodes of the networks into "supervised communities". Additionally, this constraint effectively reduces the number of different coefficients to deal with the high-dimensionality of the problem. The constraint has the form

<p style="text-align: center;">
<a href="https://www.codecogs.com/eqnedit.php?latex=B&space;=&space;ZCZ^\top" target="_blank"><img src="https://latex.codecogs.com/gif.latex?B&space;=&space;ZCZ^\top" title="B = ZCZ^\top" /></a>,</p>

where <a href="https://www.codecogs.com/eqnedit.php?latex=Z\in\{0,1\}^{n\times&space;K},\&space;\sum_{t=1}^KZ_{jk}=1\&space;\forall&space;j=1,\ldots,&space;K" target="_blank"><img src="https://latex.codecogs.com/gif.latex?Z\in\{0,1\}^{n\times&space;K},\&space;\sum_{t=1}^KZ_{jk}=1\&space;\forall&space;j=1,\ldots,&space;K" title="Z\in\{0,1\}^{n\times K},\ \sum_{t=1}^KZ_{jk}=1\ \forall j=1,\ldots, K" /></a> is a matrix with its rows indicating the community memberships for each node, and <a href="https://www.codecogs.com/eqnedit.php?latex=C\in\mathbb{R}^{K\times&space;K}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?C\in\mathbb{R}^{K\times&space;K}" title="C\in\mathbb{R}^{K\times K}" /></a> is a matrix of coefficients.

![supervised_community_detection](https://raw.githubusercontent.com/jesusdaniel/glmblock/master/img/B-ZCZT.png)


# References

Arroyo, J. and Levina, E. "Simultaneous prediction and community detection for networks with application to neuroimaging" [![arXiv shield](https://img.shields.io/badge/arXiv-2002.01645-red.svg?style=flat)](https://arxiv.org/abs/2002.01645)









# Improving Shoreline Forecasting Models with Multi-Objective Genetic Programming
by Al Najar, M., Almar, R., Bergsma, E. W., Delvit, J. M., & Wilson, D. G. (2023). This paper has been submitted for publication in Environmental Modelling & Software.

In this work, symbolic regression is performed using CGP and NSGA-II to evolve symbolic models for shoreline forecasting. The initial population of CGP individuals are initialized as ShoreFor, a hybrid shoreline forecasting model.

<!-- ![](outputs/plots/graphs/models/nf-generalmodel-graph-NSGA-ii-ind1.png)
*Caption for the example figure with the main results.* -->


<!-- ## Abstract

Given the current context of climate change and increasing population densities at coastal zones, there is an increasing need to be able to predict the development of our coasts. Recent advances in artificial intelligence allow for automatic analysis of observational data. This work makes use of Symbolic Regression, a type of Machine Learning algorithm, to evolve interpretable shoreline forecasting models. Cartesian Genetic Programming (CGP) is used in order to encode and im- prove upon ShoreFor, a shoreline prediction model. Coupled with NSGA-II, the CGP individuals are evaluated and selected during evolution according to their predictive skills at five coastal sites. This work presents a comparison between the CGP-evolved models and the base ShoreFor model. In addition to its ability to produce well-performing models, the work demonstrates the usefulness of CGP as a research tool to gain insight into the behaviors of shorelines at different points around the globe. -->

## Implementation
The Julia programming language is used throughout this work. All source code used to generate the results and figures in the paper are in the `scripts` folder. The evolved models, evolution logs and all result plots generated by the code are saved in `outputs`.

<!-- See the `README.md` files in each directory for a full description. -->

The `scripts` directory contains the following scripts:
* `create_shorefor_ind.jl`: contains the manually-defined template of the ShoreFor model... 
* `dataloader_monthly.jl`: data loading and any preprocessing.
* `evaluation_metrics.jl`: contains different functions for model evaluation.
* `evolve_5sites.jl`: launching the experiments as described in the article.
* `model_template_utils.jl`: different utilitary functions for transforming a template into a CGP individual.
* `result_analysis.jl`: evaluating evolved populations following the article.

This repository relies heavily on:
* [CartesianGeneticProgramming.jl](https://github.com/mahmoud-al-najar/CartesianGeneticProgramming.jl) for the implementations of CGP and NSGA-II, and
* [Cambrian.jl]() as the base evolutionary computation framework.

## Dependencies
This code was developed using Julia Version 1.6.1.

To install the dependencies listed in `Project.toml`, navigate to the project's main directory and activate the project's Julia environment:
```
(@v1.5) pkg> activate .
```
Then install using:
```
(CGP-ShoreFor) pkg> instantiate
```

Download [CartesianGeneticProgramming.jl](https://github.com/mahmoud-al-najar/CartesianGeneticProgramming.jl) and [Cambrian.jl](https://github.com/mahmoud-al-najar/Cambrian.jl), 

## Reproducing the results

## License

<!-- All source code is made available under a BSD 3-clause license. You can freely use and modify the code, without warranty, so long as you provide attribution to the authors. See `LICENSE.md` for the full license text. -->

<!-- The manuscript text is not open source. The authors reserve the rights to the article content, which is currently submitted for publication in the JOURNAL NAME. -->

<!-- * bullet1
* bullet2 -->

<!-- [link](https://github.com/dbader/readme-template) -->

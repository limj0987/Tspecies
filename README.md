# Tspecies

An R package for estimating species divergence times using Ks distribution variance and effective population size (Ne) correction.

## Description

Tspecies employs a pre-fitted Generalized Additive Model (GAM) to predict ancestral effective population size (*Ne*) based on synonymous substitution rates (*Ks*) and neutral mutation rates (*μ*). It addresses the inherent bias in divergence time estimation by incorporating coalescent theory, where species divergence time differs from gene divergence time by \(2*Ne*\). The package automatically corrects Ks values using multiple substitution models and population genetics principles.

**Key Features**:
- Corrects divergence time estimates using *Ks* variance
- Requires only *Ks* distributions and neutral mutation rates
- Computationally efficient workflow
- Robust to moderate demographic fluctuations

## Installation

To install the development version from GitHub:

```r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("limj0987/Tspecies")
```

## Quick Start

```
# Example dataset
Ks_values <- c(0.257, 0.202, 0.066, 0.197, 0.278, 0.089)
neutral_mutation_rate <- 8e-8  # 8×10⁻⁸ per generation

# Estimate divergence time
Tspecies(Ks = Ks_values, miu = neutral_mutation_rate)
```

## Usage

### Main Function

```
Tspecies(Ks, miu)
```

#### Arguments

- `Ks`: Numeric vector of synonymous substitution rates (d/L), where:
  - `d` = base differences between orthologous sequences
  - `L` = sequence length
- `miu`: Neutral mutation rate per generation (numeric)

#### Returns

- Numeric value representing estimated species divergence time in generations
- Console output includes:
  - Predicted effective population size (*Ne*)
  - Divergence time estimate
  - **0** if divergence time ≤0 (indicating no divergence)

## Workflow Diagram

```
[Ks Input] → [Jukes-Cantor Correction] → [Variance Calculation]
                   ↓
[μ Input] → [GAM Model Selection] → [Ne Prediction]
                   ↓
[Divergence Time Calculation: T = (Ks_mean/2μ) - 2Ne]
```

## Example Output

```
The predicted number of effective population size is around 45218 
Estimated divergence time of the species: 1.32e+06 (in generations)
```

## Model Specifications

- 23 pre-trained GAM models covering *μ* from 1e-10 to 1e-7
- Automatic model selection based on *μ* value
- Nonlinear relationship between *Ks* variance and *Ne*

## Example Used in Paper

Please download the folder "example" and run the Rscript "example.R".

## License  
This project is licensed under the **MIT License**.  
Copyright © 2025 Mijia 

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:  

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.  

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.  

For details, see [MIT License](https://opensource.org/licenses/MIT).  

# Declining food availability and habitat shifts drive community responses to marine hypoxia

The hypoxiaSSM folder contains all code necessary to run the model found in our paper.  Briefly, the code alters the baseline structure provided by the R package mizer (https://github.com/sizespectrum/mizer) by introducing oxygen dependence of several biological rates.  One can calibrate and project using any combination of oxygen dependence of each rate or set of rates and compare their performance on real data arising from the Baltic Sea.

## Abstract

Worsening marine hypoxia has had severe negative consequences for fish communities across the globe.  While individual- and population-level impacts of deoxygenation have been identified, it is unknown how they interact to drive changes in food webs.  To address this, we incorporated several major impacts of hypoxia, including declines in benthic resources, habitat shifts, increasing mortality, and changes to rates of feeding, assimilation, and reproductive efficiency, into an existing size spectrum food web modeling framework.  We used this structure to ask the following questions: which of these direct effects are most critical to capturing population and community dynamics in a representative hypoxic system, how do they interact to result in community responses to deoxygenation, and what are the potential consequences of these effects in the context of accelerating deoxygenation?  We tested the effect of different combinations of oxygen-dependent processes, driven by observed oxygen levels, on the food web model's ability to explain time series of observed somatic growth, diets, biomass, and fishery yields of commercially relevant species in the Baltic Sea.  Model results suggest that food availability is most critical to capturing observed dynamics.  There is also some evidence for oxygen-dependent habitat use and physological rates as drivers of observed dynamics.  Deoxygenation results in declining growth both of benthic and benthopelagic fish species, as the latter are unable to compensate for the loss of benthic resources by consuming more pelagic fish and resources.  Analysis of scenarios of ideal, declining, and degraded oxygen conditions show that deoxygenation results in a decline in somatic growth of predators, an altered habitat occupancy resulting in changing species interactions, and a shift in energy flow to benthopelagic predators from benthic to pelagic resources.  This may have important implications for management as oxygen declines or improves.
## Data

All data used to perform the analysis are contained within the Data folder in this repository.  Their origin is briefly described in the .keep text file in the Data folder.

## Requirements

The size spectrum models contained within the code are run in the R package mizer.  Description, installation, updates, and helpful tips are provided at the following website: https://sizespectrum.org/mizer/

## Usage

Follow these steps prior to first time usage to ensure the code runs properly on your machine:

1. Download the "hypoxiaSSM-main" zipped folder from Github in its entirety
2. Unzip and place "hypoxiaSSM-main" in your preferred directory
3. Navigate to the Code folder and open the file "depend.R" in R or RStudio
4. Confirm that all listed packages are installed on your machine; uncomment and run lines corresponding to packages you have not yet installed, then re-comment
5. Change the "mypath" variable to the file path which contains the hypoxiaSSM-main folder
6. Save and close "depend.R"

Upon all subsequent uses, we recommend running depend.R first, or pasting your own working directory over the code in each relevant file you use.

NOTE: Nearly every folder contains a text file called ".keep" which describes the files and subfolders within.  Read these for a more detailed understanding of the repository and how to use it.

## Contact us

If you have questions or concerns regarding this code, or would like help in re-formatting it for your own use, please do not hesitate to contant the corresponding author at:

elizabeth.duskey@slu.se

## License

Copyright 2022 Elizabeth Duskey, Michele Casini, Karing Limburg, and Anna Gårdmark

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
OLS (1): regress GeneticDistance Synthetic if AntibioticResistance<1 & GeneticDistance>=0 & FusionProtein<1 & YearPlasmidReceived>2005 & YearPlasmidReceived<2016 

OLS (2): regress GeneticDistance Synthetic PartLengthKB PartLengthKBXSynthetic if AntibioticResistance<1 & GeneticDistance>=0 & FusionProtein<1 & YearPlasmidReceived>2005 & YearPlasmidReceived<2016 

OLS (3): regress CrossKingdom Synthetic PartLengthKB PartLengthKBXSynthetic if AntibioticResistance<1 & GeneticDistance>=0 & FusionProtein<1 & YearPlasmidReceived>2005 & YearPlasmidReceived<2016 

Logit (4): regress CrossKingdom Synthetic PartLengthKB PartLengthKBXSynthetic if AntibioticResistance<1 & GeneticDistance>=0 & FusionProtein<1 & YearPlasmidReceived>2005 & YearPlasmidReceived<2016

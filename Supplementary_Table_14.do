OLS (5): regress GeneticDistance Synthetic if AntibioticResistance<1 & GeneticDistance>=0 & FusionProtein<1 & YearPlasmidReceived>2005 & YearPlasmidReceived<2016 & CRISPR<1

OLS (6): regress GeneticDistance Synthetic PartLengthKB PartLengthKBXSynthetic if AntibioticResistance<1 & GeneticDistance>=0 & FusionProtein<1 & YearPlasmidReceived>2005 & YearPlasmidReceived<2016 & CRISPR<1

OLS (7): regress CrossKingdom Synthetic PartLengthKB PartLengthKBXSynthetic if AntibioticResistance<1 & GeneticDistance>=0 & FusionProtein<1 & YearPlasmidReceived>2005 & YearPlasmidReceived<2016 & CRISPR<1

Logit (8): regress CrossKingdom Synthetic PartLengthKB PartLengthKBXSynthetic if AntibioticResistance<1 & GeneticDistance>=0 & FusionProtein<1 & YearPlasmidReceived>2005 & YearPlasmidReceived<2016 & CRISPR<1

# DRUIDom

**DRUIDom (DRUg Interacting Domain prediction): a computational method for predicting new drug/compound - target protein interactions for drug discovery and repurposing, via mapping ligands to structural domains**

Motivation: Predictive approaches such as virtual screening have been used in drug discovery with the objective of reducing experimental time, labor and costs. Current machine learning or network-based approaches have issues related to generalization, usability or model interpretability, especially due to complexity of the target protein structure and function, and the bias in system training datasets.

Proposed solution: Here we propose a new computational method “DRUIDom” to predict the interactions between drug candidate compounds and target proteins by utilizing the domain modularity of proteins, to overcome problems associated with current approaches. DRUIDom is composed of two methodological steps. First, small molecule compounds are statistically mapped to structural domains of their target proteins, where the binding site/region is expected to reside. This way, other proteins containing the mapped domain or domain combination become new candidate targets for the corresponding compounds. Next, a million-scale dataset of small molecule compounds, including the ones mapped to domains in the previous step, are clustered based on their molecular similarities, and their domain/protein associations has been propagated to other compounds within the same clusters.

Results: Experimentally verified bioactivity data points, obtained from public databases, are meticulously filtered to construct datasets of active/interacting and inactive/noninteracting compound-target pairs (¬¬~2.9 million data points), to be used as training data for calculating parameters of compound-domain mappings, which led to 27,032 high-confidence associations between 250 domains and 8,165 compounds, and a finalized output of 5 million new compound-protein interactions. DRUIDom is experimentally validated by the synthesis and the bioactivity analyses of compounds predicted to target LIM-Kinase proteins, which play critical roles in the regulation cell motility, cell cycle progression, and differentiation through actin filament dynamics. We showed that LIMK-inhibitor-2 and its novel derivatives significantly block the cancer cell migration through inhibition of LIMK phosphorylation and the downstream protein cofilin. One of the derivative compounds (LIMKi-2d) is identified as a promising candidate due to its action on resistant Mahlavu liver cancer cells. The results demonstrated that DRUIDom can be exploited to identify drug candidate compounds for known targets, and to predict new target proteins based on the defined compound-domain relationships.

## Drug/compound – target protein interaction prediction approach of DRUIDom:
<img src="https://user-images.githubusercontent.com/13165170/118031417-2c3eff80-b36f-11eb-9c58-d05040d2512c.png" width="893"> 

## Drug/compound - domain mapping procedure explained over toy examples:
<img src="https://user-images.githubusercontent.com/13165170/118361079-69370c00-b592-11eb-9dd2-2f666e23f452.png" width="900"> 

## Article Info & Citation

Protein Domain-Based Prediction of Drug/Compound–Target Interactions and Experimental Validation on LIM Kinases

Tunca Doğan1,2,3\*, Ece Akhan3,4, Marcus Baumann5, Altay Koyas3, Heval Atas3, Ian Baxendale6, Maria Martin7 and Rengul Cetin-Atalay3,8,\*<br>

1 Department of Computer Engineering, Hacettepe University, 06800 Ankara, Turkey<br>
2 Institute of Informatics, Hacettepe University, 06800 Ankara, Turkey<br>
3 CanSyL, Graduate School of Informatics, Middle East Technical University, 06800 Ankara, Turkey<br>
4 Center for Genomics and Rare Diseases & Biobank for Rare Diseases, Hacettepe University, 06230 Ankara, Turkey<br>
5 School of Chemistry, University College Dublin, D04 N2E2 Dublin, Ireland<br>
6 Department of Chemistry, University of Durham, DH1 3LE Durham, UK<br>
7 European Molecular Biology Laboratory, European Bioinformatics Institute (EMBL-EBI), Wellcome Genome Campus, CB10 1SD Hinxton, Cambridge, UK<br>
8 Section of Pulmonary and Critical Care Medicine, University of Chicago, Chicago IL, 60637, USA<br>
\*  Corresponding authors

If you find DRUIDom useful, please consider citing our article:

Doğan, T., Akhan, E., Baumann, M., Koyas, A., Atas, H., Baxendale, I., Martin, M., Cetin-Atalay, R. (2021). Protein domain-based prediction of drug/compound–target interactions and experimental validation on LIM kinases. *PLOS Computational Biology, 17*(11), e1009171.

DOI/URL: [10.1371/journal.pcbi.1009171](https://doi.org/10.1371/journal.pcbi.1009171)

## Programming Environment & Files

All analyses and computation in this study (available through: 'DRUIDom_Source_Code.m') were carried out in MATLAB by The MathWorks.

**Descriptions of files and folders:**

* **"DRUIDom_source_bioactivity_dataset_filtered.txt.zip":** Curated bioactivity dataset containing drug/compound - target protein interaction data points from experimental assays. This dataset has been constructed using data from ChEMBL and PubChem databases, and filtered to obtain binarized (i.e., active/interacting or inactive/noninteracting) drug/compound - target protein pairs.
* **"DRUIDom_domain-compound_mappings_finalized.txt":** Associations between single domains (from InterPro) and drugs/compounds obtained by applying the DRUIDom domain mapping procedure.
* **"DRUIDom_domainpair-compound_mappings_finalized.txt":** Associations between domain pairs (from InterPro) and drugs/compounds obtained by applying the DRUIDom domain pair mapping procedure.
* **"DRUIDom_DTI_predictions_novel_finalized.txt.zip":** Predicted drug/compound - target protein interactions that are novel (not found in the source dataset), produced by DRUIDom.
* **"DRUIDom_DTI_predictions_novel_pathwaybased.zip/z01 (2 files)":** Predicted drug/compound - target protein interactions (novel), grouped based on signaling and metabolic pathways of the corresponding target proteins.
* **"DRUIDom_Files":** A folder that contains all input and intermediary files generated during the application of DRUIDom pipeline. Large input and intermediary files (> 25 Mb) that could not be uploaded to the repository are available [here](https://drive.google.com/file/d/1aQ3t93LD6H3R4aKzVsMqSnzDTvaKsxsS/view?usp=sharing)

## Contact

For inquiries about DRUIDom please contact: tuncadogan@hacettepe.edu.tr, tuncadogan@gmail.com


## License

DRUIDom (DRUg Interacting Domains): a computational method for predicting new drug/compound - target protein interactions for drug discovery and repurposing, via mapping ligands to structural domains

Copyright (C) 2021 CanSyL

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

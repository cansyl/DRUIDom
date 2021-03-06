% DRUIDom v1.0
%
% DRUIDom (DRUg Interacting Domains): a computational method for
% predicting new drug/compound - target protein interactions for drug
% discovery and repurposing, via mapping ligands to structural domains
% where the binding site resides.
%
% Repository: https://github.com/cansyl/DRUIDom
%
% Contact: tuncadogan@hacettepe.edu.tr
%
% ------------------------------------------------------------------------
% 
% Article
%
% Protein Domain-Based Prediction of Compound–Target Interactions and 
% Experimental Validation on LIM Kinases.
%
% Authors: Tunca Dogan*, Ece Akhan, Marcus Baumann, Altay Koyas, Heval 
% Atas, Ian Baxendale, Maria Martin and Rengul Cetin-Atalay*
% * Corresponding authors
% 
% DOI: 
%
% ------------------------------------------------------------------------
%
% Copyright (C) 2021 CanSyL
% 
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License along
% with this program. If not, see http://www.gnu.org/licenses/.
%
% ------------------------------------------------------------------------
%
%
%
% ------------------------
% Source Code of the Study
% ------------------------





% ---------------------------
% Dataset preparation process
% ---------------------------


% Extracting active and inactive data points from the ExCAPE dataset and further filtering and organization operations:

% xz -d DRUIDom_Files/pubchem.chembl.dataset4publication_inchi_smiles_v2.tsv.xz
% cut -d$'\t' -f2,4,5,8,9,12 DRUIDom_Files/ExCAPE.pubchem.chembl.dataset4publication_inchi_smiles_v2.tsv > DRUIDom_Files/ExCAPE_all_datapoints.tsv
% cut -d$'\t' -f4,5 DRUIDom_Files/ExCAPE_all_datapoints.tsv > DRUIDom_Files/ExCAPE_all_datapoints_Taxonid_GeneName.tsv
% awk -F'\t' '$2~/N/' DRUIDom_Files/ExCAPE_all_datapoints.tsv > DRUIDom_Files/ExCAPE_inactive_datapoints.tsv
% awk -F'\t' '$2~/A/' DRUIDom_Files/ExCAPE_all_datapoints.tsv > DRUIDom_Files/ExCAPE_active_datapoints.tsv
% cut -d$'\t' -f4,5 DRUIDom_Files/ExCAPE_active_datapoints.tsv > DRUIDom_Files/ExCAPE_active_datapoints_Taxonid_GeneName.tsv
% cut -d$'\t' -f1,6 DRUIDom_Files/ExCAPE_active_datapoints.tsv > DRUIDom_Files/ExCAPE_active_datapoints_Compoundid_SMILES.tsv
% sed 1d DRUIDom_Files/ExCAPE_active_datapoints.tsv > DRUIDom_Files/ExCAPE_active_datapoints2.tsv (rename the file by removing 2)
% awk -F'\t' '{if($3=="" || $3 < 4.7) print $0}' DRUIDom_Files/ExCAPE_inactive_datapoints.tsv > DRUIDom_Files/ExCAPE_inactive_datapoints_filtered.tsv

% awk -F'\t' '{if(!$3) print $0}' DRUIDom_Files/ExCAPE_inactive_datapoints_filtered.tsv > DRUIDom_Files/ExCAPE_inactive_datapoints_filtered_noact.tsv
% awk -F'\t' '{if($3) print $0}' DRUIDom_Files/ExCAPE_inactive_datapoints_filtered.tsv > DRUIDom_Files/ExCAPE_inactive_datapoints_filtered_withact.tsv
% cut -d$'\t' -f4,5 DRUIDom_Files/ExCAPE_inactive_datapoints_filtered_noact.tsv > DRUIDom_Files/ExCAPE_inactive_datapoints_filtered_noact_Taxonid_GeneName.tsv
% cut -d$'\t' -f1,6 DRUIDom_Files/ExCAPE_inactive_datapoints_filtered_noact.tsv > DRUIDom_Files/ExCAPE_inactive_datapoints_filtered_noact_Compoundid_SMILES.tsv
% cut -d$'\t' -f4,5 DRUIDom_Files/ExCAPE_inactive_datapoints_filtered_withact.tsv > DRUIDom_Files/ExCAPE_inactive_datapoints_filtered_withact_Taxonid_GeneName.tsv
% cut -d$'\t' -f1,6 DRUIDom_Files/ExCAPE_inactive_datapoints_filtered_withact.tsv > DRUIDom_Files/ExCAPE_inactive_datapoints_filtered_withact_Compoundid_SMILES.tsv

% cut -d$'\t' -f1 DRUIDom_Files/DEEPScreen_ChEMBLv23_act_inact_filtered_data_original.txt > DRUIDom_Files/DEEPScreen_ChEMBLv23_act_inact_filtered_data_Targetid.txt
% cut -d$'\t' -f2 DRUIDom_Files/DEEPScreen_ChEMBLv23_act_inact_filtered_data_original.txt > DRUIDom_Files/DEEPScreen_ChEMBLv23_act_inact_filtered_data_Compoundid.txt

% (Filename starting with "DEEPScreen" contains drug/compound-target protein interaction/bioactivity data points obtained from the ChEMBL database by rigorous
% filtering operations; detailed information regarding the preparation of this dataset can be obtained from our previous publication at:
% https://pubs.rsc.org/en/content/articlehtml/2020/sc/c9sc03414e)



% Mapping UniProt accessions to targets in ExCAPE dataset:

UniProt_mapping=cell(0,3);
[UniProt_UniProtAcc,UniProt_TaxonID,UniProt_GeneName]=textread('DRUIDom_Files/Mapping_UniProtAcc_TaxonID_GeneName.tab', '%s %s %s', 'delimiter', '\t', 'headerlines', 1);
to=1;
for i=1:length(UniProt_UniProtAcc)
    x=strsplit(UniProt_GeneName{i,1},' ')';
    UniProt_mapping(to:(to+length(x)-1),1)=repmat(UniProt_UniProtAcc(i,1),length(x),1);
    UniProt_mapping(to:(to+length(x)-1),2)=repmat(UniProt_TaxonID(i,1),length(x),1);
    UniProt_mapping(to:(to+length(x)-1),3)=x;
    to=to+length(x);
end
UniProt_mapping=upper(UniProt_mapping);
save DRUIDom_Files/UniProt_mapping.mat UniProt_mapping
UniProt_mapping_23=strcat(UniProt_mapping(:,2), {' '}, UniProt_mapping(:,3));
save DRUIDom_Files/UniProt_mapping_23.mat UniProt_mapping_23

% (ExCAPE all datapoints)
[ExCAPE_all_datapoints_Taxonid,ExCAPE_all_datapoints_GeneName]=textread('DRUIDom_Files/ExCAPE_all_datapoints_Taxonid_GeneName.tsv', '%s %s', 'delimiter', '\t');
ExCAPE_all_datapoints_Taxonid(1,:)=[];
ExCAPE_all_datapoints_GeneName(1,:)=[];
ExCAPE_all_datapoints_23=strcat(ExCAPE_all_datapoints_Taxonid(:,1), {' '}, ExCAPE_all_datapoints_GeneName(:,1));
clear ExCAPE_all_datapoints_Taxonid ExCAPE_all_datapoints_GeneName
load DRUIDom_Files/UniProt_mapping_23.mat
load DRUIDom_Files/UniProt_mapping.mat
[Lia,Locb]=ismember(ExCAPE_all_datapoints_23,UniProt_mapping_23);
ExCAPE_all_datapoints_UniProtAcc=cell(length(ExCAPE_all_datapoints_23),1);
ExCAPE_all_datapoints_UniProtAcc(Lia==1,1)=UniProt_mapping(Locb(Locb>0),1);
% save DRUIDom_Files/ExCAPE_all_datapoints_UniProtAcc.mat ExCAPE_all_datapoints_UniProtAcc
% (not possible to save as mat file due to extremely large size of the variable)
t=ExCAPE_all_datapoints_UniProtAcc';
fid = fopen('DRUIDom_Files/ExCAPE_all_datapoints_UniProtAcc.txt','w');
fprintf(fid,'%s\n',t{:});
fclose(fid);

% (ExCAPE active datapoints)
[ExCAPE_active_datapoints_Taxonid,ExCAPE_active_datapoints_GeneName]=textread('DRUIDom_Files/ExCAPE_active_datapoints_Taxonid_GeneName.tsv', '%s %s', 'delimiter', '\t');
ExCAPE_active_datapoints_23=strcat(ExCAPE_active_datapoints_Taxonid(:,1), {' '}, ExCAPE_active_datapoints_GeneName(:,1));
[Lia,Locb]=ismember(ExCAPE_active_datapoints_23,UniProt_mapping_23);
ExCAPE_active_datapoints_UniProtAcc=cell(length(ExCAPE_active_datapoints_23),1);
ExCAPE_active_datapoints_UniProtAcc(Lia==1,1)=UniProt_mapping(Locb(Locb>0),1);
save DRUIDom_Files/ExCAPE_active_datapoints_UniProtAcc.mat ExCAPE_active_datapoints_UniProtAcc
t=ExCAPE_active_datapoints_UniProtAcc';
fid = fopen('DRUIDom_Files/ExCAPE_active_datapoints_UniProtAcc.txt','w');
fprintf(fid,'%s\n',t{:});
fclose(fid);
[ExCAPE_active_datapoints_Compoundid,ExCAPE_active_datapoints_SMILES]=textread('DRUIDom_Files/ExCAPE_active_datapoints_Compoundid_SMILES.tsv', '%s %s', 'delimiter', '\t');
del_ind=find(cellfun(@isempty,ExCAPE_active_datapoints_UniProtAcc));
ExCAPE_active_datapoints_UniProtAcc(del_ind,:)=[];
ExCAPE_active_datapoints_Compoundid(del_ind,:)=[];
ExCAPE_active_datapoints_SMILES(del_ind,:)=[];
save DRUIDom_Files/ExCAPE_active_datapoints_Org_Var.mat ExCAPE_active_datapoints_UniProtAcc ExCAPE_active_datapoints_Compoundid ExCAPE_active_datapoints_SMILES
ExCAPE_Org_act=[ExCAPE_active_datapoints_UniProtAcc ExCAPE_active_datapoints_Compoundid];
save DRUIDom_Files/ExCAPE_Org_act.mat ExCAPE_Org_act

% (ExCAPE inactive datapoints with activity measures)
[ExCAPE_inactive_datapoints_filtered_withact_Taxonid,ExCAPE_inactive_datapoints_filtered_withact_GeneName]=textread('DRUIDom_Files/ExCAPE_inactive_datapoints_filtered_withact_Taxonid_GeneName.tsv', '%s %s', 'delimiter', '\t');
ExCAPE_inactive_datapoints_filtered_withact_23=strcat(ExCAPE_inactive_datapoints_filtered_withact_Taxonid(:,1), {' '}, ExCAPE_inactive_datapoints_filtered_withact_GeneName(:,1));
[Lia,Locb]=ismember(ExCAPE_inactive_datapoints_filtered_withact_23,UniProt_mapping_23);
ExCAPE_inactive_datapoints_filtered_withact_UniProtAcc=cell(length(ExCAPE_inactive_datapoints_filtered_withact_23),1);
ExCAPE_inactive_datapoints_filtered_withact_UniProtAcc(Lia==1,1)=UniProt_mapping(Locb(Locb>0),1);
save DRUIDom_Files/ExCAPE_inactive_datapoints_filtered_withact_UniProtAcc.mat ExCAPE_inactive_datapoints_filtered_withact_UniProtAcc
t=ExCAPE_inactive_datapoints_filtered_withact_UniProtAcc';
fid = fopen('DRUIDom_Files/ExCAPE_inactive_datapoints_filtered_withact_UniProtAcc.txt','w');
fprintf(fid,'%s\n',t{:});
fclose(fid);
[ExCAPE_inactive_datapoints_filtered_withact_Compoundid,ExCAPE_inactive_datapoints_filtered_withact_SMILES]=textread('DRUIDom_Files/ExCAPE_inactive_datapoints_filtered_withact_Compoundid_SMILES.tsv', '%s %s', 'delimiter', '\t');
del_ind=find(cellfun(@isempty,ExCAPE_inactive_datapoints_filtered_withact_UniProtAcc));
ExCAPE_inactive_datapoints_filtered_withact_UniProtAcc(del_ind,:)=[];
ExCAPE_inactive_datapoints_filtered_withact_Compoundid(del_ind,:)=[];
ExCAPE_inactive_datapoints_filtered_withact_SMILES(del_ind,:)=[];
save DRUIDom_Files/ExCAPE_inactive_datapoints_filtered_withact_Org_Var.mat ExCAPE_inactive_datapoints_filtered_withact_UniProtAcc ExCAPE_inactive_datapoints_filtered_withact_Compoundid ExCAPE_inactive_datapoints_filtered_withact_SMILES
ExCAPE_Org_inact_withact=[ExCAPE_inactive_datapoints_filtered_withact_UniProtAcc ExCAPE_inactive_datapoints_filtered_withact_Compoundid];
save DRUIDom_Files/ExCAPE_Org_inact_withact.mat ExCAPE_Org_inact_withact

ExCAPE_Org_act_inact=[ExCAPE_Org_act;ExCAPE_Org_inact_withact];
save DRUIDom_Files/ExCAPE_Org_act_inact.mat ExCAPE_Org_act_inact
length(unique(ExCAPE_Org_act_inact(:,1)))
length(unique(ExCAPE_Org_act_inact(:,2)))

% (ExCAPE inactive datapoints without any activity measures at all)
[ExCAPE_inactive_datapoints_filtered_noact_Taxonid,ExCAPE_inactive_datapoints_filtered_noact_GeneName]=textread('DRUIDom_Files/ExCAPE_inactive_datapoints_filtered_noact_Taxonid_GeneName.tsv', '%s %s', 'delimiter', '\t');
ExCAPE_inactive_datapoints_filtered_noact_23=strcat(ExCAPE_inactive_datapoints_filtered_noact_Taxonid(:,1), {' '}, ExCAPE_inactive_datapoints_filtered_noact_GeneName(:,1));
[Lia,Locb]=ismember(ExCAPE_inactive_datapoints_filtered_noact_23,UniProt_mapping_23);
ExCAPE_inactive_datapoints_filtered_noact_UniProtAcc=cell(length(ExCAPE_inactive_datapoints_filtered_noact_23),1);
ExCAPE_inactive_datapoints_filtered_noact_UniProtAcc(Lia==1,1)=UniProt_mapping(Locb(Locb>0),1);
% save DRUIDom_Files/ExCAPE_inactive_datapoints_filtered_noact_UniProtAcc.mat ExCAPE_inactive_datapoints_filtered_noact_UniProtAcc
% (not possible to save as mat file due to extremely large size of the variable)
t=ExCAPE_inactive_datapoints_filtered_noact_UniProtAcc';
fid = fopen('DRUIDom_Files/ExCAPE_inactive_datapoints_filtered_noact_UniProtAcc.txt','w');
fprintf(fid,'%s\n',t{:});
fclose(fid);


% (Generating the ExCAPE_PubChem dataset)
[ExCAPE_Org_act_PubChem,idx]=sort(ExCAPE_Org_act(:,2));
ExCAPE_Org_act_PubChem=ExCAPE_Org_act_PubChem(1:201013,1);
ExCAPE_Org_act_PubChem_prot=ExCAPE_Org_act(idx,1);
ExCAPE_Org_act_PubChem_prot=ExCAPE_Org_act_PubChem_prot(1:201013,1);
ExCAPE_Org_act_PubChem=[ExCAPE_Org_act_PubChem_prot ExCAPE_Org_act_PubChem];
save DRUIDom_Files/ExCAPE_Org_act_PubChem.mat ExCAPE_Org_act_PubChem
[ExCAPE_Org_inact_withact_PubChem,idx]=sort(ExCAPE_Org_inact_withact(:,2));
ExCAPE_Org_inact_withact_PubChem=ExCAPE_Org_inact_withact_PubChem(1:99946,1);
ExCAPE_Org_inact_withact_PubChem_prot=ExCAPE_Org_inact_withact(idx,1);
ExCAPE_Org_inact_withact_PubChem_prot=ExCAPE_Org_inact_withact_PubChem_prot(1:99946,1);
ExCAPE_Org_inact_withact_PubChem=[ExCAPE_Org_inact_withact_PubChem_prot ExCAPE_Org_inact_withact_PubChem];
save DRUIDom_Files/ExCAPE_Org_inact_withact_PubChem.mat ExCAPE_Org_inact_withact_PubChem
ExCAPE_Org_act_inact_PubChem=[ExCAPE_Org_act_PubChem;ExCAPE_Org_inact_withact_PubChem];
save DRUIDom_Files/ExCAPE_Org_act_inact_PubChem.mat ExCAPE_Org_act_inact_PubChem
length(unique(ExCAPE_Org_act_inact_PubChem(:,1)))
length(unique(ExCAPE_Org_act_inact_PubChem(:,2)))



% Mapping UniProt accessions to targets in DEEPScreen ChEMBLv23 training dataset:

[DEEPScreen_ChEMBLv23_act_inact_Target_ChEMBLid]=textread('DRUIDom_Files/DEEPScreen_ChEMBLv23_act_inact_filtered_data_Targetid.txt', '%s');
DEEPScreen_ChEMBLv23_act_inact_Target_ChEMBLid=strrep(DEEPScreen_ChEMBLv23_act_inact_Target_ChEMBLid,'_inact','');
DEEPScreen_ChEMBLv23_act_inact_Target_ChEMBLid=strrep(DEEPScreen_ChEMBLv23_act_inact_Target_ChEMBLid,'_act','');
[DEEPScreen_ChEMBLv23_act_inact_Compound_ChEMBLid]=textread('DRUIDom_Files/DEEPScreen_ChEMBLv23_act_inact_filtered_data_Compoundid.txt', '%s', 'delimiter', '\n', 'bufsize', 500000);

[Mapping_UniProtAcc,Mapping_ChEMBLid]=textread('DRUIDom_Files/Mapping_UniProtAcc_ChEMBLid.tab', '%s %s', 'delimiter', '\t', 'headerlines', 1);
[Lia,~]=ismember(DEEPScreen_ChEMBLv23_act_inact_Target_ChEMBLid,Mapping_ChEMBLid);
del_ind=find(Lia==0);
DEEPScreen_ChEMBLv23_act_inact_Target_ChEMBLid(del_ind,:)=[];
DEEPScreen_ChEMBLv23_act_inact_Compound_ChEMBLid(del_ind,:)=[];
[~,Locb]=ismember(DEEPScreen_ChEMBLv23_act_inact_Target_ChEMBLid,Mapping_ChEMBLid);
DEEPScreen_ChEMBLv23_act_inact_Target_UniProtAcc=Mapping_UniProtAcc(Locb);

to=1;
for i=1:2:length(DEEPScreen_ChEMBLv23_act_inact_Compound_ChEMBLid)
    if isempty(DEEPScreen_ChEMBLv23_act_inact_Compound_ChEMBLid{i,1})==0
        x=strsplit(DEEPScreen_ChEMBLv23_act_inact_Compound_ChEMBLid{i,1},',')';
        DEEPScreen_ChEMBLv23_Org_act(to:(to+length(x)-1),1)=repmat(DEEPScreen_ChEMBLv23_act_inact_Target_UniProtAcc(i,1),length(x),1);
        DEEPScreen_ChEMBLv23_Org_act(to:(to+length(x)-1),2)=x;
        to=to+length(x);
    end
end
to=1;
for i=2:2:length(DEEPScreen_ChEMBLv23_act_inact_Compound_ChEMBLid)
    if isempty(DEEPScreen_ChEMBLv23_act_inact_Compound_ChEMBLid{i,1})==0
        x=strsplit(DEEPScreen_ChEMBLv23_act_inact_Compound_ChEMBLid{i,1},',')';
        DEEPScreen_ChEMBLv23_Org_inact(to:(to+length(x)-1),1)=repmat(DEEPScreen_ChEMBLv23_act_inact_Target_UniProtAcc(i,1),length(x),1);
        DEEPScreen_ChEMBLv23_Org_inact(to:(to+length(x)-1),2)=x;
        to=to+length(x);
    end
end
save DRUIDom_Files/DEEPScreen_ChEMBLv23_Org_act_raw.mat DEEPScreen_ChEMBLv23_Org_act
save DRUIDom_Files/DEEPScreen_ChEMBLv23_Org_inact_raw.mat DEEPScreen_ChEMBLv23_Org_inact

[ChEMBLv23_Compoundid,ChEMBLv23_SMILES]=textread('DRUIDom_Files/ChEMBLv23_compound_SMILES.txt', '%s %s', 'delimiter', '\t', 'headerlines', 1);
[Lia,Locb]=ismember(DEEPScreen_ChEMBLv23_Org_act(:,2),ChEMBLv23_Compoundid);
DEEPScreen_ChEMBLv23_Org_act_SMILES=cell(length(DEEPScreen_ChEMBLv23_Org_act),1);
DEEPScreen_ChEMBLv23_Org_act_SMILES(Lia==1,1)=ChEMBLv23_SMILES(Locb(Locb>0),1);
del_ind=find(Lia==0);
DEEPScreen_ChEMBLv23_Org_act(del_ind,:)=[];
DEEPScreen_ChEMBLv23_Org_act_SMILES(del_ind,:)=[];
save DRUIDom_Files/DEEPScreen_ChEMBLv23_Org_act_var.mat DEEPScreen_ChEMBLv23_Org_act DEEPScreen_ChEMBLv23_Org_act_SMILES

[Lia,Locb]=ismember(DEEPScreen_ChEMBLv23_Org_inact(:,2),ChEMBLv23_Compoundid);
DEEPScreen_ChEMBLv23_Org_inact_SMILES=cell(length(DEEPScreen_ChEMBLv23_Org_inact),1);
DEEPScreen_ChEMBLv23_Org_inact_SMILES(Lia==1,1)=ChEMBLv23_SMILES(Locb(Locb>0),1);
del_ind=find(Lia==0);
DEEPScreen_ChEMBLv23_Org_inact(del_ind,:)=[];
DEEPScreen_ChEMBLv23_Org_inact_SMILES(del_ind,:)=[];
save DRUIDom_Files/DEEPScreen_ChEMBLv23_Org_inact_var.mat DEEPScreen_ChEMBLv23_Org_inact DEEPScreen_ChEMBLv23_Org_inact_SMILES

DEEPScreen_ChEMBLv23_Org_act_inact=[DEEPScreen_ChEMBLv23_Org_act;DEEPScreen_ChEMBLv23_Org_inact];
save DRUIDom_Files/DEEPScreen_ChEMBLv23_Org_act_inact.mat DEEPScreen_ChEMBLv23_Org_act_inact
length(unique(DEEPScreen_ChEMBLv23_Org_act_inact(:,1)))
length(unique(DEEPScreen_ChEMBLv23_Org_act_inact(:,2)))



% Comparison and the merge between ExCAPE and DEEPScreen training datasets: 

% (actives dataset)
DEEPScreen_Org_act_51=strcat(DEEPScreen_ChEMBLv23_Org_act(:,1), {' '}, DEEPScreen_ChEMBLv23_Org_act(:,2));
ExCAPE_Org_act_51=strcat(ExCAPE_Org_act(:,1), {' '}, ExCAPE_Org_act(:,2));
[Lia,~]=ismember(ExCAPE_Org_act_51,DEEPScreen_Org_act_51);
ExCAPE_Org_act_onlyExCAPE=ExCAPE_Org_act(Lia==0,:);
ExCAPE_Org_act_onlyExCAPE_SMILES=ExCAPE_active_datapoints_SMILES(Lia==0,:);
save DRUIDom_Files/ExCAPE_Org_act_onlyExCAPE_var.mat ExCAPE_Org_act_onlyExCAPE ExCAPE_Org_act_onlyExCAPE_SMILES
Combined_Org_act=[DEEPScreen_ChEMBLv23_Org_act;ExCAPE_Org_act_onlyExCAPE];
Combined_Org_act_SMILES=[DEEPScreen_ChEMBLv23_Org_act_SMILES;ExCAPE_Org_act_onlyExCAPE_SMILES];
save DRUIDom_Files/Combined_Org_act_var.mat Combined_Org_act Combined_Org_act_SMILES
t=Combined_Org_act_SMILES';
fid = fopen('DRUIDom_Files/Combined_Org_act_SMILES.txt','w');
fprintf(fid,'%s\n',t{:});
fclose(fid);
Combined_Org_act_SMILES_id=[Combined_Org_act_SMILES Combined_Org_act(:,2)];
t=Combined_Org_act_SMILES_id';
fid = fopen('DRUIDom_Files/Combined_Org_act_SMILES_id.smi','w');
fprintf(fid,'%s\t%s\n',t{:});
fclose(fid);

% (inactives dataset)
DEEPScreen_Org_inact_51=strcat(DEEPScreen_ChEMBLv23_Org_inact(:,1), {' '}, DEEPScreen_ChEMBLv23_Org_inact(:,2));
ExCAPE_Org_inact_withact_51=strcat(ExCAPE_Org_inact_withact(:,1), {' '}, ExCAPE_Org_inact_withact(:,2));
[Lia,~]=ismember(ExCAPE_Org_inact_withact_51,DEEPScreen_Org_inact_51);
ExCAPE_Org_inact_withact_onlyExCAPE=ExCAPE_Org_inact_withact(Lia==0,:);
ExCAPE_Org_inact_withact_onlyExCAPE_SMILES=ExCAPE_inactive_datapoints_filtered_withact_SMILES(Lia==0,:);
save DRUIDom_Files/ExCAPE_Org_inact_withact_onlyExCAPE_var.mat ExCAPE_Org_inact_withact_onlyExCAPE ExCAPE_Org_inact_withact_onlyExCAPE_SMILES
Combined_Org_inact_withact=[DEEPScreen_ChEMBLv23_Org_inact;ExCAPE_Org_inact_withact_onlyExCAPE];
Combined_Org_inact_withact_SMILES=[DEEPScreen_ChEMBLv23_Org_inact_SMILES;ExCAPE_Org_inact_withact_onlyExCAPE_SMILES];
save DRUIDom_Files/Combined_Org_inact_withact_var.mat Combined_Org_inact_withact Combined_Org_inact_withact_SMILES
t=Combined_Org_inact_withact_SMILES';
fid = fopen('DRUIDom_Files/Combined_Org_inact_withact_SMILES.txt','w');
fprintf(fid,'%s\n',t{:});
fclose(fid);
Combined_Org_inact_withact_SMILES_id=[Combined_Org_inact_withact_SMILES Combined_Org_inact_withact(:,2)];
t=Combined_Org_inact_withact_SMILES_id';
fid = fopen('DRUIDom_Files/Combined_Org_inact_withact_SMILES_id.smi','w');
fprintf(fid,'%s\t%s\n',t{:});
fclose(fid);

Combined_Org_All_id=[Combined_Org_act(:,2);Combined_Org_inact_withact(:,2)];
[Combined_Org_All_id_unique,I,J] = unique(Combined_Org_All_id);
Combined_Org_All_SMILES=[Combined_Org_act_SMILES;Combined_Org_inact_withact_SMILES];
Combined_Org_All_SMILES_unique=Combined_Org_All_SMILES(I);
Combined_Org_All_SMILES_id_unique=[Combined_Org_All_SMILES_unique Combined_Org_All_id_unique];
t=Combined_Org_All_SMILES_id_unique';
fid = fopen('DRUIDom_Files/Combined_Org_All_SMILES_id_unique.smi','w');
fprintf(fid,'%s\t%s\n',t{:});
fclose(fid);

Combined_Org_act_51=strcat(Combined_Org_act(:,1), {' '}, Combined_Org_act(:,2));
Combined_Org_inact_withact_51=strcat(Combined_Org_inact_withact(:,1), {' '}, Combined_Org_inact_withact(:,2));
[Lia,~]=ismember(Combined_Org_act_51,Combined_Org_inact_withact_51);
sum(Lia)
% (there are very few corresponding data points between actives and inactives: 1574)

Combined_Org_act_inact=[Combined_Org_act;Combined_Org_inact_withact];
save DRUIDom_Files/Combined_Org_act_inact.mat Combined_Org_act_inact
length(unique(Combined_Org_act_inact(:,1)))
length(unique(Combined_Org_act_inact(:,2)))
Target_UniProtacc_Combined_Org_act_inact=unique(Combined_Org_act_inact(:,1));
save DRUIDom_Files/Target_UniProtacc_Combined_Org_act_inact.mat Target_UniProtacc_Combined_Org_act_inact



% Histogram of data point distributions for all targets:

length(unique(Combined_Org_act(:,1)))
length(unique(Combined_Org_act(:,2)))
[~,Locb]=ismember(Combined_Org_act(:,1),unique(Combined_Org_act(:,1)));
Hist_combined_act=hist(Locb,0.5:1:(max(Locb)+0.5));
figure;bar(log10(sort(Hist_combined_act,'descend')))
[~,Locb]=ismember(Combined_Org_act(:,2),unique(Combined_Org_act(:,2)));
Hist_combined_act_comp=sort(hist(Locb,0.5:1:(max(Locb)+0.5)),'descend')';
length(find(Hist_combined_act_comp>=5))

length(unique(Combined_Org_inact_withact(:,1)))
length(unique(Combined_Org_inact_withact(:,2)))
[~,Locb]=ismember(Combined_Org_inact_withact(:,1),unique(Combined_Org_inact_withact(:,1)));
Hist_combined_inact_withact=hist(Locb,0.5:1:(max(Locb)+0.5));
figure;bar(log10(sort(Hist_combined_inact_withact,'descend')))
[~,Locb]=ismember(Combined_Org_inact_withact(:,2),unique(Combined_Org_inact_withact(:,2)));
Hist_combined_inact_withact_comp=sort(hist(Locb,0.5:1:(max(Locb)+0.5)),'descend')';
length(find(Hist_combined_inact_withact_comp>=5))



% Generation of the Morgan Fingerprints for all compounds and pairwise tanimoto similarity search using Chemfp:

% conda create -n py2condaenv python=2
% conda activate py2-rdkit-env
% python .../chemfp-1.5/setup.py install --without-openmp
% cd .../DRUIDom
% rdkit2fps --morgan DRUIDom_Files/Combined_Org_All_SMILES_id_unique.smi -o DRUIDom_Files/Combined_Org_All_ECFP4_id_unique.fps
% split -l 206705 DRUIDom_Files/Combined_Org_All_ECFP4_id_unique.fps
% simsearch  --times --threshold 0.5 -q xaa xaa -o DRUIDom_Files/xaa_xaa_Simsearch.txt
% simsearch  --times --threshold 0.5 -q xaa xab -o DRUIDom_Files/xaa_xab_Simsearch.txt
% simsearch  --times --threshold 0.5 -q xaa xac -o DRUIDom_Files/xaa_xac_Simsearch.txt
% simsearch  --times --threshold 0.5 -q xaa xad -o DRUIDom_Files/xaa_xad_Simsearch.txt
% simsearch  --times --threshold 0.5 -q xaa xae -o DRUIDom_Files/xaa_xae_Simsearch.txt
% simsearch  --times --threshold 0.5 -q xab xab -o DRUIDom_Files/xab_xab_Simsearch.txt
% simsearch  --times --threshold 0.5 -q xab xac -o DRUIDom_Files/xab_xac_Simsearch.txt
% simsearch  --times --threshold 0.5 -q xab xad -o DRUIDom_Files/xab_xad_Simsearch.txt
% simsearch  --times --threshold 0.5 -q xab xae -o DRUIDom_Files/xab_xae_Simsearch.txt
% simsearch  --times --threshold 0.5 -q xac xac -o DRUIDom_Files/xac_xac_Simsearch.txt
% simsearch  --times --threshold 0.5 -q xac xad -o DRUIDom_Files/xac_xad_Simsearch.txt
% simsearch  --times --threshold 0.5 -q xac xae -o DRUIDom_Files/xac_xae_Simsearch.txt
% simsearch  --times --threshold 0.5 -q xad xad -o DRUIDom_Files/xad_xad_Simsearch.txt
% simsearch  --times --threshold 0.5 -q xad xae -o DRUIDom_Files/xad_xae_Simsearch.txt
% simsearch  --times --threshold 0.5 -q xae xae -o DRUIDom_Files/xae_xae_Simsearch.txt
% simsearch  --times --threshold 0.5 --NxN DRUIDom_Files/Combined_Org_All_ECFP4_id_unique.fps -o DRUIDom_Files/Combined_Org_All_Simsearch.txt



% Organizing the compound similarity files (repeat below code for all 15 parts from xaa to xae):

[Part_Simsearch]=textread('DRUIDom_Files/xaa_xaa_Simsearch.txt', '%s', 'delimiter', '\n', 'bufsize', 500000);
Part_SimS_Org=cell(10000000,3);
to=1;
for i=1:length(Part_Simsearch)
    disp(['Line number: ' num2str(i), ' / ' num2str(length(Part_Simsearch))])
    x=strsplit(Part_Simsearch{i,1},'\t')';
    if str2num(x{1,1})~=0
        Part_SimS_Org(to:(to+((length(x)-2)/2)-1),1)=repmat(x(2,1),(length(x)-2)/2,1);
        Part_SimS_Org(to:(to+((length(x)-2)/2)-1),2)=x(3:2:(end-1),1);
        Part_SimS_Org(to:(to+((length(x)-2)/2)-1),3)=x(4:2:(end),1);
        to=to+((length(x)-2)/2);
    end
end
del_ind=cellfun(@isempty,Part_SimS_Org(:,1));
Part_SimS_Org(del_ind,:)=[];
self_ind=cellfun(@isequal,Part_SimS_Org(:,1),Part_SimS_Org(:,2));
Part_SimS_Org(self_ind==1,:)=[];
t=Part_SimS_Org';
fid=fopen('DRUIDom_Files/xaa_xaa_SimS_Org.txt','w');
fprintf(fid,'%s\t%s\t%s\n',t{:});
fclose(fid);

% cat DRUIDom_Files/xaa_xaa_SimS_Org.txt DRUIDom_Files/xaa_xab_SimS_Org.txt DRUIDom_Files/xaa_xac_SimS_Org.txt DRUIDom_Files/xaa_xad_SimS_Org.txt DRUIDom_Files/xaa_xae_SimS_Org.txt DRUIDom_Files/xab_xab_SimS_Org.txt DRUIDom_Files/xab_xac_SimS_Org.txt DRUIDom_Files/xab_xad_SimS_Org.txt DRUIDom_Files/xab_xae_SimS_Org.txt DRUIDom_Files/xac_xac_SimS_Org.txt DRUIDom_Files/xac_xad_SimS_Org.txt DRUIDom_Files/xac_xae_SimS_Org.txt DRUIDom_Files/xad_xad_SimS_Org.txt DRUIDom_Files/xad_xae_SimS_Org.txt DRUIDom_Files/xae_xae_SimS_Org.txt > DRUIDom_Files/Combined_SimS_Org.txt



% Generating the combined protein groups for each ligand cluster using the selected ligand similarity threshold:

[Combined_Org_Comp1,Combined_Org_Comp2,Combined_Org_Sim]=textread('DRUIDom_Files/Combined_SimS_Org.txt', '%s %s %d', 'delimiter', '\t');
load DRUIDom_Files/Combined_Org_act_var.mat
load DRUIDom_Files/Combined_Org_inact_withact_var.mat
Combined_Org_id_unique=unique([Combined_Org_act(:,2);Combined_Org_inact_withact(:,2)]);
save DRUIDom_Files/Combined_Org_id_unique.mat Combined_Org_id_unique

thr=0.7;

ind_thr=find(Combined_Org_Sim>=thr);
Combined_Org_Sim=Combined_Org_Sim(ind_thr);
Combined_Org_Comp1=Combined_Org_Comp1(ind_thr);
Combined_Org_Comp2=Combined_Org_Comp2(ind_thr);
save DRUIDom_Files/Combined_Org_Comp1_07.mat Combined_Org_Comp1
save DRUIDom_Files/Combined_Org_Comp2_07.mat Combined_Org_Comp2
save DRUIDom_Files/Combined_Org_Sim_07.mat Combined_Org_Sim

[Lia1,Locb1]=ismember(Combined_Org_Comp1,Combined_Org_id_unique);
[Lia2,Locb2]=ismember(Combined_Org_Comp2,Combined_Org_id_unique);
save DRUIDom_Files/Locations_07.mat Locb1 Locb2
[Lia_act,Locb_act]=ismember(Combined_Org_act(:,2),Combined_Org_id_unique);
[Lia_inact,Locb_inact]=ismember(Combined_Org_inact_withact(:,2),Combined_Org_id_unique);
save DRUIDom_Files/Locations_act_inact.mat Locb_act Locb_inact
Combined_Org_id_unique_ind=(1:length(Combined_Org_id_unique))';

parpool(4)
Combined_Org_unique_act_Prot=cell(length(Combined_Org_id_unique),1);
Combined_Org_unique_inact_withact_Prot=cell(length(Combined_Org_id_unique),1);
parfor i=1:length(Combined_Org_id_unique)
    disp(['Compound number: ' num2str(i), ' / ' num2str(length(Combined_Org_id_unique))])
    Combined_Org_thrSim_temp_ind=unique([i;Locb2(Locb1==i);Locb1(Locb2==i)]);
    Combined_Org_unique_act_Prot{i,1}=unique(Combined_Org_act(ismember(Locb_act,Combined_Org_thrSim_temp_ind)==1,1));
    Combined_Org_unique_inact_withact_Prot{i,1}=unique(Combined_Org_inact_withact(ismember(Locb_inact,Combined_Org_thrSim_temp_ind)==1,1));
end
save DRUIDom_Files/Combined_Org_unique_act_Prot_07.mat Combined_Org_unique_act_Prot
save DRUIDom_Files/Combined_Org_unique_inact_withact_Prot_07.mat Combined_Org_unique_inact_withact_Prot



% Analyzing and thresholding the protein arrays:

Combined_Org_unique_Samplenum_act_inact_sup=zeros(length(Combined_Org_id_unique),3);
for i=1:length(Combined_Org_id_unique)
    disp(['Compound number: ' num2str(i), ' / ' num2str(length(Combined_Org_id_unique))])
    Combined_Org_unique_Samplenum_act_inact_sup(i,1)=length(Combined_Org_unique_act_Prot{i,1});
    Combined_Org_unique_Samplenum_act_inact_sup(i,2)=length(Combined_Org_unique_inact_withact_Prot{i,1});
    Combined_Org_unique_Samplenum_act_inact_sup(i,3)=(Combined_Org_unique_Samplenum_act_inact_sup(i,1)+Combined_Org_unique_Samplenum_act_inact_sup(i,2))/2;
end
save DRUIDom_Files/Combined_Org_unique_Samplenum_act_inact_sup.mat Combined_Org_unique_Samplenum_act_inact_sup

length(find(Combined_Org_unique_Samplenum_act_inact_sup(:,1)>=3))
length(find(Combined_Org_unique_Samplenum_act_inact_sup(:,2)>=3))
length(find(Combined_Org_unique_Samplenum_act_inact_sup(:,1)>=3 & Combined_Org_unique_Samplenum_act_inact_sup(:,2)>=3))
length(find(Combined_Org_unique_Samplenum_act_inact_sup(:,1)>=5))
length(find(Combined_Org_unique_Samplenum_act_inact_sup(:,2)>=5))
length(find(Combined_Org_unique_Samplenum_act_inact_sup(:,1)>=5 & Combined_Org_unique_Samplenum_act_inact_sup(:,2)>=5))
length(find(Combined_Org_unique_Samplenum_act_inact_sup(:,1)>=10))
length(find(Combined_Org_unique_Samplenum_act_inact_sup(:,2)>=10))
length(find(Combined_Org_unique_Samplenum_act_inact_sup(:,1)>=10 & Combined_Org_unique_Samplenum_act_inact_sup(:,2)>=10))

% (histograms of active and inactive distribution in the Prot5 dataset)
figure;
histogram(sort(2*Combined_Org_unique_Samplenum_act_inact_sup(:,3),'descend'),-0.5:1:(max(2*Combined_Org_unique_Samplenum_act_inact_sup(:,3))+0.5),'FaceAlpha',1)
set(gca,'YScale','log')
axis([0.5 350 0 1000000])
figure;
histogram(sort(Combined_Org_unique_Samplenum_act_inact_sup(:,1),'descend'),-0.5:1:(max(Combined_Org_unique_Samplenum_act_inact_sup(:,1))+0.5),'FaceColor',[0 0.7 0],'FaceAlpha',1)
set(gca,'YScale','log')
axis([0.5 350 0 1000000])
figure;
histogram(sort(Combined_Org_unique_Samplenum_act_inact_sup(:,2),'descend'),-0.5:1:(max(Combined_Org_unique_Samplenum_act_inact_sup(:,2))+0.5),'FaceColor',[0.8 0 0],'FaceAlpha',1)
set(gca,'YScale','log')
axis([0.5 350 0 1000000])


thr_act_samp=5;
thr_inact_samp=5;

Combined_Org_unique_Samplenum_act_inact_sup_5=Combined_Org_unique_Samplenum_act_inact_sup(Combined_Org_unique_Samplenum_act_inact_sup(:,1)>=thr_act_samp & Combined_Org_unique_Samplenum_act_inact_sup(:,2)>=thr_inact_samp,:);
Combined_Org_unique_act_Prot_5=Combined_Org_unique_act_Prot(Combined_Org_unique_Samplenum_act_inact_sup(:,1)>=thr_act_samp & Combined_Org_unique_Samplenum_act_inact_sup(:,2)>=thr_inact_samp,:);
Combined_Org_unique_inact_withact_Prot_5=Combined_Org_unique_inact_withact_Prot(Combined_Org_unique_Samplenum_act_inact_sup(:,1)>=thr_act_samp & Combined_Org_unique_Samplenum_act_inact_sup(:,2)>=thr_inact_samp,:);
Combined_Org_id_unique_Prot_5=Combined_Org_id_unique(Combined_Org_unique_Samplenum_act_inact_sup(:,1)>=thr_act_samp & Combined_Org_unique_Samplenum_act_inact_sup(:,2)>=thr_inact_samp,:);
save DRUIDom_Files/Combined_Org_unique_act_inact_Prot_5_var.mat Combined_Org_unique_Samplenum_act_inact_sup_5 Combined_Org_unique_act_Prot_5 Combined_Org_unique_inact_withact_Prot_5 Combined_Org_id_unique_Prot_5



% InterPro domain array generation for active and inactive sets:

% (saving a unique list of target protein accessions from the combined activity set)
Combined_Org_act_inact_Prot_unique=unique([Combined_Org_act(:,1);Combined_Org_inact_withact(:,1)]);
save DRUIDom_Files/Combined_Org_act_inact_Prot_unique.mat Combined_Org_act_inact_Prot_unique
% (UniProt id retrieval is used to obtain InterPro hits for these 3644 target proteins)


% (Organizing InterPro hits of our target proteins)
[InterPro_entriesmatched_onlytarget_UniProtAcc,InterPro_entriesmatched_onlytarget_IPRid]=textread('DRUIDom_Files/InterPro_UniProtAcc_IPRid_entriesmatched_onlytargets.tab', '%s %s', 'delimiter', '\t', 'headerlines', 1);
InterPro_entriesmatched_onlytarget_Org=cell(0,2);
to=1;
for i=1:length(InterPro_entriesmatched_onlytarget_UniProtAcc)
    disp(['Line number: ' num2str(i), ' / ' num2str(length(InterPro_entriesmatched_onlytarget_UniProtAcc))])
    if isempty(InterPro_entriesmatched_onlytarget_IPRid{i,1})==0
        x=strsplit(InterPro_entriesmatched_onlytarget_IPRid{i,1},';')';
        x(end,:)=[];
        InterPro_entriesmatched_onlytarget_Org(to:(to+(length(x))-1),1)=repmat(InterPro_entriesmatched_onlytarget_UniProtAcc(i,1),length(x),1);
        InterPro_entriesmatched_onlytarget_Org(to:(to+(length(x))-1),2)=x;
        to=to+(length(x));
    end
end
save DRUIDom_Files/InterPro_entriesmatched_onlytarget_Org.mat InterPro_entriesmatched_onlytarget_Org

% (Organizing InterPro hits of all ChEMBL targets)
[InterPro_entriesmatched_onlyChEMBL_UniProtAcc,InterPro_entriesmatched_onlyChEMBL_IPRid]=textread('DRUIDom_Files/InterPro_UniProtAcc_IPRid_entriesmatched_onlyChEMBL.tab', '%s %s', 'delimiter', '\t', 'headerlines', 1);
InterPro_entriesmatched_onlyChEMBL_Org=cell(0,2);
to=1;
for i=1:length(InterPro_entriesmatched_onlyChEMBL_UniProtAcc)
    disp(['Line number: ' num2str(i), ' / ' num2str(length(InterPro_entriesmatched_onlyChEMBL_UniProtAcc))])
    if isempty(InterPro_entriesmatched_onlyChEMBL_IPRid{i,1})==0
        x=strsplit(InterPro_entriesmatched_onlyChEMBL_IPRid{i,1},';')';
        x(end,:)=[];
        InterPro_entriesmatched_onlyChEMBL_Org(to:(to+(length(x))-1),1)=repmat(InterPro_entriesmatched_onlyChEMBL_UniProtAcc(i,1),length(x),1);
        InterPro_entriesmatched_onlyChEMBL_Org(to:(to+(length(x))-1),2)=x;
        to=to+(length(x));
    end
end
save DRUIDom_Files/InterPro_entriesmatched_onlyChEMBL_Org.mat InterPro_entriesmatched_onlyChEMBL_Org

% (Organizing InterPro hits of all human proteins)
[InterPro_entriesmatched_allhuman_UniProtAcc,InterPro_entriesmatched_allhuman_IPRid]=textread('DRUIDom_Files/InterPro_UniProtAcc_IPRid_entriesmatched_allhuman.tab', '%s %s', 'delimiter', '\t', 'headerlines', 1);
InterPro_entriesmatched_allhuman_Org=cell(0,2);
to=1;
for i=1:length(InterPro_entriesmatched_allhuman_UniProtAcc)
    disp(['Line number: ' num2str(i), ' / ' num2str(length(InterPro_entriesmatched_allhuman_UniProtAcc))])
    if isempty(InterPro_entriesmatched_allhuman_IPRid{i,1})==0
        x=strsplit(InterPro_entriesmatched_allhuman_IPRid{i,1},';')';
        x(end,:)=[];
        InterPro_entriesmatched_allhuman_Org(to:(to+(length(x))-1),1)=repmat(InterPro_entriesmatched_allhuman_UniProtAcc(i,1),length(x),1);
        InterPro_entriesmatched_allhuman_Org(to:(to+(length(x))-1),2)=x;
        to=to+(length(x));
    end
end
save DRUIDom_Files/InterPro_entriesmatched_allhuman_Org.mat InterPro_entriesmatched_allhuman_Org


% (Selecting only domain type hits from the organized protein InterPro hit files)
[InterPro_v72_Domain_id,~,~]=textread('DRUIDom_Files/InterPro_v72_Domain_entrylist.txt', '%s %s %s', 'delimiter', '\t');
InterPro_v72_Domain_id=unique(InterPro_v72_Domain_id);

[Lia,~]=ismember(InterPro_entriesmatched_onlytarget_Org(:,2),InterPro_v72_Domain_id);
InterPro_entriesmatched_onlytarget_Org_onlydom=InterPro_entriesmatched_onlytarget_Org(Lia==1,:);
save DRUIDom_Files/InterPro_entriesmatched_onlytarget_Org_onlydom.mat InterPro_entriesmatched_onlytarget_Org_onlydom

[Lia,~]=ismember(InterPro_entriesmatched_onlyChEMBL_Org(:,2),InterPro_v72_Domain_id);
InterPro_entriesmatched_onlyChEMBL_Org_onlydom=InterPro_entriesmatched_onlyChEMBL_Org(Lia==1,:);
save DRUIDom_Files/InterPro_entriesmatched_onlyChEMBL_Org_onlydom.mat InterPro_entriesmatched_onlyChEMBL_Org_onlydom

[Lia,~]=ismember(InterPro_entriesmatched_allhuman_Org(:,2),InterPro_v72_Domain_id);
InterPro_entriesmatched_allhuman_Org_onlydom=InterPro_entriesmatched_allhuman_Org(Lia==1,:);
save DRUIDom_Files/InterPro_entriesmatched_allhuman_Org_onlydom.mat InterPro_entriesmatched_allhuman_Org_onlydom


% (organizing InterPro term relation file to determine domain relations)
[InterPro_v72_ParentChildTree,~]=textread('DRUIDom_Files/InterPro_v72_ParentChildTreeFile.txt', '%s %s', 'delimiter', '::');
InterPro_v72_ParentChildTree=InterPro_v72_ParentChildTree(1:2:end,1);

% (removal of non-domain IRP entries from the hierarchy file)
for i=1:length(InterPro_v72_ParentChildTree)
    InterPro_v72_terms_from_ParentChildTree(i,1)=erase(InterPro_v72_ParentChildTree(i,1),"-");
end
[Lia,~]=ismember(InterPro_v72_terms_from_ParentChildTree,InterPro_v72_Domain_id);
InterPro_v72_ParentChildTree_Dom=InterPro_v72_ParentChildTree(Lia==1,:);
% (addition of domain entries to the table that have no hierarchical relationship with other entries)
InterPro_v72_terms_without_hierarc=setdiff(InterPro_v72_Domain_id,InterPro_v72_terms_from_ParentChildTree);
InterPro_v72_ParentChildTree_Dom=[InterPro_v72_ParentChildTree_Dom;InterPro_v72_terms_without_hierarc];
save DRUIDom_Files/InterPro_v72_ParentChildTree_var.mat InterPro_v72_ParentChildTree InterPro_v72_ParentChildTree_Dom

InterPro_v72_Dom_groups=cell(length(InterPro_v72_ParentChildTree_Dom),5);
for i=1:length(InterPro_v72_ParentChildTree_Dom)
    if isequal(InterPro_v72_ParentChildTree_Dom{i,1}(1,1),'-')==0
        Level1=erase(InterPro_v72_ParentChildTree_Dom(i,1),"-");
        InterPro_v72_Dom_groups(i,1)=erase(InterPro_v72_ParentChildTree_Dom(i,1),"-");
    elseif isequal(InterPro_v72_ParentChildTree_Dom{i,1}(1,1:3),'--I')==1
        Level2=erase(InterPro_v72_ParentChildTree_Dom(i,1),"-");
        InterPro_v72_Dom_groups(i,1)=Level1;
        InterPro_v72_Dom_groups(i,2)=Level2;
    elseif isequal(InterPro_v72_ParentChildTree_Dom{i,1}(1,1:5),'----I')==1
        Level3=erase(InterPro_v72_ParentChildTree_Dom(i,1),"-");
        InterPro_v72_Dom_groups(i,1)=Level1;
        InterPro_v72_Dom_groups(i,2)=Level2;
        InterPro_v72_Dom_groups(i,3)=Level3;
    elseif isequal(InterPro_v72_ParentChildTree_Dom{i,1}(1,1:7),'------I')==1
        Level4=erase(InterPro_v72_ParentChildTree_Dom(i,1),"-");
        InterPro_v72_Dom_groups(i,1)=Level1;
        InterPro_v72_Dom_groups(i,2)=Level2;
        InterPro_v72_Dom_groups(i,3)=Level3;
        InterPro_v72_Dom_groups(i,4)=Level4;
    elseif isequal(InterPro_v72_ParentChildTree_Dom{i,1}(1,1:9),'--------I')==1
        Level5=erase(InterPro_v72_ParentChildTree_Dom(i,1),"-");
        InterPro_v72_Dom_groups(i,1)=Level1;
        InterPro_v72_Dom_groups(i,2)=Level2;
        InterPro_v72_Dom_groups(i,3)=Level3;
        InterPro_v72_Dom_groups(i,4)=Level4;
        InterPro_v72_Dom_groups(i,5)=Level5;
    end
end
spc_ind=cellfun(@isempty,InterPro_v72_Dom_groups);
InterPro_v72_Dom_groups(spc_ind==1)=cellstr('X');
save DRUIDom_Files/InterPro_v72_Dom_groups.mat InterPro_v72_Dom_groups





% -------------------------------
% Domain-Compoung Mapping Process
% -------------------------------


% Associating compounds with single domains:

% Calculating the compound cluster (similarity based) single domain assocation scores for all compounds:

load DRUIDom_Files/Combined_Org_unique_act_inact_Prot_5_var.mat
load DRUIDom_Files/InterPro_v72_Dom_groups.mat
load DRUIDom_Files/InterPro_entriesmatched_onlytarget_Org_onlydom.mat
Domain_mapping_Prot_5.Compounds=Combined_Org_id_unique_Prot_5;
Domain_mapping_Prot_5.Domains=cell(length(Combined_Org_id_unique_Prot_5),1);
Domain_mapping_Prot_5.MapScr=cell(length(Combined_Org_id_unique_Prot_5),1);
for i=1:length(Combined_Org_id_unique_Prot_5)
    disp(['Compound number: ' num2str(i), ' / ' num2str(length(Combined_Org_id_unique_Prot_5))])
    Lia=ismember(InterPro_entriesmatched_onlytarget_Org_onlydom(:,1),Combined_Org_unique_act_Prot_5{i,1});
    InterPro_entriesmatched_temp_act=InterPro_entriesmatched_onlytarget_Org_onlydom(Lia==1,:);
    Dom_temp_act=unique(InterPro_entriesmatched_onlytarget_Org_onlydom(Lia==1,2));
    Lia2=ismember(InterPro_entriesmatched_onlytarget_Org_onlydom(:,1),Combined_Org_unique_inact_withact_Prot_5{i,1});
    InterPro_entriesmatched_temp_inact=InterPro_entriesmatched_onlytarget_Org_onlydom(Lia2==1,:);
    clear TP FN FP TN
    for j=1:length(Dom_temp_act)
        Lian=ismember(InterPro_v72_Dom_groups,Dom_temp_act{j,1});
        Dom_temp_act_group=unique(InterPro_v72_Dom_groups(sum(Lian,2)>0,:));
        Dom_temp_act_group(strncmp(Dom_temp_act_group(:,1),'X',1),:)=[];
        Dom_temp_act_group(:,strncmp(Dom_temp_act_group(1,:),'X',1))=[];
        Lia=ismember(InterPro_entriesmatched_temp_act(:,2),Dom_temp_act_group);
        TP(j,1)=length(unique(InterPro_entriesmatched_temp_act(Lia==1,1)));
        FN(j,1)=length(unique(Combined_Org_unique_act_Prot_5{i,1}))-TP(j,1);
        Lia2=ismember(InterPro_entriesmatched_temp_inact(:,2),Dom_temp_act_group);
        if sum(Lia2)~=0
            FP(j,1)=length(unique(InterPro_entriesmatched_temp_inact(Lia2==1,1)));
        else
            FP(j,1)=0;
        end
        TN(j,1)=length(unique(Combined_Org_unique_inact_withact_Prot_5{i,1}))-FP(j,1);
    end
    REC=TP./(TP+FN);
    PRE=TP./(TP+FP);
    ACCU=(TP+TN)./(TP+TN+FP+FN);
    F1=(2*TP)./(2*TP+FP+FN);
    MCC=((TP.*TN)-(FP.*FN))./(((TP+FP).*(TP+FN).*(TN+FP).*(TN+FN)).^(1/2));
    Domain_mapping_Prot_5.Domains{i,1}=Dom_temp_act;
    Domain_mapping_Prot_5.MapScr{i,1}=[TP FN FP TN REC PRE ACCU F1 MCC];
end
save DRUIDom_Files/Domain_mapping_Prot_5.mat Domain_mapping_Prot_5

% (organizing all mappings)
Domain_mapping_Prot_5_all_Org=cell(5000000,9);
to=1;
for i=1:length(Combined_Org_id_unique_Prot_5)
    disp(['Compound number: ' num2str(i), ' / ' num2str(length(Combined_Org_id_unique_Prot_5))])
    len_temp=size(Domain_mapping_Prot_5.MapScr{i,1},1);
    if len_temp~=0
        Domain_mapping_Prot_5_all_Org(to:(to+(len_temp)-1),1)=cellstr(repmat(Domain_mapping_Prot_5.Compounds{i,1},len_temp,1));
        Domain_mapping_Prot_5_all_Org(to:(to+(len_temp)-1),2)=Domain_mapping_Prot_5.Domains{i,1};
        Domain_mapping_Prot_5_all_Org(to:(to+(len_temp)-1),3:11)=num2cell(Domain_mapping_Prot_5.MapScr{i,1});
        to=to+(len_temp);
    end
end
del_ind=cellfun(@isempty,Domain_mapping_Prot_5_all_Org(:,1));
Domain_mapping_Prot_5_all_Org(del_ind,:)=[];
length(unique(Domain_mapping_Prot_5_all_Org(:,1)))
length(unique(Domain_mapping_Prot_5_all_Org(:,2)))
t=Domain_mapping_Prot_5_all_Org';
fid=fopen('DRUIDom_Files/Domain_mapping_Prot_5_all_Org.txt','w');
fprintf(fid,'%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n',t{:});
fclose(fid);
Domain_mapping_Prot_5_all_Org_part1=Domain_mapping_Prot_5_all_Org(1:1150000,:);
Domain_mapping_Prot_5_all_Org_part2=Domain_mapping_Prot_5_all_Org(1150001:2300000,:);
Domain_mapping_Prot_5_all_Org_part3=Domain_mapping_Prot_5_all_Org(2300001:end,:);
save DRUIDom_Files/Domain_mapping_Prot_5_all_Org_part1.mat Domain_mapping_Prot_5_all_Org_part1
save DRUIDom_Files/Domain_mapping_Prot_5_all_Org_part2.mat Domain_mapping_Prot_5_all_Org_part2
save DRUIDom_Files/Domain_mapping_Prot_5_all_Org_part3.mat Domain_mapping_Prot_5_all_Org_part3

% (filtering high score mappings)
Domain_mapping_Prot_5_all_MCC07=cell(0,9);
to=1;
for i=1:length(Combined_Org_id_unique_Prot_5)
    disp(['Compound number: ' num2str(i), ' / ' num2str(length(Combined_Org_id_unique_Prot_5))])
    ind_temp=find(Domain_mapping_Prot_5.MapScr{i,1}(:,9)>=0.7);
    if isempty(ind_temp)==0
        Domain_mapping_Prot_5_all_MCC07(to:(to+(length(ind_temp))-1),1)=cellstr(repmat(Domain_mapping_Prot_5.Compounds{i,1},length(ind_temp),1));
        Domain_mapping_Prot_5_all_MCC07(to:(to+(length(ind_temp))-1),2)=Domain_mapping_Prot_5.Domains{i,1}(ind_temp,1);
        Domain_mapping_Prot_5_all_MCC07(to:(to+(length(ind_temp))-1),3:11)=num2cell(Domain_mapping_Prot_5.MapScr{i,1}(ind_temp,1:9));
        to=to+(length(ind_temp));
    end
end
save DRUIDom_Files/Domain_mapping_Prot_5_all_MCC07.mat Domain_mapping_Prot_5_all_MCC07
length(unique(Domain_mapping_Prot_5_all_MCC07(:,1)))
length(unique(Domain_mapping_Prot_5_all_MCC07(:,2)))

Domain_mapping_Prot_5_all_Acc08=cell(0,9);
to=1;
for i=1:length(Combined_Org_id_unique_Prot_5)
    disp(['Compound number: ' num2str(i), ' / ' num2str(length(Combined_Org_id_unique_Prot_5))])
    ind_temp=find(Domain_mapping_Prot_5.MapScr{i,1}(:,7)>=0.8);
    if isempty(ind_temp)==0
        Domain_mapping_Prot_5_all_Acc08(to:(to+(length(ind_temp))-1),1)=cellstr(repmat(Domain_mapping_Prot_5.Compounds{i,1},length(ind_temp),1));
        Domain_mapping_Prot_5_all_Acc08(to:(to+(length(ind_temp))-1),2)=Domain_mapping_Prot_5.Domains{i,1}(ind_temp,1);
        Domain_mapping_Prot_5_all_Acc08(to:(to+(length(ind_temp))-1),3:11)=num2cell(Domain_mapping_Prot_5.MapScr{i,1}(ind_temp,1:9));
        to=to+(length(ind_temp));
    end
end
save DRUIDom_Files/Domain_mapping_Prot_5_all_Acc08.mat Domain_mapping_Prot_5_all_Acc08
length(unique(Domain_mapping_Prot_5_all_Acc08(:,1)))
length(unique(Domain_mapping_Prot_5_all_Acc08(:,2)))

Domain_mapping_Prot_5_all_F108=cell(0,9);
to=1;
for i=1:length(Combined_Org_id_unique_Prot_5)
    disp(['Compound number: ' num2str(i), ' / ' num2str(length(Combined_Org_id_unique_Prot_5))])
    ind_temp=find(Domain_mapping_Prot_5.MapScr{i,1}(:,8)>=0.8);
    if isempty(ind_temp)==0
        Domain_mapping_Prot_5_all_F108(to:(to+(length(ind_temp))-1),1)=cellstr(repmat(Domain_mapping_Prot_5.Compounds{i,1},length(ind_temp),1));
        Domain_mapping_Prot_5_all_F108(to:(to+(length(ind_temp))-1),2)=Domain_mapping_Prot_5.Domains{i,1}(ind_temp,1);
        Domain_mapping_Prot_5_all_F108(to:(to+(length(ind_temp))-1),3:11)=num2cell(Domain_mapping_Prot_5.MapScr{i,1}(ind_temp,1:9));
        to=to+(length(ind_temp));
    end
end
save DRUIDom_Files/Domain_mapping_Prot_5_all_F108.mat Domain_mapping_Prot_5_all_F108
length(unique(Domain_mapping_Prot_5_all_F108(:,1)))
length(unique(Domain_mapping_Prot_5_all_F108(:,2)))

% (according to the performance analysis at the end, multiple filters will be used: F1>=0.5, ACCU>=0.5, REC>=0.5, PRE>=0.5)
Domain_mapping_Prot_5_all_MulFil=cell(4000000,9);
to=1;
for i=1:length(Domain_mapping_Prot_5.Compounds)
    disp(['Compound number: ' num2str(i), ' / ' num2str(length(Domain_mapping_Prot_5.Compounds))])
    ind_temp=find(Domain_mapping_Prot_5.MapScr{i,1}(:,5)>=0.5 & Domain_mapping_Prot_5.MapScr{i,1}(:,6)>=0.5 & Domain_mapping_Prot_5.MapScr{i,1}(:,7)>=0.5 & Domain_mapping_Prot_5.MapScr{i,1}(:,8)>=0.5 & Domain_mapping_Prot_5.MapScr{i,1}(:,9)>=0);
    if isempty(ind_temp)==0
        Domain_mapping_Prot_5_all_MulFil(to:(to+(length(ind_temp))-1),1)=cellstr(repmat(Domain_mapping_Prot_5.Compounds{i,1},length(ind_temp),1));
        Domain_mapping_Prot_5_all_MulFil(to:(to+(length(ind_temp))-1),2)=Domain_mapping_Prot_5.Domains{i,1}(ind_temp,1);
        Domain_mapping_Prot_5_all_MulFil(to:(to+(length(ind_temp))-1),3:11)=num2cell(Domain_mapping_Prot_5.MapScr{i,1}(ind_temp,1:9));
        to=to+(length(ind_temp));
    end
end
del_ind=cellfun(@isempty,Domain_mapping_Prot_5_all_MulFil(:,1));
Domain_mapping_Prot_5_all_MulFil(del_ind,:)=[];
save DRUIDom_Files/Domain_mapping_Prot_5_all_MulFil.mat Domain_mapping_Prot_5_all_MulFil
length(unique(Domain_mapping_Prot_5_all_MulFil(:,1)))
length(unique(Domain_mapping_Prot_5_all_MulFil(:,2)))

% (manually analyzing the mappings generated by thresholding various performance metrics)
figure;hist(sort(cell2mat(Domain_mapping_Prot_5_all_F108(:,end)),'descend'),0:0.05:1)
figure;hist(sort(cell2mat(Domain_mapping_Prot_5_all_Acc08(:,end)),'descend'),0:0.05:1)
min(cell2mat(Domain_mapping_Prot_5_all_MCC07(:,7)))
min(cell2mat(Domain_mapping_Prot_5_all_MCC07(:,8)))
min(cell2mat(Domain_mapping_Prot_5_all_MCC07(:,9)))
min(cell2mat(Domain_mapping_Prot_5_all_MCC07(:,10)))
min(cell2mat(Domain_mapping_Prot_5_all_F108(:,9)))
Domain_mapping_Prot_5_all_Acc08(ismember(Domain_mapping_Prot_5_all_Acc08(:,2),'IPR000719')==1,:)
LIMK_kinasedom=Domain_mapping_Prot_5_all_MulFil(ismember(Domain_mapping_Prot_5_all_MulFil(:,2),'IPR000719')==1,:);
Domain_mapping_Prot_5_all_Acc08(ismember(Domain_mapping_Prot_5_all_Acc08(:,2),'IPR001478')==1,:)
Domain_mapping_Prot_5.MapScr{ismember(Domain_mapping_Prot_5.Compounds,'11223992')==1,1}(ismember(Domain_mapping_Prot_5.Domains{ismember(Domain_mapping_Prot_5.Compounds,'11223992')==1,1},'IPR039192')==1,:)



% Calculating the single compound (non-similarity based) single domain assocation scores for all compounds:

load DRUIDom_Files/Combined_Org_act_var.mat
load DRUIDom_Files/Combined_Org_inact_withact_var.mat
load DRUIDom_Files/Combined_Org_id_unique.mat
load DRUIDom_Files/Locations_act_inact.mat

% (generating the combined protein groups for each ligand)
parpool(4)
Combined_Org_unique_act_Prot_noSim=cell(length(Combined_Org_id_unique),1);
Combined_Org_unique_inact_withact_Prot_noSim=cell(length(Combined_Org_id_unique),1);
parfor i=1:length(Combined_Org_id_unique)
    disp(['Compound number: ' num2str(i), ' / ' num2str(length(Combined_Org_id_unique))])
    Combined_Org_thrSim_temp_ind=i;
    Combined_Org_unique_act_Prot_noSim{i,1}=unique(Combined_Org_act(ismember(Locb_act,Combined_Org_thrSim_temp_ind)==1,1));
    Combined_Org_unique_inact_withact_Prot_noSim{i,1}=unique(Combined_Org_inact_withact(ismember(Locb_inact,Combined_Org_thrSim_temp_ind)==1,1));
end
save DRUIDom_Files/Combined_Org_unique_act_Prot_noSim.mat Combined_Org_unique_act_Prot_noSim
save DRUIDom_Files/Combined_Org_unique_inact_withact_Prot_noSim.mat Combined_Org_unique_inact_withact_Prot_noSim

% (analyzing and thresholding the protein arrays)
Combined_Org_unique_Samplenum_act_inact_noSim_sup=zeros(length(Combined_Org_id_unique),3);
for i=1:length(Combined_Org_id_unique)
    disp(['Compound number: ' num2str(i), ' / ' num2str(length(Combined_Org_id_unique))])
    Combined_Org_unique_Samplenum_act_inact_noSim_sup(i,1)=length(Combined_Org_unique_act_Prot_noSim{i,1});
    Combined_Org_unique_Samplenum_act_inact_noSim_sup(i,2)=length(Combined_Org_unique_inact_withact_Prot_noSim{i,1});
    Combined_Org_unique_Samplenum_act_inact_noSim_sup(i,3)=(Combined_Org_unique_Samplenum_act_inact_noSim_sup(i,1)+Combined_Org_unique_Samplenum_act_inact_noSim_sup(i,2))/2;
end
save DRUIDom_Files/Combined_Org_unique_Samplenum_act_inact_noSim_sup.mat Combined_Org_unique_Samplenum_act_inact_noSim_sup

length(find(Combined_Org_unique_Samplenum_act_inact_noSim_sup(:,1)>=1))
length(find(Combined_Org_unique_Samplenum_act_inact_noSim_sup(:,2)>=1))
length(find(Combined_Org_unique_Samplenum_act_inact_noSim_sup(:,1)>=1 & Combined_Org_unique_Samplenum_act_inact_noSim_sup(:,2)>=1))
length(find(Combined_Org_unique_Samplenum_act_inact_noSim_sup(:,1)>=3))
length(find(Combined_Org_unique_Samplenum_act_inact_noSim_sup(:,2)>=3))
length(find(Combined_Org_unique_Samplenum_act_inact_noSim_sup(:,1)>=3 & Combined_Org_unique_Samplenum_act_inact_noSim_sup(:,2)>=3))
length(find(Combined_Org_unique_Samplenum_act_inact_noSim_sup(:,1)>=5))
length(find(Combined_Org_unique_Samplenum_act_inact_noSim_sup(:,2)>=5))
length(find(Combined_Org_unique_Samplenum_act_inact_noSim_sup(:,1)>=5 & Combined_Org_unique_Samplenum_act_inact_noSim_sup(:,2)>=5))
length(find(Combined_Org_unique_Samplenum_act_inact_noSim_sup(:,1)>=10))
length(find(Combined_Org_unique_Samplenum_act_inact_noSim_sup(:,2)>=10))
length(find(Combined_Org_unique_Samplenum_act_inact_noSim_sup(:,1)>=10 & Combined_Org_unique_Samplenum_act_inact_noSim_sup(:,2)>=10))

% (histograms of active and inactive distribution in the noSim dataset)
figure;
histogram(sort(2*Combined_Org_unique_Samplenum_act_inact_noSim_sup(:,3),'descend'),-0.5:1:(max(2*Combined_Org_unique_Samplenum_act_inact_noSim_sup(:,3))+0.5),'FaceAlpha',1)
set(gca,'YScale','log')
axis([0.5 250 0 1000000])
figure;
histogram(sort(Combined_Org_unique_Samplenum_act_inact_noSim_sup(:,1),'descend'),-0.5:1:(max(Combined_Org_unique_Samplenum_act_inact_noSim_sup(:,1))+0.5),'FaceColor',[0 0.7 0],'FaceAlpha',1)
set(gca,'YScale','log')
axis([0.5 250 0 1000000])
figure;
histogram(sort(Combined_Org_unique_Samplenum_act_inact_noSim_sup(:,2),'descend'),-0.5:1:(max(Combined_Org_unique_Samplenum_act_inact_noSim_sup(:,2))+0.5),'FaceColor',[0.8 0 0],'FaceAlpha',1)
set(gca,'YScale','log')
axis([0.5 250 0 1000000])


thr_act_samp=3;
thr_inact_samp=3;

Combined_Org_unique_Samplenum_act_inact_noSim_sup_3=Combined_Org_unique_Samplenum_act_inact_noSim_sup(Combined_Org_unique_Samplenum_act_inact_noSim_sup(:,1)>=thr_act_samp & Combined_Org_unique_Samplenum_act_inact_noSim_sup(:,2)>=thr_inact_samp,:);
Combined_Org_unique_act_Prot_noSim_3=Combined_Org_unique_act_Prot_noSim(Combined_Org_unique_Samplenum_act_inact_noSim_sup(:,1)>=thr_act_samp & Combined_Org_unique_Samplenum_act_inact_noSim_sup(:,2)>=thr_inact_samp,:);
Combined_Org_unique_inact_withact_Prot_noSim_3=Combined_Org_unique_inact_withact_Prot_noSim(Combined_Org_unique_Samplenum_act_inact_noSim_sup(:,1)>=thr_act_samp & Combined_Org_unique_Samplenum_act_inact_noSim_sup(:,2)>=thr_inact_samp,:);
Combined_Org_id_unique_Prot_noSim_3=Combined_Org_id_unique(Combined_Org_unique_Samplenum_act_inact_noSim_sup(:,1)>=thr_act_samp & Combined_Org_unique_Samplenum_act_inact_noSim_sup(:,2)>=thr_inact_samp,:);
save DRUIDom_Files/Combined_Org_unique_act_inact_Prot_noSim_3_var.mat Combined_Org_unique_Samplenum_act_inact_noSim_sup_3 Combined_Org_unique_act_Prot_noSim_3 Combined_Org_unique_inact_withact_Prot_noSim_3 Combined_Org_id_unique_Prot_noSim_3

% (calculating the single domain assocation scores for all compounds)
load DRUIDom_Files/Combined_Org_unique_act_inact_Prot_noSim_3_var.mat
load DRUIDom_Files/InterPro_v72_Dom_groups.mat
load DRUIDom_Files/InterPro_entriesmatched_onlytarget_Org_onlydom.mat
Domain_mapping_Prot_noSim_3.Compounds=Combined_Org_id_unique_Prot_noSim_3;
Domain_mapping_Prot_noSim_3.Domains=cell(length(Combined_Org_id_unique_Prot_noSim_3),1);
Domain_mapping_Prot_noSim_3.MapScr=cell(length(Combined_Org_id_unique_Prot_noSim_3),1);
for i=1:length(Combined_Org_id_unique_Prot_noSim_3)
    disp(['Compound number: ' num2str(i), ' / ' num2str(length(Combined_Org_id_unique_Prot_noSim_3))])
    Lia=ismember(InterPro_entriesmatched_onlytarget_Org_onlydom(:,1),Combined_Org_unique_act_Prot_noSim_3{i,1});
    InterPro_entriesmatched_temp_act=InterPro_entriesmatched_onlytarget_Org_onlydom(Lia==1,:);
    Dom_temp_act=unique(InterPro_entriesmatched_onlytarget_Org_onlydom(Lia==1,2));
    Lia2=ismember(InterPro_entriesmatched_onlytarget_Org_onlydom(:,1),Combined_Org_unique_inact_withact_Prot_noSim_3{i,1});
    InterPro_entriesmatched_temp_inact=InterPro_entriesmatched_onlytarget_Org_onlydom(Lia2==1,:);
    clear TP FN FP TN
    if isempty(Dom_temp_act)==0
        for j=1:length(Dom_temp_act)
            Lian=ismember(InterPro_v72_Dom_groups,Dom_temp_act{j,1});
            Dom_temp_act_group=unique(InterPro_v72_Dom_groups(sum(Lian,2)>0,:));
            Dom_temp_act_group(strncmp(Dom_temp_act_group(:,1),'X',1),:)=[];
            Dom_temp_act_group(:,strncmp(Dom_temp_act_group(1,:),'X',1))=[];
            Lia=ismember(InterPro_entriesmatched_temp_act(:,2),Dom_temp_act_group);
            TP(j,1)=length(unique(InterPro_entriesmatched_temp_act(Lia==1,1)));
            FN(j,1)=length(unique(Combined_Org_unique_act_Prot_noSim_3{i,1}))-TP(j,1);
            Lia2=ismember(InterPro_entriesmatched_temp_inact(:,2),Dom_temp_act_group);
            if sum(Lia2)~=0
                FP(j,1)=length(unique(InterPro_entriesmatched_temp_inact(Lia2==1,1)));
            else
                FP(j,1)=0;
            end
            TN(j,1)=length(unique(Combined_Org_unique_inact_withact_Prot_noSim_3{i,1}))-FP(j,1);
        end
        REC=TP./(TP+FN);
        PRE=TP./(TP+FP);
        ACCU=(TP+TN)./(TP+TN+FP+FN);
        F1=(2*TP)./(2*TP+FP+FN);
        MCC=((TP.*TN)-(FP.*FN))./(((TP+FP).*(TP+FN).*(TN+FP).*(TN+FN)).^(1/2));
        Domain_mapping_Prot_noSim_3.Domains{i,1}=Dom_temp_act;
        Domain_mapping_Prot_noSim_3.MapScr{i,1}=[TP FN FP TN REC PRE ACCU F1 MCC];
    end
end
save DRUIDom_Files/Domain_mapping_Prot_noSim_3.mat Domain_mapping_Prot_noSim_3

% (organizing all mappings)
Domain_mapping_Prot_noSim_3_all_Org=cell(5000000,9);
to=1;
for i=1:length(Combined_Org_id_unique_Prot_noSim_3)
    disp(['Compound number: ' num2str(i), ' / ' num2str(length(Combined_Org_id_unique_Prot_noSim_3))])
    len_temp=size(Domain_mapping_Prot_noSim_3.MapScr{i,1},1);
    if len_temp~=0
        Domain_mapping_Prot_noSim_3_all_Org(to:(to+(len_temp)-1),1)=cellstr(repmat(Domain_mapping_Prot_noSim_3.Compounds{i,1},len_temp,1));
        Domain_mapping_Prot_noSim_3_all_Org(to:(to+(len_temp)-1),2)=Domain_mapping_Prot_noSim_3.Domains{i,1};
        Domain_mapping_Prot_noSim_3_all_Org(to:(to+(len_temp)-1),3:11)=num2cell(Domain_mapping_Prot_noSim_3.MapScr{i,1});
        to=to+(len_temp);
    end
end
del_ind=cellfun(@isempty,Domain_mapping_Prot_noSim_3_all_Org(:,1));
Domain_mapping_Prot_noSim_3_all_Org(del_ind,:)=[];
length(unique(Domain_mapping_Prot_noSim_3_all_Org(:,1)))
length(unique(Domain_mapping_Prot_noSim_3_all_Org(:,2)))
t=Domain_mapping_Prot_noSim_3_all_Org';
fid=fopen('DRUIDom_Files/Domain_mapping_Prot_noSim_3_all_Org.txt','w');
fprintf(fid,'%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n',t{:});
fclose(fid);
save DRUIDom_Files/Domain_mapping_Prot_noSim_3_all_Org.mat Domain_mapping_Prot_noSim_3_all_Org

% (filtering high score mappings, using the same filers are previously: MCC>=0, F1>=0.5, ACCU>=0.5, REC>=0.5, PRE>=0.5)
Domain_mapping_Prot_noSim_3_all_MulFil=cell(4000000,9);
to=1;
for i=1:length(Domain_mapping_Prot_noSim_3.Compounds)
    disp(['Compound number: ' num2str(i), ' / ' num2str(length(Domain_mapping_Prot_noSim_3.Compounds))])
    if isempty(Domain_mapping_Prot_noSim_3.MapScr{i,1})==0
        ind_temp=find(Domain_mapping_Prot_noSim_3.MapScr{i,1}(:,5)>=0.5 & Domain_mapping_Prot_noSim_3.MapScr{i,1}(:,6)>=0.5 & Domain_mapping_Prot_noSim_3.MapScr{i,1}(:,7)>=0.5 & Domain_mapping_Prot_noSim_3.MapScr{i,1}(:,8)>=0.5 & Domain_mapping_Prot_noSim_3.MapScr{i,1}(:,9)>=0);
    end
    if isempty(ind_temp)==0
        Domain_mapping_Prot_noSim_3_all_MulFil(to:(to+(length(ind_temp))-1),1)=cellstr(repmat(Domain_mapping_Prot_noSim_3.Compounds{i,1},length(ind_temp),1));
        Domain_mapping_Prot_noSim_3_all_MulFil(to:(to+(length(ind_temp))-1),2)=Domain_mapping_Prot_noSim_3.Domains{i,1}(ind_temp,1);
        Domain_mapping_Prot_noSim_3_all_MulFil(to:(to+(length(ind_temp))-1),3:11)=num2cell(Domain_mapping_Prot_noSim_3.MapScr{i,1}(ind_temp,1:9));
        to=to+(length(ind_temp));
    end
    ind_temp=[];
end
del_ind=cellfun(@isempty,Domain_mapping_Prot_noSim_3_all_MulFil(:,1));
Domain_mapping_Prot_noSim_3_all_MulFil(del_ind,:)=[];
save DRUIDom_Files/Domain_mapping_Prot_noSim_3_all_MulFil.mat Domain_mapping_Prot_noSim_3_all_MulFil
length(Domain_mapping_Prot_noSim_3_all_MulFil)
length(unique(Domain_mapping_Prot_noSim_3_all_MulFil(:,1)))
length(unique(Domain_mapping_Prot_noSim_3_all_MulFil(:,2)))



% Merging compound similarity based domain associations with noSim associations:

load DRUIDom_Files/Domain_mapping_Prot_5_all_MulFil.mat
load DRUIDom_Files/Domain_mapping_Prot_noSim_3_all_MulFil.mat
% (extracting mutual associations and selecting maximum performance according to MCC for merging)
Domain_mapping_Prot_5_all_MulFil_23=strcat(Domain_mapping_Prot_5_all_MulFil(:,1), {' '}, Domain_mapping_Prot_5_all_MulFil(:,2));
Domain_mapping_Prot_noSim_3_all_MulFil_23=strcat(Domain_mapping_Prot_noSim_3_all_MulFil(:,1), {' '}, Domain_mapping_Prot_noSim_3_all_MulFil(:,2));
[~,I1,I2]=intersect(Domain_mapping_Prot_5_all_MulFil_23,Domain_mapping_Prot_noSim_3_all_MulFil_23);

C1=(cell2mat(Domain_mapping_Prot_5_all_MulFil(I1,end))>cell2mat(Domain_mapping_Prot_noSim_3_all_MulFil(I2,end)));
mer1=Domain_mapping_Prot_5_all_MulFil(I1(C1==1),:);
mer2=Domain_mapping_Prot_noSim_3_all_MulFil(I2(C1==0),:);
I1inv=setdiff((1:length(Domain_mapping_Prot_5_all_MulFil)),I1)';
I2inv=setdiff((1:length(Domain_mapping_Prot_noSim_3_all_MulFil)),I2)';

% (saving the merged and filtered domain accosiation file)
Domain_mapping_Prot_Sim5_noSim3_MulFil_merged=[Domain_mapping_Prot_5_all_MulFil(I1inv,:);Domain_mapping_Prot_noSim_3_all_MulFil(I2inv,:);mer1;mer2];
save DRUIDom_Files/Domain_mapping_Prot_Sim5_noSim3_MulFil_merged.mat Domain_mapping_Prot_Sim5_noSim3_MulFil_merged
t=Domain_mapping_Prot_Sim5_noSim3_MulFil_merged';
fid=fopen('DRUIDom_Files/Domain_mapping_Prot_Sim5_noSim3_MulFil_merged_finalized.txt','w');
fprintf(fid,'%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n',t{:});
fclose(fid);
length(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged)
length(unique(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged(:,1)))
length(unique(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged(:,2)))
% (final domain association stats: 27,032 mappings between 8,165 compounds and 250 InterPro domains)
fid=fopen('DRUIDom_finalized_domain_compound_mappings.txt','w');
fprintf(fid,'%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n',t{:});
fclose(fid);



% Associating compounds with domain pairs (where it performs better
% compared to mappings to single domains in the corresponding pairs):

% (generating an array for InterPro entry pairs in the same hierarchy, so that these pairs can be subtracted from all domain pairs since these pairs are not 2 domains hits but 2 related entries hit to the same position)
load DRUIDom_Files/InterPro_v72_Dom_groups.mat
InterPro_domain_relation_pairs=cell(0,1);
to=1;
for i=1:length(InterPro_v72_Dom_groups)
    if isequal(InterPro_v72_Dom_groups{i,2},'X')~=1
        temp=InterPro_v72_Dom_groups(i,ismember(InterPro_v72_Dom_groups(i,:),({'X'}))==0);
        n=length(temp);
        a=fliplr(fullfact([n n]));
        a(~diff(a')',:)=[];
        temp_pair=temp(a);
        x=strcat(temp_pair(:,1), {'-'}, temp_pair(:,2));
        InterPro_domain_relation_pairs(to:(to+(size(x,1))-1),1)=x;
        to=to+(size(x,1));
    end
end
InterPro_domain_relation_pairs=unique(InterPro_domain_relation_pairs);
save DRUIDom_Files/InterPro_domain_relation_pairs.mat InterPro_domain_relation_pairs

% (generating the 'bag of domain' domain pairs for our targets)
InterPro_entriesmatched_onlytarget_Org_bagofdom_onlydompairs=cell(0,2);
to=1;
Prot_unique=unique(InterPro_entriesmatched_onlytarget_Org_onlydom(:,1));
for i=1:length(Prot_unique)
    disp(['Line number: ' num2str(i), ' / ' num2str(length(Prot_unique))])
    dom_list_temp=InterPro_entriesmatched_onlytarget_Org_onlydom(ismember(InterPro_entriesmatched_onlytarget_Org_onlydom(:,1),Prot_unique(i,1))==1,2);
    n=length(dom_list_temp);
    if n>1
        a=fliplr(fullfact([n n]));
        a(~diff(a')',:)=[];
        dompair_list_temp=dom_list_temp(a);
        x=strcat(dompair_list_temp(:,1), {'-'}, dompair_list_temp(:,2));
        del_ind=ismember(x,InterPro_domain_relation_pairs);
        x(del_ind==1,:)=[];
        if isempty(x)==0
            InterPro_entriesmatched_onlytarget_Org_bagofdom_onlydompairs(to:(to+(size(x,1))-1),1)=repmat(Prot_unique(i,1),size(x,1),1);
            InterPro_entriesmatched_onlytarget_Org_bagofdom_onlydompairs(to:(to+(size(x,1))-1),2)=x;
            to=to+(size(x,1));
        end
    end
end
save DRUIDom_Files/InterPro_entriesmatched_onlytarget_Org_bagofdom_onlydompairs.mat InterPro_entriesmatched_onlytarget_Org_bagofdom_onlydompairs

% (generating the 'bag of domain' domain pairs for all human proteins)
InterPro_entriesmatched_allhuman_Org_bagofdom_onlydompairs=cell(0,2);
to=1;
Prot_unique=unique(InterPro_entriesmatched_allhuman_Org_onlydom(:,1));
for i=1:length(Prot_unique)
    disp(['Line number: ' num2str(i), ' / ' num2str(length(Prot_unique))])
    dom_list_temp=InterPro_entriesmatched_allhuman_Org_onlydom(ismember(InterPro_entriesmatched_allhuman_Org_onlydom(:,1),Prot_unique(i,1))==1,2);
    n=length(dom_list_temp);
    if n>1
        a=fliplr(fullfact([n n]));
        a(~diff(a')',:)=[];
        dompair_list_temp=dom_list_temp(a);
        x=strcat(dompair_list_temp(:,1), {'-'}, dompair_list_temp(:,2));
        del_ind=ismember(x,InterPro_domain_relation_pairs);
        x(del_ind==1,:)=[];
        if isempty(x)==0
            InterPro_entriesmatched_allhuman_Org_bagofdom_onlydompairs(to:(to+(size(x,1))-1),1)=repmat(Prot_unique(i,1),size(x,1),1);
            InterPro_entriesmatched_allhuman_Org_bagofdom_onlydompairs(to:(to+(size(x,1))-1),2)=x;
            to=to+(size(x,1));
        end
    end
end
save DRUIDom_Files/InterPro_entriesmatched_allhuman_Org_bagofdom_onlydompairs.mat InterPro_entriesmatched_allhuman_Org_bagofdom_onlydompairs



% Calculating the compound cluster (similarity based) domain pair association scores for all compounds:

load DRUIDom_Files/InterPro_entriesmatched_onlytarget_Org_bagofdom_onlydompairs.mat
load DRUIDom_Files/Combined_Org_unique_act_inact_Prot_5_var.mat
DomainPair_mapping_Prot_5.Compounds=Combined_Org_id_unique_Prot_5;
DomainPair_mapping_Prot_5.Domains=cell(length(Combined_Org_id_unique_Prot_5),1);
DomainPair_mapping_Prot_5.MapScr=cell(length(Combined_Org_id_unique_Prot_5),1);
for i=1:length(Combined_Org_id_unique_Prot_5)
    disp(['Compound number: ' num2str(i), ' / ' num2str(length(Combined_Org_id_unique_Prot_5))])
    Lia=ismember(InterPro_entriesmatched_onlytarget_Org_bagofdom_onlydompairs(:,1),Combined_Org_unique_act_Prot_5{i,1});
    InterPro_entriesmatched_temp_act=InterPro_entriesmatched_onlytarget_Org_bagofdom_onlydompairs(Lia==1,:);
    Dom_temp_act=unique(InterPro_entriesmatched_onlytarget_Org_bagofdom_onlydompairs(Lia==1,2));
    Lia2=ismember(InterPro_entriesmatched_onlytarget_Org_bagofdom_onlydompairs(:,1),Combined_Org_unique_inact_withact_Prot_5{i,1});
    InterPro_entriesmatched_temp_inact=InterPro_entriesmatched_onlytarget_Org_bagofdom_onlydompairs(Lia2==1,:);
    clear TP FN FP TN
    if isempty(Dom_temp_act)==0
        for j=1:length(Dom_temp_act)
            Lia=ismember(InterPro_entriesmatched_temp_act(:,2),Dom_temp_act(j,1));
            TP(j,1)=length(unique(InterPro_entriesmatched_temp_act(Lia==1,1)));
            FN(j,1)=length(unique(Combined_Org_unique_act_Prot_5{i,1}))-TP(j,1);
            Lia2=ismember(InterPro_entriesmatched_temp_inact(:,2),Dom_temp_act(j,1));
            if sum(Lia2)~=0
                FP(j,1)=length(unique(InterPro_entriesmatched_temp_inact(Lia2==1,1)));
            else
                FP(j,1)=0;
            end
            TN(j,1)=length(unique(Combined_Org_unique_inact_withact_Prot_5{i,1}))-FP(j,1);
        end
        REC=TP./(TP+FN);
        PRE=TP./(TP+FP);
        ACCU=(TP+TN)./(TP+TN+FP+FN);
        F1=(2*TP)./(2*TP+FP+FN);
        MCC=((TP.*TN)-(FP.*FN))./(((TP+FP).*(TP+FN).*(TN+FP).*(TN+FN)).^(1/2));
        DomainPair_mapping_Prot_5.Domains{i,1}=Dom_temp_act;
        DomainPair_mapping_Prot_5.MapScr{i,1}=[TP FN FP TN REC PRE ACCU F1 MCC];
    end
end
save DRUIDom_Files/DomainPair_mapping_Prot_5.mat DomainPair_mapping_Prot_5

% (organizing all compound-domain pair mappings)
DomainPair_mapping_Prot_5_all_Org=cell(25000000,9);
to=1;
for i=1:length(Combined_Org_id_unique_Prot_5)
    disp(['Compound number: ' num2str(i), ' / ' num2str(length(Combined_Org_id_unique_Prot_5))])
    len_temp=size(DomainPair_mapping_Prot_5.MapScr{i,1},1);
    if len_temp~=0
        DomainPair_mapping_Prot_5_all_Org(to:(to+(len_temp)-1),1)=cellstr(repmat(DomainPair_mapping_Prot_5.Compounds{i,1},len_temp,1));
        DomainPair_mapping_Prot_5_all_Org(to:(to+(len_temp)-1),2)=DomainPair_mapping_Prot_5.Domains{i,1};
        DomainPair_mapping_Prot_5_all_Org(to:(to+(len_temp)-1),3:11)=num2cell(DomainPair_mapping_Prot_5.MapScr{i,1});
        to=to+(len_temp);
    end
end
del_ind=cellfun(@isempty,DomainPair_mapping_Prot_5_all_Org(:,1));
DomainPair_mapping_Prot_5_all_Org(del_ind,:)=[];
length(unique(DomainPair_mapping_Prot_5_all_Org(:,1)))
length(unique(DomainPair_mapping_Prot_5_all_Org(:,2)))
t=DomainPair_mapping_Prot_5_all_Org';
fid=fopen('DRUIDom_Files/DomainPair_mapping_Prot_5_all_Org.txt','w');
fprintf(fid,'%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n',t{:});
fclose(fid);
DomainPair_mapping_Prot_5_all_Org_part1=DomainPair_mapping_Prot_5_all_Org(1:1400000,:);
DomainPair_mapping_Prot_5_all_Org_part2=DomainPair_mapping_Prot_5_all_Org(1400001:2800000,:);
DomainPair_mapping_Prot_5_all_Org_part3=DomainPair_mapping_Prot_5_all_Org(2800001:4200000,:);
DomainPair_mapping_Prot_5_all_Org_part4=DomainPair_mapping_Prot_5_all_Org(4200001:5400000,:);
DomainPair_mapping_Prot_5_all_Org_part5=DomainPair_mapping_Prot_5_all_Org(5400001:6600000,:);
DomainPair_mapping_Prot_5_all_Org_part6=DomainPair_mapping_Prot_5_all_Org(6600001:7800000,:);
DomainPair_mapping_Prot_5_all_Org_part7=DomainPair_mapping_Prot_5_all_Org(7800001:9000000,:);
DomainPair_mapping_Prot_5_all_Org_part8=DomainPair_mapping_Prot_5_all_Org(9000001:end,:);
save DRUIDom_Files/DomainPair_mapping_Prot_5_all_Org_part1.mat DomainPair_mapping_Prot_5_all_Org_part1
save DRUIDom_Files/DomainPair_mapping_Prot_5_all_Org_part2.mat DomainPair_mapping_Prot_5_all_Org_part2
save DRUIDom_Files/DomainPair_mapping_Prot_5_all_Org_part3.mat DomainPair_mapping_Prot_5_all_Org_part3
save DRUIDom_Files/DomainPair_mapping_Prot_5_all_Org_part4.mat DomainPair_mapping_Prot_5_all_Org_part4
save DRUIDom_Files/DomainPair_mapping_Prot_5_all_Org_part5.mat DomainPair_mapping_Prot_5_all_Org_part5
save DRUIDom_Files/DomainPair_mapping_Prot_5_all_Org_part6.mat DomainPair_mapping_Prot_5_all_Org_part6
save DRUIDom_Files/DomainPair_mapping_Prot_5_all_Org_part7.mat DomainPair_mapping_Prot_5_all_Org_part7
save DRUIDom_Files/DomainPair_mapping_Prot_5_all_Org_part8.mat DomainPair_mapping_Prot_5_all_Org_part8

% (filtering high score mappings using the same filters as previously: MCC>=0, F1>=0.5, ACCU>=0.5, REC>=0.5, PRE>=0.5)
DomainPair_mapping_Prot_5_all_MulFil=cell(4000000,9);
to=1;
for i=1:length(DomainPair_mapping_Prot_5.Compounds)
    disp(['Compound number: ' num2str(i), ' / ' num2str(length(DomainPair_mapping_Prot_5.Compounds))])
    if isempty(DomainPair_mapping_Prot_5.MapScr{i,1})==0
        ind_temp=find(DomainPair_mapping_Prot_5.MapScr{i,1}(:,5)>=0.5 & DomainPair_mapping_Prot_5.MapScr{i,1}(:,6)>=0.5 & DomainPair_mapping_Prot_5.MapScr{i,1}(:,7)>=0.5 & DomainPair_mapping_Prot_5.MapScr{i,1}(:,8)>=0.5 & DomainPair_mapping_Prot_5.MapScr{i,1}(:,9)>=0);
        if isempty(ind_temp)==0
            DomainPair_mapping_Prot_5_all_MulFil(to:(to+(length(ind_temp))-1),1)=cellstr(repmat(DomainPair_mapping_Prot_5.Compounds{i,1},length(ind_temp),1));
            DomainPair_mapping_Prot_5_all_MulFil(to:(to+(length(ind_temp))-1),2)=DomainPair_mapping_Prot_5.Domains{i,1}(ind_temp,1);
            DomainPair_mapping_Prot_5_all_MulFil(to:(to+(length(ind_temp))-1),3:11)=num2cell(DomainPair_mapping_Prot_5.MapScr{i,1}(ind_temp,1:9));
            to=to+(length(ind_temp));
        end
    end
end
del_ind=cellfun(@isempty,DomainPair_mapping_Prot_5_all_MulFil(:,1));
DomainPair_mapping_Prot_5_all_MulFil(del_ind,:)=[];
save DRUIDom_Files/DomainPair_mapping_Prot_5_all_MulFil.mat DomainPair_mapping_Prot_5_all_MulFil
length(unique(DomainPair_mapping_Prot_5_all_MulFil(:,1)))
length(unique(DomainPair_mapping_Prot_5_all_MulFil(:,2)))



% Calculating the single compound (non-similarity based) domain pair assocation scores for all compounds:

load DRUIDom_Files/Combined_Org_unique_act_inact_Prot_noSim_3_var.mat
load DRUIDom_Files/InterPro_entriesmatched_onlytarget_Org_bagofdom_onlydompairs.mat
DomainPair_mapping_Prot_noSim_3.Compounds=Combined_Org_id_unique_Prot_noSim_3;
DomainPair_mapping_Prot_noSim_3.Domains=cell(length(Combined_Org_id_unique_Prot_noSim_3),1);
DomainPair_mapping_Prot_noSim_3.MapScr=cell(length(Combined_Org_id_unique_Prot_noSim_3),1);
for i=1:length(Combined_Org_id_unique_Prot_noSim_3)
    disp(['Compound number: ' num2str(i), ' / ' num2str(length(Combined_Org_id_unique_Prot_noSim_3))])
    Lia=ismember(InterPro_entriesmatched_onlytarget_Org_bagofdom_onlydompairs(:,1),Combined_Org_unique_act_Prot_noSim_3{i,1});
    InterPro_entriesmatched_temp_act=InterPro_entriesmatched_onlytarget_Org_bagofdom_onlydompairs(Lia==1,:);
    Dom_temp_act=unique(InterPro_entriesmatched_onlytarget_Org_bagofdom_onlydompairs(Lia==1,2));
    Lia2=ismember(InterPro_entriesmatched_onlytarget_Org_bagofdom_onlydompairs(:,1),Combined_Org_unique_inact_withact_Prot_noSim_3{i,1});
    InterPro_entriesmatched_temp_inact=InterPro_entriesmatched_onlytarget_Org_bagofdom_onlydompairs(Lia2==1,:);
    clear TP FN FP TN
    if isempty(Dom_temp_act)==0
        for j=1:length(Dom_temp_act)
            Lia=ismember(InterPro_entriesmatched_temp_act(:,2),Dom_temp_act(j,1));
            TP(j,1)=length(unique(InterPro_entriesmatched_temp_act(Lia==1,1)));
            FN(j,1)=length(unique(Combined_Org_unique_act_Prot_noSim_3{i,1}))-TP(j,1);
            Lia2=ismember(InterPro_entriesmatched_temp_inact(:,2),Dom_temp_act(j,1));
            if sum(Lia2)~=0
                FP(j,1)=length(unique(InterPro_entriesmatched_temp_inact(Lia2==1,1)));
            else
                FP(j,1)=0;
            end
            TN(j,1)=length(unique(Combined_Org_unique_inact_withact_Prot_noSim_3{i,1}))-FP(j,1);
        end
        REC=TP./(TP+FN);
        PRE=TP./(TP+FP);
        ACCU=(TP+TN)./(TP+TN+FP+FN);
        F1=(2*TP)./(2*TP+FP+FN);
        MCC=((TP.*TN)-(FP.*FN))./(((TP+FP).*(TP+FN).*(TN+FP).*(TN+FN)).^(1/2));
        DomainPair_mapping_Prot_noSim_3.Domains{i,1}=Dom_temp_act;
        DomainPair_mapping_Prot_noSim_3.MapScr{i,1}=[TP FN FP TN REC PRE ACCU F1 MCC];
    end
end
save DRUIDom_Files/DomainPair_mapping_Prot_noSim_3.mat DomainPair_mapping_Prot_noSim_3

% (organizing all mappings)
DomainPair_mapping_Prot_noSim_3_all_Org=cell(5000000,9);
to=1;
for i=1:length(Combined_Org_id_unique_Prot_noSim_3)
    disp(['Compound number: ' num2str(i), ' / ' num2str(length(Combined_Org_id_unique_Prot_noSim_3))])
    len_temp=size(DomainPair_mapping_Prot_noSim_3.MapScr{i,1},1);
    if len_temp~=0
        DomainPair_mapping_Prot_noSim_3_all_Org(to:(to+(len_temp)-1),1)=cellstr(repmat(DomainPair_mapping_Prot_noSim_3.Compounds{i,1},len_temp,1));
        DomainPair_mapping_Prot_noSim_3_all_Org(to:(to+(len_temp)-1),2)=DomainPair_mapping_Prot_noSim_3.Domains{i,1};
        DomainPair_mapping_Prot_noSim_3_all_Org(to:(to+(len_temp)-1),3:11)=num2cell(DomainPair_mapping_Prot_noSim_3.MapScr{i,1});
        to=to+(len_temp);
    end
end
del_ind=cellfun(@isempty,DomainPair_mapping_Prot_noSim_3_all_Org(:,1));
DomainPair_mapping_Prot_noSim_3_all_Org(del_ind,:)=[];
length(unique(DomainPair_mapping_Prot_noSim_3_all_Org(:,1)))
length(unique(DomainPair_mapping_Prot_noSim_3_all_Org(:,2)))
t=DomainPair_mapping_Prot_noSim_3_all_Org';
fid=fopen('DRUIDom_Files/DomainPair_mapping_Prot_noSim_3_all_Org.txt','w');
fprintf(fid,'%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n',t{:});
fclose(fid);
save DRUIDom_Files/DomainPair_mapping_Prot_noSim_3_all_Org.mat DomainPair_mapping_Prot_noSim_3_all_Org

% (filtering high score mappings using the same filters as previously: MCC>=0, F1>=0.5, ACCU>=0.5, REC>=0.5, PRE>=0.5)
DomainPair_mapping_Prot_noSim_3_all_MulFil=cell(4000000,9);
to=1;
for i=1:length(DomainPair_mapping_Prot_noSim_3.Compounds)
    disp(['Compound number: ' num2str(i), ' / ' num2str(length(DomainPair_mapping_Prot_noSim_3.Compounds))])
    if isempty(DomainPair_mapping_Prot_noSim_3.MapScr{i,1})==0
        ind_temp=find(DomainPair_mapping_Prot_noSim_3.MapScr{i,1}(:,5)>=0.5 & DomainPair_mapping_Prot_noSim_3.MapScr{i,1}(:,6)>=0.5 & DomainPair_mapping_Prot_noSim_3.MapScr{i,1}(:,7)>=0.5 & DomainPair_mapping_Prot_noSim_3.MapScr{i,1}(:,8)>=0.5 & DomainPair_mapping_Prot_noSim_3.MapScr{i,1}(:,9)>=0);
    end
    if isempty(ind_temp)==0
        DomainPair_mapping_Prot_noSim_3_all_MulFil(to:(to+(length(ind_temp))-1),1)=cellstr(repmat(DomainPair_mapping_Prot_noSim_3.Compounds{i,1},length(ind_temp),1));
        DomainPair_mapping_Prot_noSim_3_all_MulFil(to:(to+(length(ind_temp))-1),2)=DomainPair_mapping_Prot_noSim_3.Domains{i,1}(ind_temp,1);
        DomainPair_mapping_Prot_noSim_3_all_MulFil(to:(to+(length(ind_temp))-1),3:11)=num2cell(DomainPair_mapping_Prot_noSim_3.MapScr{i,1}(ind_temp,1:9));
        to=to+(length(ind_temp));
    end
    ind_temp=[];
end
del_ind=cellfun(@isempty,DomainPair_mapping_Prot_noSim_3_all_MulFil(:,1));
DomainPair_mapping_Prot_noSim_3_all_MulFil(del_ind,:)=[];
save DRUIDom_Files/DomainPair_mapping_Prot_noSim_3_all_MulFil.mat DomainPair_mapping_Prot_noSim_3_all_MulFil
length(unique(DomainPair_mapping_Prot_noSim_3_all_MulFil(:,1)))
length(unique(DomainPair_mapping_Prot_noSim_3_all_MulFil(:,2)))



% Merging compound similarity based domain pair associations with noSim associations:

load DRUIDom_Files/DomainPair_mapping_Prot_5_all_MulFil.mat
load DRUIDom_Files/DomainPair_mapping_Prot_noSim_3_all_MulFil.mat
% (extracting mutual associations and selecting maximum performance according to MCC for merging)
DomainPair_mapping_Prot_5_all_MulFil_23=strcat(DomainPair_mapping_Prot_5_all_MulFil(:,1), {' '}, DomainPair_mapping_Prot_5_all_MulFil(:,2));
DomainPair_mapping_Prot_noSim_3_all_MulFil_23=strcat(DomainPair_mapping_Prot_noSim_3_all_MulFil(:,1), {' '}, DomainPair_mapping_Prot_noSim_3_all_MulFil(:,2));
[~,I1,I2]=intersect(DomainPair_mapping_Prot_5_all_MulFil_23,DomainPair_mapping_Prot_noSim_3_all_MulFil_23);

C1=(cell2mat(DomainPair_mapping_Prot_5_all_MulFil(I1,end))>cell2mat(DomainPair_mapping_Prot_noSim_3_all_MulFil(I2,end)));
mer1=DomainPair_mapping_Prot_5_all_MulFil(I1(C1==1),:);
mer2=DomainPair_mapping_Prot_noSim_3_all_MulFil(I2(C1==0),:);
I1inv=setdiff((1:length(DomainPair_mapping_Prot_5_all_MulFil)),I1)';
I2inv=setdiff((1:length(DomainPair_mapping_Prot_noSim_3_all_MulFil)),I2)';

% (merging the filtered domain pair association file and removing duplicate rows that correspond to the same domain pair)
DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged=[DomainPair_mapping_Prot_5_all_MulFil(I1inv,:);DomainPair_mapping_Prot_noSim_3_all_MulFil(I2inv,:);mer1;mer2];
DomainPair_merged_temp=cell(0,11);
comp_uniq=unique(DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged(:,1));
for i=1:length(comp_uniq)
   ind=find(ismember(DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged(:,1),comp_uniq(i,1))==1);
   all_temp=DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged(ind,:);
   dom_temp=DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged(ind,2);
   clear dom_temp_on
   for k=1:size(dom_temp,1)
       dom_temp_on(k,1:2)=strsplit(dom_temp{k,1},'-')';
   end
   dom_temp_rev=strcat(dom_temp_on(:,2), {'-'}, dom_temp_on(:,1));
   del_temp=[];
   [~,Locb]=ismember(dom_temp_rev,dom_temp);
   for l=1:length(Locb)
       if Locb(Locb(l))==l
           if Locb(l)~=l
               Locb(l)=0;
           end
       end
   end
   Locb=Locb(Locb>0);
   DomainPair_merged_temp=[DomainPair_merged_temp;all_temp(Locb,:)];
end
DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged=DomainPair_merged_temp;
save DRUIDom_Files/DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged.mat DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged
t=DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged';
fid=fopen('DRUIDom_Files/DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged_finalized.txt','w');
fprintf(fid,'%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n',t{:});
fclose(fid);
length(DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged)
length(unique(DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged(:,1)))
length(unique(DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged(:,2)))
% (final domain pair association stats: 3,721 mappings between 1,456 compounds and 270 InterPro domain pairs)



% Filtering single domain and domain pair mappings by comparing and removing the corresponding single or pair mapping to the same compound when the performance is inferior:

load DRUIDom_Files/Domain_mapping_Prot_Sim5_noSim3_MulFil_merged.mat
load DRUIDom_Files/DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged.mat
del_realind1=[];del_realind2=[];Examples_DomPair_better=cell(0,11);
[IntSec,I1,I2]=intersect(DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged(:,1),Domain_mapping_Prot_Sim5_noSim3_MulFil_merged(:,1));
to1=0;to2=0;toex=0;
for i=1:length(IntSec)
    realind1=find(ismember(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged(:,1),IntSec(i))==1);
    Dom_temp=Domain_mapping_Prot_Sim5_noSim3_MulFil_merged(realind1,:);
    Dom_temp_on=Dom_temp(:,2);
    realind2=find(ismember(DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged(:,1),IntSec(i))==1);
    DomPair_temp=DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged(realind2,:);
    clear DomPair_temp_on
    for j=1:size(DomPair_temp,1) 
        DomPair_temp_on=strsplit(DomPair_temp{j,2},'-')';
        [Lia,~]=ismember(Dom_temp_on,DomPair_temp_on);
        Dom_temp_max_F1=max(cell2mat(Dom_temp(Lia==1,(end-1))));
        Dom_temp_max_MCC=max(cell2mat(Dom_temp(Lia==1,end)));
        if sum(Lia)>0
            if cell2mat(DomPair_temp(j,end))>Dom_temp_max_MCC && cell2mat(DomPair_temp(j,(end-1)))>Dom_temp_max_F1
                del_realind1=[del_realind1;realind1(Lia==1)];
                Examples_DomPair_better=[Examples_DomPair_better;Dom_temp(Lia==1,:)];
                Examples_DomPair_better=[Examples_DomPair_better;DomPair_temp(j,:)];
            else
                del_realind2=[del_realind2;realind2(j)];
            end
        end
    end
end
del_realind1=unique(del_realind1);
del_realind2=unique(del_realind2);
Domain_mapping_Prot_Sim5_noSim3_MulFil_merged_final=Domain_mapping_Prot_Sim5_noSim3_MulFil_merged;
Domain_mapping_Prot_Sim5_noSim3_MulFil_merged_final(del_realind1,:)=[];
save DRUIDom_Files/Domain_mapping_Prot_Sim5_noSim3_MulFil_merged_final.mat Domain_mapping_Prot_Sim5_noSim3_MulFil_merged_final
length(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged_final)
length(unique(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged_final(:,1)))
length(unique(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged_final(:,2)))
% (final single domain association stats: 27,021 mappings between 8,164 compounds and 250 InterPro domains)

DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged_final=DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged;
DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged_final(del_realind2,:)=[];
save DRUIDom_Files/DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged_final.mat DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged_final
length(DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged_final)
length(unique(DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged_final(:,1)))
length(unique(DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged_final(:,2)))
% (final domain pair association stats: 22 mappings between 10 compounds and 12 InterPro domain pairs)
t=DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged_final';
fid=fopen('DRUIDom_finalized_domainpair_compound_mappings.txt','w');
fprintf(fid,'%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n',t{:});
fclose(fid);


save DRUIDom_Files/Examples_DomPair_better.mat Examples_DomPair_better





% ---------------------------------------------------------
% Predicting DTIs based on the compound-domain associations
% ---------------------------------------------------------


% (merging the compound-domain association file with the protein domain annotation file by eliminating duplicates and by selecting the pairs with the maximum performance)
load DRUIDom_Files/InterPro_entriesmatched_allhuman_Org_onlydom.mat
load DRUIDom_Files/Domain_mapping_Prot_Sim5_noSim3_MulFil_merged_final.mat
h=cell(0,3);
Dom_map_unique=unique(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged_final(:,2));
for i=1:length(Dom_map_unique)
    disp(['Domain association number: ' num2str(i), ' / ' num2str(length(Dom_map_unique))])
    [Lia,~]=ismember(InterPro_entriesmatched_allhuman_Org_onlydom(:,2),Dom_map_unique(i,1));
    prot_temp=InterPro_entriesmatched_allhuman_Org_onlydom(Lia==1,1);
    [Lia2,Locb2]=ismember(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged_final(:,2),Dom_map_unique(i,1));
    comp_temp=Domain_mapping_Prot_Sim5_noSim3_MulFil_merged_final(Lia2==1,1);
    score_temp=Domain_mapping_Prot_Sim5_noSim3_MulFil_merged_final(Lia2==1,end);
    [A,B]=meshgrid(prot_temp,comp_temp);
    c=cat(2,A',B');
    d=reshape(c,[],2);
    [C,D]=meshgrid(prot_temp,score_temp);
    e=cat(2,C',D');
    f=reshape(e,[],2);
    d=[d(:,[2 1]) f(:,2)];
    if isempty(d)==0 && isempty(h)==0
        h_23=strcat(h(:,1), {' '}, h(:,2));
        d_23=strcat(d(:,1), {' '}, d(:,2));
        [~,I1,I2]=intersect(h_23,d_23);
        C1=(cell2mat(h(I1,end))>cell2mat(d(I2,end)));
        mer1=h(I1(C1==1),:);
        mer2=d(I2(C1==0),:);
        I1inv=setdiff((1:size(h,1)),I1)';
        I2inv=setdiff((1:size(d,1)),I2)';
        h=[h(I1inv,:);d(I2inv,:);mer1;mer2];
    elseif isempty(d)==0 && isempty(h)==1
        h=d;
    end
end
DomMap_DTI_Predictions_finalized=h;
save DRUIDom_Files/DomMap_DTI_Predictions_finalized.mat DomMap_DTI_Predictions_finalized
length(unique(DomMap_DTI_Predictions_finalized(:,1)))
length(unique(DomMap_DTI_Predictions_finalized(:,2)))
t=DomMap_DTI_Predictions_finalized';
fid=fopen('DRUIDom_Files/DomMap_DTI_Predictions_finalized.txt','w');
fprintf(fid,'%s\t%s\t%.4f\n',t{:});
fclose(fid);

% (extracting novel DTIs from the finalized prediction set)
load DRUIDom_Files/Combined_Org_act_var.mat
load DRUIDom_Files/DomMap_DTI_Predictions_finalized.mat
Combined_Org_act_23=strcat(Combined_Org_act(:,2), {' '}, Combined_Org_act(:,1));
DomMap_DTI_Predictions_finalized_23=strcat(DomMap_DTI_Predictions_finalized(:,1), {' '}, DomMap_DTI_Predictions_finalized(:,2));
[Lia,~]=ismember(DomMap_DTI_Predictions_finalized_23,Combined_Org_act_23);
sum(Lia)
% (40,442 predictions were shared with the actives in the training set)
load DRUIDom_Files/Combined_Org_inact_withact_var.mat
Combined_Org_inact_withact_23=strcat(Combined_Org_inact_withact(:,2), {' '}, Combined_Org_inact_withact(:,1));
[Lia,~]=ismember(DomMap_DTI_Predictions_finalized_23,Combined_Org_inact_withact_23);
sum(Lia)
% (6,887 predictions were shared and contradicted with the inactives in the training set)

[Lia,~]=ismember(DomMap_DTI_Predictions_finalized_23,Combined_Org_act_23);
DomMap_DTI_Predictions_finalized_novelpred=DomMap_DTI_Predictions_finalized(Lia==0,:);
save DRUIDom_Files/DomMap_DTI_Predictions_finalized_novelpred.mat DomMap_DTI_Predictions_finalized_novelpred
t=DomMap_DTI_Predictions_finalized_novelpred';
fid=fopen('DRUIDom_Files/DomMap_DTI_Predictions_finalized_novelpred.txt','w');
fprintf(fid,'%s\t%s\t%.4f\n',t{:});
fclose(fid);
length(DomMap_DTI_Predictions_finalized_novelpred)
length(unique(DomMap_DTI_Predictions_finalized_novelpred(:,1)))
length(unique(DomMap_DTI_Predictions_finalized_novelpred(:,2)))
% (single domain mapping based final prediction stats: 3,672,076 novel DTIs between 8,158 compounds and 5,563 proteins)



% Predicting DTIs based on the compound-domain pair associations:

% (merging the compound-domain association file with the protein domain annotation file by eliminating duplicates and by selecting the pairs with the maximum performance)
load DRUIDom_Files/InterPro_DA_based_entriesmatched_allhuman_Org_onlydom.mat
load DRUIDom_Files/DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged_final.mat
h=cell(0,3);
DomPair_map_unique=unique(DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged_final(:,2));
for i=1:length(DomPair_map_unique)
    disp(['Domain association number: ' num2str(i), ' / ' num2str(length(DomPair_map_unique))])
    [Lia,~]=ismember(InterPro_entriesmatched_allhuman_Org_bagofdom_onlydompairs(:,2),DomPair_map_unique(i,1));
    prot_temp=InterPro_entriesmatched_allhuman_Org_bagofdom_onlydompairs(Lia==1,1);
    [Lia2,Locb2]=ismember(DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged_final(:,2),DomPair_map_unique(i,1));
    comp_temp=DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged_final(Lia2==1,1);
    score_temp=DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged_final(Lia2==1,end);
    [A,B]=meshgrid(prot_temp,comp_temp);
    c=cat(2,A',B');
    d=reshape(c,[],2);
    [C,D]=meshgrid(prot_temp,score_temp);
    e=cat(2,C',D');
    f=reshape(e,[],2);
    d=[d(:,[2 1]) f(:,2)];
    if isempty(d)==0 && isempty(h)==0
        h_23=strcat(h(:,1), {' '}, h(:,2));
        d_23=strcat(d(:,1), {' '}, d(:,2));
        [~,I1,I2]=intersect(h_23,d_23);
        C1=(cell2mat(h(I1,end))>cell2mat(d(I2,end)));
        mer1=h(I1(C1==1),:);
        mer2=d(I2(C1==0),:);
        I1inv=setdiff((1:size(h,1)),I1)';
        I2inv=setdiff((1:size(d,1)),I2)';
        h=[h(I1inv,:);d(I2inv,:);mer1;mer2];
    elseif isempty(d)==0 && isempty(h)==1
        h=d;
    end
end
DomMap_DomPair_DTI_Predictions_finalized=h;
save DRUIDom_Files/DomMap_DomPair_DTI_Predictions_finalized.mat DomMap_DomPair_DTI_Predictions_finalized
length(unique(DomMap_DomPair_DTI_Predictions_finalized(:,1)))
length(unique(DomMap_DomPair_DTI_Predictions_finalized(:,2)))
t=DomMap_DomPair_DTI_Predictions_finalized';
fid=fopen('DRUIDom_Files/DomMap_DomPair_DTI_Predictions_finalized.txt','w');
fprintf(fid,'%s\t%s\t%.4f\n',t{:});
fclose(fid);

% (extracting novel DTIs from the finalized prediction set)
load DRUIDom_Files/Combined_Org_act_var.mat
load DRUIDom_Files/DomMap_DomPair_DTI_Predictions_finalized.mat
Combined_Org_act_23=strcat(Combined_Org_act(:,2), {' '}, Combined_Org_act(:,1));
DomMap_DomPair_DTI_Predictions_finalized_23=strcat(DomMap_DomPair_DTI_Predictions_finalized(:,1), {' '}, DomMap_DomPair_DTI_Predictions_finalized(:,2));
[Lia,~]=ismember(DomMap_DomPair_DTI_Predictions_finalized_23,Combined_Org_act_23);
sum(Lia)
% (27 predictions were shared with the actives in the training set)
load DRUIDom_Files/Combined_Org_inact_withact_var.mat
Combined_Org_inact_withact_23=strcat(Combined_Org_inact_withact(:,2), {' '}, Combined_Org_inact_withact(:,1));
[Lia,~]=ismember(DomMap_DomPair_DTI_Predictions_finalized_23,Combined_Org_inact_withact_23);
sum(Lia)
% (1 predictions were shared and contradicted with the inactives in the training set)

[Lia,~]=ismember(DomMap_DomPair_DTI_Predictions_finalized_23,Combined_Org_act_23);
DomMap_DomPair_DTI_Predictions_finalized_novelpred=DomMap_DomPair_DTI_Predictions_finalized(Lia==0,:);
save DRUIDom_Files/DomMap_DomPair_DTI_Predictions_finalized_novelpred.mat DomMap_DomPair_DTI_Predictions_finalized_novelpred
t=DomMap_DomPair_DTI_Predictions_finalized_novelpred';
fid=fopen('DRUIDom_Files/DomMap_DomPair_DTI_Predictions_finalized_novelpred.txt','w');
fprintf(fid,'%s\t%s\t%.4f\n',t{:});
fclose(fid);
length(DomMap_DomPair_DTI_Predictions_finalized_novelpred)
length(unique(DomMap_DomPair_DTI_Predictions_finalized_novelpred(:,1)))
length(unique(DomMap_DomPair_DTI_Predictions_finalized_novelpred(:,2)))
% (domain pair mapping based final prediction stats: 631 novel DTIs between 9 compounds and 286 proteins)



% Merging the single domain and domain pair association based novel DTI prediction files:

load DRUIDom_Files/DomMap_DTI_Predictions_finalized_novelpred.mat
load DRUIDom_Files/DomMap_DomPair_DTI_Predictions_finalized_novelpred.mat
DomMap_DTI_Predictions_finalized_novelpred_23=strcat(DomMap_DTI_Predictions_finalized_novelpred(:,1), {' '}, DomMap_DTI_Predictions_finalized_novelpred(:,2));
DomMap_DomPair_DTI_Predictions_finalized_novelpred_23=strcat(DomMap_DomPair_DTI_Predictions_finalized_novelpred(:,1), {' '}, DomMap_DomPair_DTI_Predictions_finalized_novelpred(:,2));
[ins,ia,ib]=intersect(DomMap_DTI_Predictions_finalized_novelpred_23,DomMap_DomPair_DTI_Predictions_finalized_novelpred_23);
find(cell2mat(DomMap_DTI_Predictions_finalized_novelpred(ia,3))>cell2mat(DomMap_DomPair_DTI_Predictions_finalized_novelpred(ib,3)))
% (none of the single domain based preditions scored better compared to domain pair mappings -most of them had the exact same score-, so all intersections can be selected from the domain pair association based predictions)
% (there is one example where the predictions score was improved due to domain pair mapping: compound id: CHEMBL450519, protein acc: P06239, score of singe domain and domain pair mappings: 0.7746 & 0.8433)
DomMap_DTI_Predictions_finalized_novelpred(ia,:)=[];
DomMap_SingDom_DomPair_merged_DTI_Pred_final_novelpred=[DomMap_DTI_Predictions_finalized_novelpred;DomMap_DomPair_DTI_Predictions_finalized_novelpred];
DomMap_SingDom_DomPair_merged_DTI_Pred_final_novelpred=sortrows(DomMap_SingDom_DomPair_merged_DTI_Pred_final_novelpred);
save DRUIDom_Files/DomMap_SingDom_DomPair_merged_DTI_Pred_final_novelpred.mat DomMap_SingDom_DomPair_merged_DTI_Pred_final_novelpred
t=DomMap_SingDom_DomPair_merged_DTI_Pred_final_novelpred';
fid=fopen('DRUIDom_Files/DomMap_SingDom_DomPair_merged_DTI_Pred_final_novelpred.txt','w');
fprintf(fid,'%s\t%s\t%.4f\n',t{:});
fclose(fid);
length(DomMap_SingDom_DomPair_merged_DTI_Pred_final_novelpred)
length(unique(DomMap_SingDom_DomPair_merged_DTI_Pred_final_novelpred(:,1)))
length(unique(DomMap_SingDom_DomPair_merged_DTI_Pred_final_novelpred(:,2)))
% (finalized merged prediction stats: 3,672,220 novel DTIs between 8,163 compounds and 5,563 proteins)


% Propagating the predictions to other compounds based on molecular similarity:

% (similarity threshold is selected as 0.8)

load DRUIDom_Files/Combined_Org_Sim_07.mat
load DRUIDom_Files/Combined_Org_Comp1_07.mat
load DRUIDom_Files/Combined_Org_Comp2_07.mat
ind=find(Combined_Org_Sim>=0.8);
Combined_Org_Comp1_08=Combined_Org_Comp1(ind);
Combined_Org_Comp2_08=Combined_Org_Comp2(ind);
Combined_Org_Sim_08=Combined_Org_Sim(ind);
clear Combined_Org_Comp1 Combined_Org_Comp2 Combined_Org_Sim
save DRUIDom_Files/Combined_Org_Comp1_08.mat Combined_Org_Comp1_08
save DRUIDom_Files/Combined_Org_Comp2_08.mat Combined_Org_Comp2_08
save DRUIDom_Files/Combined_Org_Sim_08.mat Combined_Org_Sim_08


load DRUIDom_Files/DomMap_SingDom_DomPair_merged_DTI_Pred_final_novelpred.mat
pred_scores_cell=DomMap_SingDom_DomPair_merged_DTI_Pred_final_novelpred(:,3);
load DRUIDom_Files/Combined_Org_Sim_08.mat
load DRUIDom_Files/Combined_Org_Comp1_08.mat
load DRUIDom_Files/Combined_Org_Comp2_08.mat
Combined_Org_Sim_08_cell=num2cell(Combined_Org_Sim_08);
pred_comp_unique=unique(DomMap_SingDom_DomPair_merged_DTI_Pred_final_novelpred(:,1));
h=cell(round((length(Combined_Org_Comp1_08)/length(unique(Combined_Org_Comp1_08)))*(length(DomMap_SingDom_DomPair_merged_DTI_Pred_final_novelpred)/length(unique(DomMap_SingDom_DomPair_merged_DTI_Pred_final_novelpred(:,1))))*(length(pred_comp_unique))*1.5),3);
to=1;
for i=1:length(pred_comp_unique)
    disp(['Compound number: ' num2str(i), ' / ' num2str(length(pred_comp_unique))])
    cent_comp=pred_comp_unique(i,1);
    
    ind1=find(ismember(DomMap_SingDom_DomPair_merged_DTI_Pred_final_novelpred(:,1),cent_comp)==1);
    prot_temp=DomMap_SingDom_DomPair_merged_DTI_Pred_final_novelpred(ind1,2);
    prot_score_temp=pred_scores_cell(ind1,1);
    
    comp_temp1=Combined_Org_Comp2_08(ismember(Combined_Org_Comp1_08,cent_comp)==1,1);
    comp_score_temp1=Combined_Org_Sim_08_cell(ismember(Combined_Org_Comp1_08,cent_comp)==1,1);
    comp_temp2=Combined_Org_Comp1_08(ismember(Combined_Org_Comp2_08,cent_comp)==1,1);
    comp_score_temp2=Combined_Org_Sim_08_cell(ismember(Combined_Org_Comp2_08,cent_comp)==1,1);
    comp_temp=[comp_temp1;comp_temp2];
    comp_score_temp=[comp_score_temp1;comp_score_temp2];
    
    if isempty(comp_temp)==0
        [A,B]=meshgrid(prot_temp,comp_temp);
        c=cat(2,A',B');
        d=reshape(c,[],2);
        [C,D]=meshgrid(prot_temp,comp_score_temp);
        e=cat(2,C',D');
        f=reshape(e,[],2);
        [E,F]=meshgrid(prot_score_temp,comp_temp);
        j=cat(2,F',E');
        k=reshape(j,[],2); 
        d=[d(:,[2 1]) num2cell(cell2mat(f(:,2)).*cell2mat(k(:,2)))];
        
        h(to:(to+size(d,1)-1),1:3)=d;
        to=to+size(d,1);
%         if isempty(d)==0 && isempty(h)==0
%             h_23=strcat(h(:,1), {' '}, h(:,2));
%             d_23=strcat(d(:,1), {' '}, d(:,2));
%             [~,I1,I2]=intersect(h_23,d_23);
%             C1=(cell2mat(h(I1,end))>cell2mat(d(I2,end)));
%             mer1=h(I1(C1==1),:);
%             mer2=d(I2(C1==0),:);
%             I1inv=setdiff((1:size(h,1)),I1)';
%             I2inv=setdiff((1:size(d,1)),I2)';
%             h=[h(I1inv,:);d(I2inv,:);mer1;mer2];
%         elseif isempty(d)==0 && isempty(h)==1
%             h=d;
%         end
    end  
end
h_23=strcat(h(:,1), {' '}, h(:,2));
h_23_unique=unique(h_23);
[~,Locb]=ismember(h_23,h_23_unique);
h_23_unique_maxscr=zeros(length(h_23_unique),1);
h_23_unique_sep=split(h_23_unique,' ');
scr_all=cell2mat(h(:,3));
for i=3776363:5082192
    disp(['Unique pair number: ' num2str(i), ' / ' num2str(length(h_23_unique_maxscr))])
    h_23_unique_maxscr(i,1)=max(scr_all(Locb==i,1));
end
DomMap_SingDom_DomPair_merged_DTI_Pred_propagat=[h_23_unique_sep num2cell(h_23_unique_maxscr)];

load DRUIDom_Files/Combined_Org_act_var.mat
Combined_Org_act_23=strcat(Combined_Org_act(:,2), {' '}, Combined_Org_act(:,1));
DomMap_SingDom_DomPair_merged_DTI_Pred_propagat_23=strcat(DomMap_SingDom_DomPair_merged_DTI_Pred_propagat(:,1), {' '}, DomMap_SingDom_DomPair_merged_DTI_Pred_propagat(:,2));
[Lia,~]=ismember(DomMap_SingDom_DomPair_merged_DTI_Pred_propagat_23,Combined_Org_act_23);
sum(Lia)
% (26,533 predictions were shared with the actives in the training set)
DomMap_SingDom_DomPair_merged_DTI_Pred_propagat=DomMap_SingDom_DomPair_merged_DTI_Pred_propagat(Lia==0,:);

load DRUIDom_Files/Combined_Org_inact_withact_var.mat
Combined_Org_inact_withact_23=strcat(Combined_Org_inact_withact(:,2), {' '}, Combined_Org_inact_withact(:,1));
DomMap_SingDom_DomPair_merged_DTI_Pred_propagat_23=strcat(DomMap_SingDom_DomPair_merged_DTI_Pred_propagat(:,1), {' '}, DomMap_SingDom_DomPair_merged_DTI_Pred_propagat(:,2));
[Lia,~]=ismember(DomMap_SingDom_DomPair_merged_DTI_Pred_propagat_23,Combined_Org_inact_withact_23);
sum(Lia)
% (4,818 predictions were shared and contradicted with the inactives in the training set)
DomMap_SingDom_DomPair_merged_DTI_Pred_propagat_novelpred_final=DomMap_SingDom_DomPair_merged_DTI_Pred_propagat(Lia==0,:);

save DRUIDom_Files/DomMap_SingDom_DomPair_merged_DTI_Pred_propagat_novelpred_final.mat DomMap_SingDom_DomPair_merged_DTI_Pred_propagat_novelpred_final
length(unique(DomMap_SingDom_DomPair_merged_DTI_Pred_propagat_novelpred_final(:,1)))
length(unique(DomMap_SingDom_DomPair_merged_DTI_Pred_propagat_novelpred_final(:,2)))
t=DomMap_SingDom_DomPair_merged_DTI_Pred_propagat_novelpred_final';
fid=fopen('DRUIDom_Files/DomMap_SingDom_DomPair_merged_DTI_Pred_propagat_novelpred_final.txt','w');
fprintf(fid,'%s\t%s\t%.4f\n',t{:});
fclose(fid);
% (finalized propagated merged prediction stats: 5,050,841 novel DTIs between 10,944 compounds and 5,461 proteins)
fid=fopen('DRUIDom_finalized_novel_DTI_predictions.txt','w');
fprintf(fid,'%s\t%s\t%.4f\n',t{:});
fclose(fid);





% --------------------------------------------------------------
% Data Analysis Steps for Performance Calculation and Comparison
% --------------------------------------------------------------


% Domain mapping performance analysis using structurally known binding regions
% of proteins from PDB using the mappings of the InteracDome study (folder: /Performance_analysis):

[MappingPI_Pfam,MappingPI_InterPro,~]=textread('DRUIDom_Files/Performance_analysis/pfam_v32_interpro_mapping.txt', '%s %s %s', 'delimiter', '\t', 'bufsize', 250000);
[MappingCP_ChEMBLid,MappingCP_PDBligandid]=textread('DRUIDom_Files/Performance_analysis/ChEMBL_PDBligand_mapping.txt', '%s %s', 'delimiter', '\t', 'headerlines', 1);
[MappingPP_PDBligandid,MappingPP_PubChemid]=textread('DRUIDom_Files/Performance_analysis/PubChem_PDBligand_mapping.txt', '%s %s', 'delimiter', '\t', 'headerlines', 1);

[InteracDome_repNR_Pfamid,~,InteracDome_repNR_PDBligid,~,~,~,~,~,~,~,~,~]=textread('DRUIDom_Files/Performance_analysis/InteracDome_v0.3-representableNR.tsv', '%s %s %s %s %s %s %s %s %s %s %s %s', 'delimiter', '\t', 'headerlines', 5, 'bufsize', 250000);
InteracDome_repNR_Pfamid=cellfun(@(x) x(1:7), InteracDome_repNR_Pfamid, 'un', 0);
InteracDome_repNR_pairs=unique(strcat(InteracDome_repNR_Pfamid, {' '}, InteracDome_repNR_PDBligid));
[InteracDome_rep_Pfamid,~,InteracDome_rep_PDBligid,~,~,~,~,~,~,~,~,~]=textread('DRUIDom_Files/Performance_analysis/InteracDome_v0.3-representable.tsv', '%s %s %s %s %s %s %s %s %s %s %s %s', 'delimiter', '\t', 'headerlines', 5, 'bufsize', 250000);
InteracDome_rep_Pfamid=cellfun(@(x) x(1:7), InteracDome_rep_Pfamid, 'un', 0);
InteracDome_rep_pairs=unique(strcat(InteracDome_rep_Pfamid, {' '}, InteracDome_rep_PDBligid));
InteracDome_rep_pairs_diff=setdiff(InteracDome_rep_pairs,InteracDome_repNR_pairs);

[InterPro_v72_Domain_id,~,~]=textread('DRUIDom_Files/Performance_analysis/InterPro_v72_Domain_entrylist.txt', '%s %s %s', 'delimiter', '\t');
[PDB_druglike_compounds]=textread('DRUIDom_Files/Performance_analysis/PDB_druglike_compounds_ids.txt', '%s', 'delimiter', '\t');

% (calculating the maximum possible coverage in terms of domains and ligands)
InterPro_v72_Domain_id=unique(InterPro_v72_Domain_id);
[Lia,~]=ismember(MappingPI_InterPro,InterPro_v72_Domain_id);
MappingPI_InterPro_domain=MappingPI_InterPro(Lia==1);
MappingPI_Pfam_domain=MappingPI_Pfam(Lia==1);
InteracDome_repNR_Pfamid_all_unique=unique(InteracDome_repNR_Pfamid);
[Lia,~]=ismember(InteracDome_repNR_Pfamid_all_unique,MappingPI_Pfam_domain);
InteracDome_repNR_Pfamid_InterProDomainMap_unique=InteracDome_repNR_Pfamid_all_unique(Lia==1);
[Lia,~]=ismember(InteracDome_repNR_Pfamid,InteracDome_repNR_Pfamid_InterProDomainMap_unique);
InteracDome_repNR_Pfamid_InterProDomainMap=InteracDome_repNR_Pfamid(Lia==1);
InteracDome_repNR_PDBligid_InterProDomainMap=InteracDome_repNR_PDBligid(Lia==1);
MappingPP_CP_PDBligandid=unique([MappingPP_PDBligandid;MappingCP_PDBligandid]);
[Lia,~]=ismember(InteracDome_repNR_PDBligid_InterProDomainMap,MappingPP_CP_PDBligandid);
InteracDome_repNR_PDBligid_InterProDomainMap_sharedligand=InteracDome_repNR_PDBligid_InterProDomainMap(Lia==1);
InteracDome_repNR_Pfamid_InterProDomainMap_sharedligand=InteracDome_repNR_Pfamid_InterProDomainMap(Lia==1);
Max_Cov_dom=length(unique(InteracDome_repNR_Pfamid_InterProDomainMap_sharedligand));
Max_Cov_lig=length(unique(InteracDome_repNR_PDBligid_InterProDomainMap_sharedligand));

% (calculating the performance for the recovery of InteracDome pairs, together with coverage)
load DRUIDom_Files/Domain_mapping_Prot_5.mat
load DRUIDom_Files/Domain_mapping_Prot_noSim_3.mat

ConfThresArray_REC=0:0.1:1;
ConfThresArray_PRE=0:0.1:1;
ConfThresArray_ACCU=0:0.1:1;
ConfThresArray_F1=0:0.1:1;
ConfThresArray_MCC=0:0.1:1;
Perf_results_DomainMapping=zeros(length(ConfThresArray_ACCU)*length(ConfThresArray_MCC),21);
too=0;
%(previous finalized mapping parameters were: MCC>=0.5, F1>=0.7, ACCU>=0.8, REC>=0.7, PRE>=0.7)
for j=1:length(ConfThresArray_MCC)
    for k=1:length(ConfThresArray_ACCU)
        too=too+1;
        disp(['Threshold array index: ' num2str(too), ' / ' num2str(length(ConfThresArray_ACCU)*length(ConfThresArray_MCC))])
        % (generating the compound similarity based mappings according to the selected confidence threshold)
        Domain_mapping_Prot_5_all_MulFil=cell(4000000,9);
        to=1;
        for i=1:length(Domain_mapping_Prot_5.Compounds)
            disp(['Compound number: ' num2str(i), ' / ' num2str(length(Domain_mapping_Prot_5.Compounds))])
            ind_temp=find(Domain_mapping_Prot_5.MapScr{i,1}(:,5)>=ConfThresArray_REC(1,k) & Domain_mapping_Prot_5.MapScr{i,1}(:,6)>=ConfThresArray_PRE(1,k) & Domain_mapping_Prot_5.MapScr{i,1}(:,7)>=ConfThresArray_ACCU(1,k) & Domain_mapping_Prot_5.MapScr{i,1}(:,8)>=ConfThresArray_F1(1,k) & Domain_mapping_Prot_5.MapScr{i,1}(:,9)>=ConfThresArray_MCC(1,j));
            if isempty(ind_temp)==0
                Domain_mapping_Prot_5_all_MulFil(to:(to+(length(ind_temp))-1),1)=cellstr(repmat(Domain_mapping_Prot_5.Compounds{i,1},length(ind_temp),1));
                Domain_mapping_Prot_5_all_MulFil(to:(to+(length(ind_temp))-1),2)=Domain_mapping_Prot_5.Domains{i,1}(ind_temp,1);
                Domain_mapping_Prot_5_all_MulFil(to:(to+(length(ind_temp))-1),3:11)=num2cell(Domain_mapping_Prot_5.MapScr{i,1}(ind_temp,1:9));
                to=to+(length(ind_temp));
            end
        end
        del_ind=cellfun(@isempty,Domain_mapping_Prot_5_all_MulFil(:,1));
        Domain_mapping_Prot_5_all_MulFil(del_ind,:)=[];
        % (generating the compound non-similarity based mappings according to the selected confidence threshold)
        Domain_mapping_Prot_noSim_3_all_MulFil=cell(4000000,9);
        to=1;
        for i=1:length(Domain_mapping_Prot_noSim_3.Compounds)
            disp(['Compound number: ' num2str(i), ' / ' num2str(length(Domain_mapping_Prot_noSim_3.Compounds))])
            if isempty(Domain_mapping_Prot_noSim_3.MapScr{i,1})==0
                ind_temp=find(Domain_mapping_Prot_noSim_3.MapScr{i,1}(:,5)>=ConfThresArray_REC(1,k) & Domain_mapping_Prot_noSim_3.MapScr{i,1}(:,6)>=ConfThresArray_PRE(1,k) & Domain_mapping_Prot_noSim_3.MapScr{i,1}(:,7)>=ConfThresArray_ACCU(1,k) & Domain_mapping_Prot_noSim_3.MapScr{i,1}(:,8)>=ConfThresArray_F1(1,k) & Domain_mapping_Prot_noSim_3.MapScr{i,1}(:,9)>=ConfThresArray_MCC(1,j));
            end
            if isempty(ind_temp)==0
                Domain_mapping_Prot_noSim_3_all_MulFil(to:(to+(length(ind_temp))-1),1)=cellstr(repmat(Domain_mapping_Prot_noSim_3.Compounds{i,1},length(ind_temp),1));
                Domain_mapping_Prot_noSim_3_all_MulFil(to:(to+(length(ind_temp))-1),2)=Domain_mapping_Prot_noSim_3.Domains{i,1}(ind_temp,1);
                Domain_mapping_Prot_noSim_3_all_MulFil(to:(to+(length(ind_temp))-1),3:11)=num2cell(Domain_mapping_Prot_noSim_3.MapScr{i,1}(ind_temp,1:9));
                to=to+(length(ind_temp));
            end
            ind_temp=[];
        end
        del_ind=cellfun(@isempty,Domain_mapping_Prot_noSim_3_all_MulFil(:,1));
        Domain_mapping_Prot_noSim_3_all_MulFil(del_ind,:)=[];
        % (merging the similarity and non-similarity based mappings)
        Domain_mapping_Prot_5_all_MulFil_23=strcat(Domain_mapping_Prot_5_all_MulFil(:,1), {' '}, Domain_mapping_Prot_5_all_MulFil(:,2));
        Domain_mapping_Prot_noSim_3_all_MulFil_23=strcat(Domain_mapping_Prot_noSim_3_all_MulFil(:,1), {' '}, Domain_mapping_Prot_noSim_3_all_MulFil(:,2));
        [IntSec,I1,I2]=intersect(Domain_mapping_Prot_5_all_MulFil_23,Domain_mapping_Prot_noSim_3_all_MulFil_23);
        C1=(cell2mat(Domain_mapping_Prot_5_all_MulFil(I1,end))>cell2mat(Domain_mapping_Prot_noSim_3_all_MulFil(I2,end)));
        mer1=Domain_mapping_Prot_5_all_MulFil(I1(C1==1),:);
        mer2=Domain_mapping_Prot_noSim_3_all_MulFil(I2(C1==0),:);
        I1inv=setdiff((1:length(Domain_mapping_Prot_5_all_MulFil)),I1)';
        I2inv=setdiff((1:length(Domain_mapping_Prot_noSim_3_all_MulFil)),I2)';
        Domain_mapping_Prot_Sim5_noSim3_MulFil_merged=[Domain_mapping_Prot_5_all_MulFil(I1inv,:);Domain_mapping_Prot_noSim_3_all_MulFil(I2inv,:);mer1;mer2];
        
        [Lia,Locb]=ismember(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged(:,1),MappingPP_PubChemid);
        Domain_mapping_PDB_PDBligandid1=MappingPP_PDBligandid(Locb(Locb>0));
        Domain_mapping_PDB_InterPro1=Domain_mapping_Prot_Sim5_noSim3_MulFil_merged(Lia==1,2);
        [Lia,Locb]=ismember(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged(:,1),MappingCP_ChEMBLid);
        Domain_mapping_PDB_PDBligandid2=MappingCP_PDBligandid(Locb(Locb>0));
        Domain_mapping_PDB_InterPro2=Domain_mapping_Prot_Sim5_noSim3_MulFil_merged(Lia==1,2);
        Domain_mapping_PDB_PDBligandid=[Domain_mapping_PDB_PDBligandid1;Domain_mapping_PDB_PDBligandid2];
        Domain_mapping_PDB_InterPro=[Domain_mapping_PDB_InterPro1;Domain_mapping_PDB_InterPro2];
        [Lia,Locb]=ismember(Domain_mapping_PDB_InterPro,MappingPI_InterPro);
        Domain_mapping_PDB_Pfam_Pfamid=MappingPI_Pfam(Locb(Locb>0));
        Domain_mapping_PDB_Pfam_PDBligandid=Domain_mapping_PDB_PDBligandid(Lia==1);
        %(above 2 variables are all Pfam-PDBligand mappings created with the methodology proposed in this project)
        
        [Lia,~]=ismember(Domain_mapping_PDB_Pfam_PDBligandid,InteracDome_repNR_PDBligid);
        Mapping_InteracDome_repNR_PDBligandid=Domain_mapping_PDB_Pfam_PDBligandid(Lia==1);
        Mapping_InteracDome_repNR_Pfamid=Domain_mapping_PDB_Pfam_Pfamid(Lia==1);
        [Lia,~]=ismember(Mapping_InteracDome_repNR_Pfamid,InteracDome_repNR_Pfamid);
        Mapping_InteracDome_repNR_final_Pfamid=Mapping_InteracDome_repNR_Pfamid(Lia==1);
        Mapping_InteracDome_repNR_final_PDBligandid=Mapping_InteracDome_repNR_PDBligandid(Lia==1);
        %(above 2 variables are the Pfam-PDBligand mappings, where Pfams and PDBligands intersect with the InteracDome Pfams and PDBligands)
        
        [Lia,~]=ismember(InteracDome_repNR_Pfamid,Mapping_InteracDome_repNR_final_Pfamid);
        Shared_InteracDome_repNR_Pfamid=InteracDome_repNR_Pfamid(Lia==1);
        Shared_InteracDome_repNR_PDBligid=InteracDome_repNR_PDBligid(Lia==1);
        [Lia,~]=ismember(Shared_InteracDome_repNR_PDBligid,Mapping_InteracDome_repNR_final_PDBligandid);
        Shared2_InteracDome_repNR_Pfamid=Shared_InteracDome_repNR_Pfamid(Lia==1);
        Shared2_InteracDome_repNR_PDBligid=Shared_InteracDome_repNR_PDBligid(Lia==1);
        %(above 2 variables are the InteracDome mappings, where Pfams and PDBligands intersect with the Pfams and PDBligands in the mappings created in this study)
        
        %(numbers, coverage and extra coverage calculations)
        Num_map=length(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged);
        Num_dom=length(unique(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged(:,2)));
        Num_lig=length(unique(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged(:,1)));
        [Lia,~]=ismember(unique(InteracDome_repNR_Pfamid_InterProDomainMap_unique),Domain_mapping_PDB_Pfam_Pfamid);
        Cov_dom=length(find(Lia==1))/Max_Cov_dom;
        [Lia,~]=ismember(unique(InteracDome_repNR_PDBligid_InterProDomainMap_sharedligand),Domain_mapping_PDB_PDBligandid);
        Cov_lig=length(find(Lia==1))/Max_Cov_lig;
        ExCov_dom=(Num_dom/Max_Cov_dom)-Cov_dom;
        ExCov_lig=(Num_lig/Max_Cov_lig)-Cov_lig;
        %(performance calculation)
        InteracDome_perf_cal_pairs=unique(strcat(Shared2_InteracDome_repNR_Pfamid, {' '}, Shared2_InteracDome_repNR_PDBligid));
        DomMap_perf_cal_pairs=unique(strcat(Mapping_InteracDome_repNR_final_Pfamid, {' '}, Mapping_InteracDome_repNR_final_PDBligandid));
        TP=length(intersect(InteracDome_perf_cal_pairs,DomMap_perf_cal_pairs));
        [Lia,~]=ismember(DomMap_perf_cal_pairs,InteracDome_rep_pairs_diff);
        FP=length(DomMap_perf_cal_pairs)-TP-sum(Lia);
        FN=length(InteracDome_perf_cal_pairs)-TP;
        TN=length(unique(Shared2_InteracDome_repNR_PDBligid))*length(unique(Shared2_InteracDome_repNR_Pfamid))-TP-FN-FP-sum(Lia);
        MCC=((TP*TN)-(FP*FN))/(((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))^(1/2));
        ACCU=(TP+TN)/(TP+TN+FP+FN);
        F1=(2*TP)/(2*TP+FP+FN);
        REC=TP/(TP+FN);
        PRE=TP/(TP+FP);
        
        Perf_results_DomainMapping(too,1:21)=[ConfThresArray_REC(1,k) ConfThresArray_PRE(1,k) ConfThresArray_ACCU(1,k) ConfThresArray_F1(1,k) ConfThresArray_MCC(1,j) Num_map Num_dom Num_lig Cov_dom Cov_lig ExCov_dom ExCov_lig TP FP FN TN REC PRE ACCU F1 MCC];
    end
end
save DRUIDom_Files/Performance_analysis/Perf_results_DomainMapping.mat Perf_results_DomainMapping


% Calculation of the statistics for InteracDome:

length(InteracDome_repNR_Pfamid)
length(unique(InteracDome_repNR_Pfamid))
length(unique(InteracDome_repNR_PDBligid))

[Lia,~]=ismember(InteracDome_repNR_PDBligid,PDB_druglike_compounds);
InteracDome_repNR_PDBligid_druglikesmallmol=InteracDome_repNR_PDBligid(Lia==1);
InteracDome_repNR_Pfamid_druglikesmallmol=InteracDome_repNR_Pfamid(Lia==1);
length(InteracDome_repNR_Pfamid_druglikesmallmol)
length(unique(InteracDome_repNR_Pfamid_druglikesmallmol))
length(unique(InteracDome_repNR_PDBligid_druglikesmallmol))

save DRUIDom_Files/Performance_analysis/Performance_analysis_all_variables.mat


% Calculation of domain mapping statistics:

subs = findgroups(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged(:, 2));
dom_uniq = unique(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged(:, 2));
dom_freq_hist=hist(subs,0.5:1:250.5);
dom_freq_hist_sort=sort(dom_freq_hist,'descend');
sum(dom_freq_hist_sort(1,1:10))/sum(dom_freq_hist_sort)
dom_most_freq_10=dom_uniq(dom_freq_hist>1080);


% Calculation of protein & family statistics in the source bioactivty datasetand the output predictions dataset:

load DRUIDom_Files/Combined_Org_act_inact.mat
load DRUIDom_Files/Combined_Org_act_var.mat
load DRUIDom_Files/Combined_Org_inact_withact_var.mat
load DRUIDom_Files/DomMap_SingDom_DomPair_merged_DTI_Pred_final_novelpred.mat
Pred_prot_acc=DomMap_SingDom_DomPair_merged_DTI_Pred_final_novelpred(:,2);

length(unique(Combined_Org_act_inact(:,1)))
length(unique(Combined_Org_act(:,1)))
length(unique(Combined_Org_inact_withact(:,1)))
length(unique(Pred_prot_acc))

[Enzyme_UniProt_acc,~,~,~]=textread('DRUIDom_Files/Protein_family_datasets/Enzyme_uniprot_reviewed_human.tab', '%s %s %s %s', 'headerlines', 1, 'delimiter', '\t');
[GPCR_UniProt_acc,~,~]=textread('DRUIDom_Files/Protein_family_datasets/GPCR_uniprot_reviewed_human.tab', '%s %s %s', 'headerlines', 1, 'delimiter', '\t');
[IonChannel_UniProt_acc,~,~]=textread('DRUIDom_Files/Protein_family_datasets/IonChannel_uniprot_reviewed_human.tab', '%s %s %s', 'headerlines', 1, 'delimiter', '\t');
[TranscriptionFactor_UniProt_acc,~,~]=textread('DRUIDom_Files/Protein_family_datasets/TranscriptionFactor_uniprot_reviewed_human.tab', '%s %s %s', 'headerlines', 1, 'delimiter', '\t');
[OtherFamilies_UniProt_acc,~,~]=textread('DRUIDom_Files/Protein_family_datasets/OtherFamilies_uniprot_reviewed_human.tab', '%s %s %s', 'headerlines', 1, 'delimiter', '\t');

Freq_act=zeros(5,1);
Freq_uniq_prot_act=zeros(5,1);
[Lia,~]=ismember(Combined_Org_act(:,1),Enzyme_UniProt_acc);
Freq_act(1,1)=sum(Lia);
Freq_uniq_prot_act(1,1)=length(unique(Combined_Org_act(Lia,1)));
[Lia,~]=ismember(Combined_Org_act(:,1),GPCR_UniProt_acc);
Freq_act(2,1)=sum(Lia);
Freq_uniq_prot_act(2,1)=length(unique(Combined_Org_act(Lia,1)));
[Lia,~]=ismember(Combined_Org_act(:,1),IonChannel_UniProt_acc);
Freq_act(3,1)=sum(Lia);
Freq_uniq_prot_act(3,1)=length(unique(Combined_Org_act(Lia,1)));
[Lia,~]=ismember(Combined_Org_act(:,1),TranscriptionFactor_UniProt_acc);
Freq_act(4,1)=sum(Lia);
Freq_uniq_prot_act(4,1)=length(unique(Combined_Org_act(Lia,1)));
[Lia,~]=ismember(Combined_Org_act(:,1),OtherFamilies_UniProt_acc);
Freq_act(5,1)=sum(Lia);
Freq_uniq_prot_act(5,1)=length(unique(Combined_Org_act(Lia,1)));
for i=1:5
    Freq_act(i,2)=Freq_act(i,1)/sum(Freq_act(:,1));
end
for i=1:5
    Freq_uniq_prot_act(i,2)=Freq_uniq_prot_act(i,1)/sum(Freq_uniq_prot_act(:,1));
end

Freq_inact=zeros(5,1);
Freq_uniq_prot_inact=zeros(5,1);
[Lia,~]=ismember(Combined_Org_inact_withact(:,1),Enzyme_UniProt_acc);
Freq_inact(1,1)=sum(Lia);
Freq_uniq_prot_inact(1,1)=length(unique(Combined_Org_inact_withact(Lia,1)));
[Lia,~]=ismember(Combined_Org_inact_withact(:,1),GPCR_UniProt_acc);
Freq_inact(2,1)=sum(Lia);
Freq_uniq_prot_inact(2,1)=length(unique(Combined_Org_inact_withact(Lia,1)));
[Lia,~]=ismember(Combined_Org_inact_withact(:,1),IonChannel_UniProt_acc);
Freq_inact(3,1)=sum(Lia);
Freq_uniq_prot_inact(3,1)=length(unique(Combined_Org_inact_withact(Lia,1)));
[Lia,~]=ismember(Combined_Org_inact_withact(:,1),TranscriptionFactor_UniProt_acc);
Freq_inact(4,1)=sum(Lia);
Freq_uniq_prot_inact(4,1)=length(unique(Combined_Org_inact_withact(Lia,1)));
[Lia,~]=ismember(Combined_Org_inact_withact(:,1),OtherFamilies_UniProt_acc);
Freq_inact(5,1)=sum(Lia);
Freq_uniq_prot_inact(5,1)=length(unique(Combined_Org_inact_withact(Lia,1)));
for i=1:5
    Freq_inact(i,2)=Freq_inact(i,1)/sum(Freq_inact(:,1));
end
for i=1:5
    Freq_uniq_prot_inact(i,2)=Freq_uniq_prot_inact(i,1)/sum(Freq_uniq_prot_inact(:,1));
end

Freq_pred_act=zeros(5,1);
Freq_pred_uniq_prot_act=zeros(5,1);
[Lia,~]=ismember(Pred_prot_acc,Enzyme_UniProt_acc);
Freq_pred_act(1,1)=sum(Lia);
Freq_pred_uniq_prot_act(1,1)=length(unique(Pred_prot_acc(Lia,1)));
[Lia,~]=ismember(Pred_prot_acc,GPCR_UniProt_acc);
Freq_pred_act(2,1)=sum(Lia);
Freq_pred_uniq_prot_act(2,1)=length(unique(Pred_prot_acc(Lia,1)));
[Lia,~]=ismember(Pred_prot_acc,IonChannel_UniProt_acc);
Freq_pred_act(3,1)=sum(Lia);
Freq_pred_uniq_prot_act(3,1)=length(unique(Pred_prot_acc(Lia,1)));
[Lia,~]=ismember(Pred_prot_acc,TranscriptionFactor_UniProt_acc);
Freq_pred_act(4,1)=sum(Lia);
Freq_pred_uniq_prot_act(4,1)=length(unique(Pred_prot_acc(Lia,1)));
[Lia,Locb]=ismember(Pred_prot_acc,OtherFamilies_UniProt_acc);
Freq_pred_act(5,1)=sum(Lia);
Freq_pred_uniq_prot_act(5,1)=length(unique(Pred_prot_acc(Lia,1)));
for i=1:5
    Freq_pred_act(i,2)=Freq_pred_act(i,1)/sum(Freq_pred_act(:,1));
end
for i=1:5
    Freq_pred_uniq_prot_act(i,2)=Freq_pred_uniq_prot_act(i,1)/sum(Freq_pred_uniq_prot_act(:,1));
end


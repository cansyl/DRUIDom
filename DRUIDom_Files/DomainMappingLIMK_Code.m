% Dataset preparation process:

% Extracting active and inactive data points from the ExCAPE dataset and further filtering and organization operations:

% xz -d pubchem.chembl.dataset4publication_inchi_smiles_v2.tsv.xz
% cut -d$'\t' -f2,4,5,8,9,12 ExCAPE.pubchem.chembl.dataset4publication_inchi_smiles_v2.tsv > ExCAPE_all_datapoints.tsv
% cut -d$'\t' -f4,5 ExCAPE_all_datapoints.tsv > ExCAPE_all_datapoints_Taxonid_GeneName.tsv
% awk -F'\t' '$2~/N/' ExCAPE_all_datapoints.tsv > ExCAPE_inactive_datapoints.tsv
% awk -F'\t' '$2~/A/' ExCAPE_all_datapoints.tsv > ExCAPE_active_datapoints.tsv
% cut -d$'\t' -f4,5 ExCAPE_active_datapoints.tsv > ExCAPE_active_datapoints_Taxonid_GeneName.tsv
% cut -d$'\t' -f1,6 ExCAPE_active_datapoints.tsv > ExCAPE_active_datapoints_Compoundid_SMILES.tsv
% sed 1d ExCAPE_active_datapoints.tsv > ExCAPE_active_datapoints2.tsv (rename the file by removing 2)
% awk -F'\t' '{if($3=="" || $3 < 4.7) print $0}' ExCAPE_inactive_datapoints.tsv > ExCAPE_inactive_datapoints_filtered.tsv

% awk -F'\t' '{if(!$3) print $0}' ExCAPE_inactive_datapoints_filtered.tsv > ExCAPE_inactive_datapoints_filtered_noact.tsv
% awk -F'\t' '{if($3) print $0}' ExCAPE_inactive_datapoints_filtered.tsv > ExCAPE_inactive_datapoints_filtered_withact.tsv
% cut -d$'\t' -f4,5 ExCAPE_inactive_datapoints_filtered_noact.tsv > ExCAPE_inactive_datapoints_filtered_noact_Taxonid_GeneName.tsv
% cut -d$'\t' -f1,6 ExCAPE_inactive_datapoints_filtered_noact.tsv > ExCAPE_inactive_datapoints_filtered_noact_Compoundid_SMILES.tsv
% cut -d$'\t' -f4,5 ExCAPE_inactive_datapoints_filtered_withact.tsv > ExCAPE_inactive_datapoints_filtered_withact_Taxonid_GeneName.tsv
% cut -d$'\t' -f1,6 ExCAPE_inactive_datapoints_filtered_withact.tsv > ExCAPE_inactive_datapoints_filtered_withact_Compoundid_SMILES.tsv

% cut -d$'\t' -f1 DEEPScreen_ChEMBLv23_act_inact_filtered_data_original.txt > DEEPScreen_ChEMBLv23_act_inact_filtered_data_Targetid.txt
% cut -d$'\t' -f2 DEEPScreen_ChEMBLv23_act_inact_filtered_data_original.txt > DEEPScreen_ChEMBLv23_act_inact_filtered_data_Compoundid.txt
% 

% Mapping UniProt accessions to targets in ExCAPE dataset:

UniProt_mapping=cell(0,3);
[UniProt_UniProtAcc,UniProt_TaxonID,UniProt_GeneName]=textread('Mapping_UniProtAcc_TaxonID_GeneName.tab', '%s %s %s', 'delimiter', '\t', 'headerlines', 1);
to=1;
for i=1:length(UniProt_UniProtAcc)
    x=strsplit(UniProt_GeneName{i,1},' ')';
    UniProt_mapping(to:(to+length(x)-1),1)=repmat(UniProt_UniProtAcc(i,1),length(x),1);
    UniProt_mapping(to:(to+length(x)-1),2)=repmat(UniProt_TaxonID(i,1),length(x),1);
    UniProt_mapping(to:(to+length(x)-1),3)=x;
    to=to+length(x);
end
UniProt_mapping=upper(UniProt_mapping);
save UniProt_mapping.mat UniProt_mapping
UniProt_mapping_23=strcat(UniProt_mapping(:,2), {' '}, UniProt_mapping(:,3));
save UniProt_mapping_23.mat UniProt_mapping_23

% (ExCAPE all datapoints)
[ExCAPE_all_datapoints_Taxonid,ExCAPE_all_datapoints_GeneName]=textread('ExCAPE_all_datapoints_Taxonid_GeneName.tsv', '%s %s', 'delimiter', '\t');
ExCAPE_all_datapoints_Taxonid(1,:)=[];
ExCAPE_all_datapoints_GeneName(1,:)=[];
ExCAPE_all_datapoints_23=strcat(ExCAPE_all_datapoints_Taxonid(:,1), {' '}, ExCAPE_all_datapoints_GeneName(:,1));
clear ExCAPE_all_datapoints_Taxonid ExCAPE_all_datapoints_GeneName
load UniProt_mapping_23.mat
load UniProt_mapping.mat
[Lia,Locb]=ismember(ExCAPE_all_datapoints_23,UniProt_mapping_23);
ExCAPE_all_datapoints_UniProtAcc=cell(length(ExCAPE_all_datapoints_23),1);
ExCAPE_all_datapoints_UniProtAcc(Lia==1,1)=UniProt_mapping(Locb(Locb>0),1);
% save ExCAPE_all_datapoints_UniProtAcc.mat ExCAPE_all_datapoints_UniProtAcc
% (not possible to save as mat file due to extremely large size of the variable)
t=ExCAPE_all_datapoints_UniProtAcc';
fid = fopen('ExCAPE_all_datapoints_UniProtAcc.txt','w');
fprintf(fid,'%s\n',t{:});
fclose(fid);

% (ExCAPE active datapoints)
[ExCAPE_active_datapoints_Taxonid,ExCAPE_active_datapoints_GeneName]=textread('ExCAPE_active_datapoints_Taxonid_GeneName.tsv', '%s %s', 'delimiter', '\t');
ExCAPE_active_datapoints_23=strcat(ExCAPE_active_datapoints_Taxonid(:,1), {' '}, ExCAPE_active_datapoints_GeneName(:,1));
[Lia,Locb]=ismember(ExCAPE_active_datapoints_23,UniProt_mapping_23);
ExCAPE_active_datapoints_UniProtAcc=cell(length(ExCAPE_active_datapoints_23),1);
ExCAPE_active_datapoints_UniProtAcc(Lia==1,1)=UniProt_mapping(Locb(Locb>0),1);
save ExCAPE_active_datapoints_UniProtAcc.mat ExCAPE_active_datapoints_UniProtAcc
t=ExCAPE_active_datapoints_UniProtAcc';
fid = fopen('ExCAPE_active_datapoints_UniProtAcc.txt','w');
fprintf(fid,'%s\n',t{:});
fclose(fid);
[ExCAPE_active_datapoints_Compoundid,ExCAPE_active_datapoints_SMILES]=textread('ExCAPE_active_datapoints_Compoundid_SMILES.tsv', '%s %s', 'delimiter', '\t');
del_ind=find(cellfun(@isempty,ExCAPE_active_datapoints_UniProtAcc));
ExCAPE_active_datapoints_UniProtAcc(del_ind,:)=[];
ExCAPE_active_datapoints_Compoundid(del_ind,:)=[];
ExCAPE_active_datapoints_SMILES(del_ind,:)=[];
save ExCAPE_active_datapoints_Org_Var.mat ExCAPE_active_datapoints_UniProtAcc ExCAPE_active_datapoints_Compoundid ExCAPE_active_datapoints_SMILES
ExCAPE_Org_act=[ExCAPE_active_datapoints_UniProtAcc ExCAPE_active_datapoints_Compoundid];
save ExCAPE_Org_act.mat ExCAPE_Org_act

% (ExCAPE inactive datapoints with activity measures)
[ExCAPE_inactive_datapoints_filtered_withact_Taxonid,ExCAPE_inactive_datapoints_filtered_withact_GeneName]=textread('ExCAPE_inactive_datapoints_filtered_withact_Taxonid_GeneName.tsv', '%s %s', 'delimiter', '\t');
ExCAPE_inactive_datapoints_filtered_withact_23=strcat(ExCAPE_inactive_datapoints_filtered_withact_Taxonid(:,1), {' '}, ExCAPE_inactive_datapoints_filtered_withact_GeneName(:,1));
[Lia,Locb]=ismember(ExCAPE_inactive_datapoints_filtered_withact_23,UniProt_mapping_23);
ExCAPE_inactive_datapoints_filtered_withact_UniProtAcc=cell(length(ExCAPE_inactive_datapoints_filtered_withact_23),1);
ExCAPE_inactive_datapoints_filtered_withact_UniProtAcc(Lia==1,1)=UniProt_mapping(Locb(Locb>0),1);
save ExCAPE_inactive_datapoints_filtered_withact_UniProtAcc.mat ExCAPE_inactive_datapoints_filtered_withact_UniProtAcc
t=ExCAPE_inactive_datapoints_filtered_withact_UniProtAcc';
fid = fopen('ExCAPE_inactive_datapoints_filtered_withact_UniProtAcc.txt','w');
fprintf(fid,'%s\n',t{:});
fclose(fid);
[ExCAPE_inactive_datapoints_filtered_withact_Compoundid,ExCAPE_inactive_datapoints_filtered_withact_SMILES]=textread('ExCAPE_inactive_datapoints_filtered_withact_Compoundid_SMILES.tsv', '%s %s', 'delimiter', '\t');
del_ind=find(cellfun(@isempty,ExCAPE_inactive_datapoints_filtered_withact_UniProtAcc));
ExCAPE_inactive_datapoints_filtered_withact_UniProtAcc(del_ind,:)=[];
ExCAPE_inactive_datapoints_filtered_withact_Compoundid(del_ind,:)=[];
ExCAPE_inactive_datapoints_filtered_withact_SMILES(del_ind,:)=[];
save ExCAPE_inactive_datapoints_filtered_withact_Org_Var.mat ExCAPE_inactive_datapoints_filtered_withact_UniProtAcc ExCAPE_inactive_datapoints_filtered_withact_Compoundid ExCAPE_inactive_datapoints_filtered_withact_SMILES
ExCAPE_Org_inact_withact=[ExCAPE_inactive_datapoints_filtered_withact_UniProtAcc ExCAPE_inactive_datapoints_filtered_withact_Compoundid];
save ExCAPE_Org_inact_withact.mat ExCAPE_Org_inact_withact

ExCAPE_Org_act_inact=[ExCAPE_Org_act;ExCAPE_Org_inact_withact];
save ExCAPE_Org_act_inact.mat ExCAPE_Org_act_inact
length(unique(ExCAPE_Org_act_inact(:,1)))
length(unique(ExCAPE_Org_act_inact(:,2)))

% (ExCAPE inactive datapoints without any activity measures at all)
[ExCAPE_inactive_datapoints_filtered_noact_Taxonid,ExCAPE_inactive_datapoints_filtered_noact_GeneName]=textread('ExCAPE_inactive_datapoints_filtered_noact_Taxonid_GeneName.tsv', '%s %s', 'delimiter', '\t');
ExCAPE_inactive_datapoints_filtered_noact_23=strcat(ExCAPE_inactive_datapoints_filtered_noact_Taxonid(:,1), {' '}, ExCAPE_inactive_datapoints_filtered_noact_GeneName(:,1));
[Lia,Locb]=ismember(ExCAPE_inactive_datapoints_filtered_noact_23,UniProt_mapping_23);
ExCAPE_inactive_datapoints_filtered_noact_UniProtAcc=cell(length(ExCAPE_inactive_datapoints_filtered_noact_23),1);
ExCAPE_inactive_datapoints_filtered_noact_UniProtAcc(Lia==1,1)=UniProt_mapping(Locb(Locb>0),1);
% save ExCAPE_inactive_datapoints_filtered_noact_UniProtAcc.mat ExCAPE_inactive_datapoints_filtered_noact_UniProtAcc
% (not possible to save as mat file due to extremely large size of the variable)
t=ExCAPE_inactive_datapoints_filtered_noact_UniProtAcc';
fid = fopen('ExCAPE_inactive_datapoints_filtered_noact_UniProtAcc.txt','w');
fprintf(fid,'%s\n',t{:});
fclose(fid);


% (Generating the ExCAPE_PubChem dataset)
[~,idx]=sort(ExCAPE_Org_act(:,2));
[ExCAPE_Org_act_PubChem,idx]=sort(ExCAPE_Org_act(:,2));
ExCAPE_Org_act_PubChem=ExCAPE_Org_act_PubChem(1:201013,1);
ExCAPE_Org_act_PubChem_prot=ExCAPE_Org_act(idx,1);
ExCAPE_Org_act_PubChem_prot=ExCAPE_Org_act_PubChem_prot(1:201013,1);
ExCAPE_Org_act_PubChem=[ExCAPE_Org_act_PubChem_prot ExCAPE_Org_act_PubChem];
save ExCAPE_Org_act_PubChem.mat ExCAPE_Org_act_PubChem
[~,idx]=sort(ExCAPE_Org_inact_withact(:,2));
[ExCAPE_Org_inact_withact_PubChem,idx]=sort(ExCAPE_Org_inact_withact(:,2));
ExCAPE_Org_inact_withact_PubChem=ExCAPE_Org_inact_withact_PubChem(1:99946,1);
ExCAPE_Org_inact_withact_PubChem_prot=ExCAPE_Org_inact_withact(idx,1);
ExCAPE_Org_inact_withact_PubChem_prot=ExCAPE_Org_inact_withact_PubChem_prot(1:99946,1);
ExCAPE_Org_inact_withact_PubChem=[ExCAPE_Org_inact_withact_PubChem_prot ExCAPE_Org_inact_withact_PubChem];
save ExCAPE_Org_inact_withact_PubChem.mat ExCAPE_Org_inact_withact_PubChem
ExCAPE_Org_act_inact_PubChem=[ExCAPE_Org_act_PubChem;ExCAPE_Org_inact_withact_PubChem];
save ExCAPE_Org_act_inact_PubChem.mat ExCAPE_Org_act_inact_PubChem
length(unique(ExCAPE_Org_act_inact_PubChem(:,1)))
length(unique(ExCAPE_Org_act_inact_PubChem(:,2)))



% Mapping UniProt accessions to targets in DEEPScreen ChEMBLv23 training dataset:

[DEEPScreen_ChEMBLv23_act_inact_Target_ChEMBLid]=textread('DEEPScreen_ChEMBLv23_act_inact_filtered_data_Targetid.txt', '%s');
DEEPScreen_ChEMBLv23_act_inact_Target_ChEMBLid=strrep(DEEPScreen_ChEMBLv23_act_inact_Target_ChEMBLid,'_inact','');
DEEPScreen_ChEMBLv23_act_inact_Target_ChEMBLid=strrep(DEEPScreen_ChEMBLv23_act_inact_Target_ChEMBLid,'_act','');
[DEEPScreen_ChEMBLv23_act_inact_Compound_ChEMBLid]=textread('DEEPScreen_ChEMBLv23_act_inact_filtered_data_Compoundid.txt', '%s', 'delimiter', '\n', 'bufsize', 500000);

[Mapping_UniProtAcc,Mapping_ChEMBLid]=textread('Mapping_UniProtAcc_ChEMBLid.tab', '%s %s', 'delimiter', '\t', 'headerlines', 1);
[Lia,Locb]=ismember(DEEPScreen_ChEMBLv23_act_inact_Target_ChEMBLid,Mapping_ChEMBLid);
del_ind=find(Lia==0);
DEEPScreen_ChEMBLv23_act_inact_Target_ChEMBLid(del_ind,:)=[];
DEEPScreen_ChEMBLv23_act_inact_Compound_ChEMBLid(del_ind,:)=[];
[Lia,Locb]=ismember(DEEPScreen_ChEMBLv23_act_inact_Target_ChEMBLid,Mapping_ChEMBLid);
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
save DEEPScreen_ChEMBLv23_Org_act_raw.mat DEEPScreen_ChEMBLv23_Org_act
save DEEPScreen_ChEMBLv23_Org_inact_raw.mat DEEPScreen_ChEMBLv23_Org_inact

[ChEMBLv23_Compoundid,ChEMBLv23_SMILES]=textread('ChEMBLv23_compound_SMILES.txt', '%s %s', 'delimiter', '\t', 'headerlines', 1);
[Lia,Locb]=ismember(DEEPScreen_ChEMBLv23_Org_act(:,2),ChEMBLv23_Compoundid);
DEEPScreen_ChEMBLv23_Org_act_SMILES=cell(length(DEEPScreen_ChEMBLv23_Org_act),1);
DEEPScreen_ChEMBLv23_Org_act_SMILES(Lia==1,1)=ChEMBLv23_SMILES(Locb(Locb>0),1);
del_ind=find(Lia==0);
DEEPScreen_ChEMBLv23_Org_act(del_ind,:)=[];
DEEPScreen_ChEMBLv23_Org_act_SMILES(del_ind,:)=[];
save DEEPScreen_ChEMBLv23_Org_act_var.mat DEEPScreen_ChEMBLv23_Org_act DEEPScreen_ChEMBLv23_Org_act_SMILES

[Lia,Locb]=ismember(DEEPScreen_ChEMBLv23_Org_inact(:,2),ChEMBLv23_Compoundid);
DEEPScreen_ChEMBLv23_Org_inact_SMILES=cell(length(DEEPScreen_ChEMBLv23_Org_inact),1);
DEEPScreen_ChEMBLv23_Org_inact_SMILES(Lia==1,1)=ChEMBLv23_SMILES(Locb(Locb>0),1);
del_ind=find(Lia==0);
DEEPScreen_ChEMBLv23_Org_inact(del_ind,:)=[];
DEEPScreen_ChEMBLv23_Org_inact_SMILES(del_ind,:)=[];
save DEEPScreen_ChEMBLv23_Org_inact_var.mat DEEPScreen_ChEMBLv23_Org_inact DEEPScreen_ChEMBLv23_Org_inact_SMILES

DEEPScreen_ChEMBLv23_Org_act_inact=[DEEPScreen_ChEMBLv23_Org_act;DEEPScreen_ChEMBLv23_Org_inact];
save DEEPScreen_ChEMBLv23_Org_act_inact.mat DEEPScreen_ChEMBLv23_Org_act_inact
length(unique(DEEPScreen_ChEMBLv23_Org_act_inact(:,1)))
length(unique(DEEPScreen_ChEMBLv23_Org_act_inact(:,2)))



% Comparison and the merge between ExCAPE and DEEPScreen training datasets: 

% (actives dataset)
DEEPScreen_Org_act_51=strcat(DEEPScreen_ChEMBLv23_Org_act(:,1), {' '}, DEEPScreen_ChEMBLv23_Org_act(:,2));
ExCAPE_Org_act_51=strcat(ExCAPE_Org_act(:,1), {' '}, ExCAPE_Org_act(:,2));
[Lia,Locb]=ismember(ExCAPE_Org_act_51,DEEPScreen_Org_act_51);
ExCAPE_Org_act_onlyExCAPE=ExCAPE_Org_act(Lia==0,:);
ExCAPE_Org_act_onlyExCAPE_SMILES=ExCAPE_active_datapoints_SMILES(Lia==0,:);
save ExCAPE_Org_act_onlyExCAPE_var.mat ExCAPE_Org_act_onlyExCAPE ExCAPE_Org_act_onlyExCAPE_SMILES
Combined_Org_act=[DEEPScreen_ChEMBLv23_Org_act;ExCAPE_Org_act_onlyExCAPE];
Combined_Org_act_SMILES=[DEEPScreen_ChEMBLv23_Org_act_SMILES;ExCAPE_Org_act_onlyExCAPE_SMILES];
save Combined_Org_act_var.mat Combined_Org_act Combined_Org_act_SMILES
t=Combined_Org_act_SMILES';
fid = fopen('Combined_Org_act_SMILES.txt','w');
fprintf(fid,'%s\n',t{:});
fclose(fid);
Combined_Org_act_SMILES_id=[Combined_Org_act_SMILES Combined_Org_act(:,2)];
t=Combined_Org_act_SMILES_id';
fid = fopen('Combined_Org_act_SMILES_id.smi','w');
fprintf(fid,'%s\t%s\n',t{:});
fclose(fid);

% (inactives dataset)
DEEPScreen_Org_inact_51=strcat(DEEPScreen_ChEMBLv23_Org_inact(:,1), {' '}, DEEPScreen_ChEMBLv23_Org_inact(:,2));
ExCAPE_Org_inact_withact_51=strcat(ExCAPE_Org_inact_withact(:,1), {' '}, ExCAPE_Org_inact_withact(:,2));
[Lia,Locb]=ismember(ExCAPE_Org_inact_withact_51,DEEPScreen_Org_inact_51);
ExCAPE_Org_inact_withact_onlyExCAPE=ExCAPE_Org_inact_withact(Lia==0,:);
ExCAPE_Org_inact_withact_onlyExCAPE_SMILES=ExCAPE_inactive_datapoints_filtered_withact_SMILES(Lia==0,:);
save ExCAPE_Org_inact_withact_onlyExCAPE_var.mat ExCAPE_Org_inact_withact_onlyExCAPE ExCAPE_Org_inact_withact_onlyExCAPE_SMILES
Combined_Org_inact_withact=[DEEPScreen_ChEMBLv23_Org_inact;ExCAPE_Org_inact_withact_onlyExCAPE];
Combined_Org_inact_withact_SMILES=[DEEPScreen_ChEMBLv23_Org_inact_SMILES;ExCAPE_Org_inact_withact_onlyExCAPE_SMILES];
save Combined_Org_inact_withact_var.mat Combined_Org_inact_withact Combined_Org_inact_withact_SMILES
t=Combined_Org_inact_withact_SMILES';
fid = fopen('Combined_Org_inact_withact_SMILES.txt','w');
fprintf(fid,'%s\n',t{:});
fclose(fid);
Combined_Org_inact_withact_SMILES_id=[Combined_Org_inact_withact_SMILES Combined_Org_inact_withact(:,2)];
t=Combined_Org_inact_withact_SMILES_id';
fid = fopen('Combined_Org_inact_withact_SMILES_id.smi','w');
fprintf(fid,'%s\t%s\n',t{:});
fclose(fid);

Combined_Org_All_id=[Combined_Org_act(:,2);Combined_Org_inact_withact(:,2)];
[Combined_Org_All_id_unique,I,J] = unique(Combined_Org_All_id);
Combined_Org_All_SMILES=[Combined_Org_act_SMILES;Combined_Org_inact_withact_SMILES];
Combined_Org_All_SMILES_unique=Combined_Org_All_SMILES(I);
Combined_Org_All_SMILES_id_unique=[Combined_Org_All_SMILES_unique Combined_Org_All_id_unique];
t=Combined_Org_All_SMILES_id_unique';
fid = fopen('Combined_Org_All_SMILES_id_unique.smi','w');
fprintf(fid,'%s\t%s\n',t{:});
fclose(fid);

Combined_Org_act_51=strcat(Combined_Org_act(:,1), {' '}, Combined_Org_act(:,2));
Combined_Org_inact_withact_51=strcat(Combined_Org_inact_withact(:,1), {' '}, Combined_Org_inact_withact(:,2));
[Lia,Locb]=ismember(Combined_Org_act_51,Combined_Org_inact_withact_51);
sum(Lia)
% (there are very few corresponding data points between actives and inactives: 1574)

Combined_Org_act_inact=[Combined_Org_act;Combined_Org_inact_withact];
save Combined_Org_act_inact.mat Combined_Org_act_inact
length(unique(Combined_Org_act_inact(:,1)))
length(unique(Combined_Org_act_inact(:,2)))
Target_UniProtacc_Combined_Org_act_inact=unique(Combined_Org_act_inact(:,1));
save Target_UniProtacc_Combined_Org_act_inact.mat Target_UniProtacc_Combined_Org_act_inact



% Histogram of data point distributions for all targets:

length(unique(Combined_Org_act(:,1)))
length(unique(Combined_Org_act(:,2)))
[Lia,Locb]=ismember(Combined_Org_act(:,1),unique(Combined_Org_act(:,1)));
Hist_combined_act=hist(Locb,0.5:1:(max(Locb)+0.5));
figure;bar(log10(sort(Hist_combined_act,'descend')))
[Lia,Locb]=ismember(Combined_Org_act(:,2),unique(Combined_Org_act(:,2)));
Hist_combined_act_comp=sort(hist(Locb,0.5:1:(max(Locb)+0.5)),'descend')';
length(find(Hist_combined_act_comp>=5))

length(unique(Combined_Org_inact_withact(:,1)))
length(unique(Combined_Org_inact_withact(:,2)))
[Lia,Locb]=ismember(Combined_Org_inact_withact(:,1),unique(Combined_Org_inact_withact(:,1)));
Hist_combined_inact_withact=hist(Locb,0.5:1:(max(Locb)+0.5));
figure;bar(log10(sort(Hist_combined_inact_withact,'descend')))
[Lia,Locb]=ismember(Combined_Org_inact_withact(:,2),unique(Combined_Org_inact_withact(:,2)));
Hist_combined_inact_withact_comp=sort(hist(Locb,0.5:1:(max(Locb)+0.5)),'descend')';
length(find(Hist_combined_inact_withact_comp>=5))



% Generation of the Morgan Fingerprints for all compounds and pairwise tanimoto similarity search using Chemfp:

% conda create -n py2condaenv python=2
% conda activate py2-rdkit-env
% python /Users/tuncadogan/Downloads/chemfp-1.5/setup.py install --without-openmp
% cd /Users/tuncadogan/Documents/METU-II_BIN_Research_Projects/DomainMappingLIMK_Study_Data_Analysis
% rdkit2fps --morgan Combined_Org_All_SMILES_id_unique.smi -o Combined_Org_All_ECFP4_id_unique.fps
% split -l 206705 Combined_Org_All_ECFP4_id_unique.fps
% simsearch  --times --threshold 0.5 -q xaa xaa -o xaa_xaa_Simsearch.txt
% simsearch  --times --threshold 0.5 -q xaa xab -o xaa_xab_Simsearch.txt
% simsearch  --times --threshold 0.5 -q xaa xac -o xaa_xac_Simsearch.txt
% simsearch  --times --threshold 0.5 -q xaa xad -o xaa_xad_Simsearch.txt
% simsearch  --times --threshold 0.5 -q xaa xae -o xaa_xae_Simsearch.txt
% simsearch  --times --threshold 0.5 -q xab xab -o xab_xab_Simsearch.txt
% simsearch  --times --threshold 0.5 -q xab xac -o xab_xac_Simsearch.txt
% simsearch  --times --threshold 0.5 -q xab xad -o xab_xad_Simsearch.txt
% simsearch  --times --threshold 0.5 -q xab xae -o xab_xae_Simsearch.txt
% simsearch  --times --threshold 0.5 -q xac xac -o xac_xac_Simsearch.txt
% simsearch  --times --threshold 0.5 -q xac xad -o xac_xad_Simsearch.txt
% simsearch  --times --threshold 0.5 -q xac xae -o xac_xae_Simsearch.txt
% simsearch  --times --threshold 0.5 -q xad xad -o xad_xad_Simsearch.txt
% simsearch  --times --threshold 0.5 -q xad xae -o xad_xae_Simsearch.txt
% simsearch  --times --threshold 0.5 -q xae xae -o xae_xae_Simsearch.txt
% simsearch  --times --threshold 0.5 --NxN Combined_Org_All_ECFP4_id_unique.fps -o Combined_Org_All_Simsearch.txt



% Organizing the compound similarity files (repeat below code for all 15 parts from xaa to xae):

[Part_Simsearch]=textread('xaa_xaa_Simsearch.txt', '%s', 'delimiter', '\n', 'bufsize', 500000);
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
fid=fopen('xaa_xaa_SimS_Org.txt','w');
fprintf(fid,'%s\t%s\t%s\n',t{:});
fclose(fid);

% cat xaa_xaa_SimS_Org.txt xaa_xab_SimS_Org.txt xaa_xac_SimS_Org.txt xaa_xad_SimS_Org.txt xaa_xae_SimS_Org.txt xab_xab_SimS_Org.txt xab_xac_SimS_Org.txt xab_xad_SimS_Org.txt xab_xae_SimS_Org.txt xac_xac_SimS_Org.txt xac_xad_SimS_Org.txt xac_xae_SimS_Org.txt xad_xad_SimS_Org.txt xad_xae_SimS_Org.txt xae_xae_SimS_Org.txt > Combined_SimS_Org.txt



% Generating the combined protein groups for each ligand cluster using the selected ligand similarity threshold:

[Combined_Org_Comp1,Combined_Org_Comp2,Combined_Org_Sim]=textread('Combined_SimS_Org.txt', '%s %s %d', 'delimiter', '\t');
load Combined_Org_act_var.mat
load Combined_Org_inact_withact_var.mat
Combined_Org_id_unique=unique([Combined_Org_act(:,2);Combined_Org_inact_withact(:,2)]);
save Combined_Org_id_unique.mat Combined_Org_id_unique

thr=0.7;

ind_thr=find(Combined_Org_Sim>=thr);
Combined_Org_Sim=Combined_Org_Sim(ind_thr);
Combined_Org_Comp1=Combined_Org_Comp1(ind_thr);
Combined_Org_Comp2=Combined_Org_Comp2(ind_thr);
save Combined_Org_Comp1_07.mat Combined_Org_Comp1
save Combined_Org_Comp2_07.mat Combined_Org_Comp2
save Combined_Org_Sim_07.mat Combined_Org_Sim

[Lia1,Locb1]=ismember(Combined_Org_Comp1,Combined_Org_id_unique);
[Lia2,Locb2]=ismember(Combined_Org_Comp2,Combined_Org_id_unique);
save Locations_07.mat Locb1 Locb2
[Lia_act,Locb_act]=ismember(Combined_Org_act(:,2),Combined_Org_id_unique);
[Lia_inact,Locb_inact]=ismember(Combined_Org_inact_withact(:,2),Combined_Org_id_unique);
save Locations_act_inact.mat Locb_act Locb_inact
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
save Combined_Org_unique_act_Prot_07.mat Combined_Org_unique_act_Prot
save Combined_Org_unique_inact_withact_Prot_07.mat Combined_Org_unique_inact_withact_Prot



% Analyzing and thresholding the protein arrays:

Combined_Org_unique_Samplenum_act_inact_sup=zeros(length(Combined_Org_id_unique),3);
for i=1:length(Combined_Org_id_unique)
    disp(['Compound number: ' num2str(i), ' / ' num2str(length(Combined_Org_id_unique))])
    Combined_Org_unique_Samplenum_act_inact_sup(i,1)=length(Combined_Org_unique_act_Prot{i,1});
    Combined_Org_unique_Samplenum_act_inact_sup(i,2)=length(Combined_Org_unique_inact_withact_Prot{i,1});
    Combined_Org_unique_Samplenum_act_inact_sup(i,3)=(Combined_Org_unique_Samplenum_act_inact_sup(i,1)+Combined_Org_unique_Samplenum_act_inact_sup(i,2))/2;
end
save Combined_Org_unique_Samplenum_act_inact_sup.mat Combined_Org_unique_Samplenum_act_inact_sup

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
save Combined_Org_unique_act_inact_Prot_5_var.mat Combined_Org_unique_Samplenum_act_inact_sup_5 Combined_Org_unique_act_Prot_5 Combined_Org_unique_inact_withact_Prot_5 Combined_Org_id_unique_Prot_5



% InterPro domain array generation for active and inactive sets:

% (saving a unique list of target protein accessions from the combined activity set)
Combined_Org_act_inact_Prot_unique=unique([Combined_Org_act(:,1);Combined_Org_inact_withact(:,1)]);
save Combined_Org_act_inact_Prot_unique.mat Combined_Org_act_inact_Prot_unique
% (UniProt id retrieval is used to obtain InterPro hits for these 3644 target proteins)


% (Organizing InterPro hits of our target proteins)
[InterPro_entriesmatched_onlytarget_UniProtAcc,InterPro_entriesmatched_onlytarget_IPRid]=textread('InterPro_UniProtAcc_IPRid_entriesmatched_onlytargets.tab', '%s %s', 'delimiter', '\t', 'headerlines', 1);
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
save InterPro_entriesmatched_onlytarget_Org.mat InterPro_entriesmatched_onlytarget_Org

% (Organizing InterPro hits of all ChEMBL targets)
[InterPro_entriesmatched_onlyChEMBL_UniProtAcc,InterPro_entriesmatched_onlyChEMBL_IPRid]=textread('InterPro_UniProtAcc_IPRid_entriesmatched_onlyChEMBL.tab', '%s %s', 'delimiter', '\t', 'headerlines', 1);
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
save InterPro_entriesmatched_onlyChEMBL_Org.mat InterPro_entriesmatched_onlyChEMBL_Org

% (Organizing InterPro hits of all human proteins)
[InterPro_entriesmatched_allhuman_UniProtAcc,InterPro_entriesmatched_allhuman_IPRid]=textread('InterPro_UniProtAcc_IPRid_entriesmatched_allhuman.tab', '%s %s', 'delimiter', '\t', 'headerlines', 1);
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
save InterPro_entriesmatched_allhuman_Org.mat InterPro_entriesmatched_allhuman_Org


% (Selecting only domain type hits from the organized protein InterPro hit files)
[InterPro_v72_Domain_id,~,~]=textread('InterPro_v72_Domain_entrylist.txt', '%s %s %s', 'delimiter', '\t');
InterPro_v72_Domain_id=unique(InterPro_v72_Domain_id);

[Lia,Locb]=ismember(InterPro_entriesmatched_onlytarget_Org(:,2),InterPro_v72_Domain_id);
InterPro_entriesmatched_onlytarget_Org_onlydom=InterPro_entriesmatched_onlytarget_Org(Lia==1,:);
save InterPro_entriesmatched_onlytarget_Org_onlydom.mat InterPro_entriesmatched_onlytarget_Org_onlydom

[Lia,Locb]=ismember(InterPro_entriesmatched_onlyChEMBL_Org(:,2),InterPro_v72_Domain_id);
InterPro_entriesmatched_onlyChEMBL_Org_onlydom=InterPro_entriesmatched_onlyChEMBL_Org(Lia==1,:);
save InterPro_entriesmatched_onlyChEMBL_Org_onlydom.mat InterPro_entriesmatched_onlyChEMBL_Org_onlydom

[Lia,Locb]=ismember(InterPro_entriesmatched_allhuman_Org(:,2),InterPro_v72_Domain_id);
InterPro_entriesmatched_allhuman_Org_onlydom=InterPro_entriesmatched_allhuman_Org(Lia==1,:);
save InterPro_entriesmatched_allhuman_Org_onlydom.mat InterPro_entriesmatched_allhuman_Org_onlydom


% (organizing InterPro term relation file to determine domain relations)
[InterPro_v72_ParentChildTree,~]=textread('InterPro_v72_ParentChildTreeFile.txt', '%s %s', 'delimiter', '::');
InterPro_v72_ParentChildTree=InterPro_v72_ParentChildTree(1:2:end,1);

% (removal of non-domain IRP entries from the hierarchy file)
for i=1:length(InterPro_v72_ParentChildTree)
    InterPro_v72_terms_from_ParentChildTree(i,1)=erase(InterPro_v72_ParentChildTree(i,1),"-");
end
[Lia,Locb]=ismember(InterPro_v72_terms_from_ParentChildTree,InterPro_v72_Domain_id);
InterPro_v72_ParentChildTree_Dom=InterPro_v72_ParentChildTree(Lia==1,:);
% (addition of domain entries to the table that have no hierarchical relationship with other entries)
InterPro_v72_terms_without_hierarc=setdiff(InterPro_v72_Domain_id,InterPro_v72_terms_from_ParentChildTree);
InterPro_v72_ParentChildTree_Dom=[InterPro_v72_ParentChildTree_Dom;InterPro_v72_terms_without_hierarc];
save InterPro_v72_ParentChildTree_var.mat InterPro_v72_ParentChildTree InterPro_v72_ParentChildTree_Dom

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
save InterPro_v72_Dom_groups.mat InterPro_v72_Dom_groups



% Associating compounds with single domains:

% Calculating the compound cluster (similarity based) single domain assocation scores for all compounds:

load Combined_Org_unique_act_inact_Prot_5_var.mat
load InterPro_v72_Dom_groups.mat
load InterPro_entriesmatched_onlytarget_Org_onlydom.mat
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
save Domain_mapping_Prot_5.mat Domain_mapping_Prot_5

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
del_ind=find(cellfun(@isempty,Domain_mapping_Prot_5_all_Org(:,1)));
Domain_mapping_Prot_5_all_Org(del_ind,:)=[];
length(unique(Domain_mapping_Prot_5_all_Org(:,1)))
length(unique(Domain_mapping_Prot_5_all_Org(:,2)))
t=Domain_mapping_Prot_5_all_Org';
fid=fopen('Domain_mapping_Prot_5_all_Org.txt','w');
fprintf(fid,'%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n',t{:});
fclose(fid);
Domain_mapping_Prot_5_all_Org_part1=Domain_mapping_Prot_5_all_Org(1:1150000,:);
Domain_mapping_Prot_5_all_Org_part2=Domain_mapping_Prot_5_all_Org(1150001:2300000,:);
Domain_mapping_Prot_5_all_Org_part3=Domain_mapping_Prot_5_all_Org(2300001:end,:);
save Domain_mapping_Prot_5_all_Org_part1.mat Domain_mapping_Prot_5_all_Org_part1
save Domain_mapping_Prot_5_all_Org_part2.mat Domain_mapping_Prot_5_all_Org_part2
save Domain_mapping_Prot_5_all_Org_part3.mat Domain_mapping_Prot_5_all_Org_part3

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
save Domain_mapping_Prot_5_all_MCC07.mat Domain_mapping_Prot_5_all_MCC07
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
save Domain_mapping_Prot_5_all_Acc08.mat Domain_mapping_Prot_5_all_Acc08
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
save Domain_mapping_Prot_5_all_F108.mat Domain_mapping_Prot_5_all_F108
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
save Domain_mapping_Prot_5_all_MulFil.mat Domain_mapping_Prot_5_all_MulFil
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

load Combined_Org_act_var.mat
load Combined_Org_inact_withact_var.mat
load Combined_Org_id_unique.mat
load Locations_act_inact.mat

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
save Combined_Org_unique_act_Prot_noSim.mat Combined_Org_unique_act_Prot_noSim
save Combined_Org_unique_inact_withact_Prot_noSim.mat Combined_Org_unique_inact_withact_Prot_noSim

% (analyzing and thresholding the protein arrays)
Combined_Org_unique_Samplenum_act_inact_noSim_sup=zeros(length(Combined_Org_id_unique),3);
for i=1:length(Combined_Org_id_unique)
    disp(['Compound number: ' num2str(i), ' / ' num2str(length(Combined_Org_id_unique))])
    Combined_Org_unique_Samplenum_act_inact_noSim_sup(i,1)=length(Combined_Org_unique_act_Prot_noSim{i,1});
    Combined_Org_unique_Samplenum_act_inact_noSim_sup(i,2)=length(Combined_Org_unique_inact_withact_Prot_noSim{i,1});
    Combined_Org_unique_Samplenum_act_inact_noSim_sup(i,3)=(Combined_Org_unique_Samplenum_act_inact_noSim_sup(i,1)+Combined_Org_unique_Samplenum_act_inact_noSim_sup(i,2))/2;
end
save Combined_Org_unique_Samplenum_act_inact_noSim_sup.mat Combined_Org_unique_Samplenum_act_inact_noSim_sup

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
save Combined_Org_unique_act_inact_Prot_noSim_3_var.mat Combined_Org_unique_Samplenum_act_inact_noSim_sup_3 Combined_Org_unique_act_Prot_noSim_3 Combined_Org_unique_inact_withact_Prot_noSim_3 Combined_Org_id_unique_Prot_noSim_3

% (calculating the single domain assocation scores for all compounds)
load Combined_Org_unique_act_inact_Prot_noSim_3_var.mat
load InterPro_v72_Dom_groups.mat
load InterPro_entriesmatched_onlytarget_Org_onlydom.mat
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
save Domain_mapping_Prot_noSim_3.mat Domain_mapping_Prot_noSim_3

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
del_ind=find(cellfun(@isempty,Domain_mapping_Prot_noSim_3_all_Org(:,1)));
Domain_mapping_Prot_noSim_3_all_Org(del_ind,:)=[];
length(unique(Domain_mapping_Prot_noSim_3_all_Org(:,1)))
length(unique(Domain_mapping_Prot_noSim_3_all_Org(:,2)))
t=Domain_mapping_Prot_noSim_3_all_Org';
fid=fopen('Domain_mapping_Prot_noSim_3_all_Org.txt','w');
fprintf(fid,'%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n',t{:});
fclose(fid);
save Domain_mapping_Prot_noSim_3_all_Org.mat Domain_mapping_Prot_noSim_3_all_Org

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
save Domain_mapping_Prot_noSim_3_all_MulFil.mat Domain_mapping_Prot_noSim_3_all_MulFil
length(Domain_mapping_Prot_noSim_3_all_MulFil)
length(unique(Domain_mapping_Prot_noSim_3_all_MulFil(:,1)))
length(unique(Domain_mapping_Prot_noSim_3_all_MulFil(:,2)))



% Merging compound similarity based domain associations with noSim associations:

load Domain_mapping_Prot_5_all_MulFil.mat
load Domain_mapping_Prot_noSim_3_all_MulFil.mat
% (extracting mutual associations and selecting maximum performance according to MCC for merging)
Domain_mapping_Prot_5_all_MulFil_23=strcat(Domain_mapping_Prot_5_all_MulFil(:,1), {' '}, Domain_mapping_Prot_5_all_MulFil(:,2));
Domain_mapping_Prot_noSim_3_all_MulFil_23=strcat(Domain_mapping_Prot_noSim_3_all_MulFil(:,1), {' '}, Domain_mapping_Prot_noSim_3_all_MulFil(:,2));
[IntSec,I1,I2]=intersect(Domain_mapping_Prot_5_all_MulFil_23,Domain_mapping_Prot_noSim_3_all_MulFil_23);

C1=(cell2mat(Domain_mapping_Prot_5_all_MulFil(I1,end))>cell2mat(Domain_mapping_Prot_noSim_3_all_MulFil(I2,end)));
mer1=Domain_mapping_Prot_5_all_MulFil(I1(C1==1),:);
mer2=Domain_mapping_Prot_noSim_3_all_MulFil(I2(C1==0),:);
I1inv=setdiff((1:length(Domain_mapping_Prot_5_all_MulFil)),I1)';
I2inv=setdiff((1:length(Domain_mapping_Prot_noSim_3_all_MulFil)),I2)';

% (saving the merged and filtered domain accosiation file)
Domain_mapping_Prot_Sim5_noSim3_MulFil_merged=[Domain_mapping_Prot_5_all_MulFil(I1inv,:);Domain_mapping_Prot_noSim_3_all_MulFil(I2inv,:);mer1;mer2];
save Domain_mapping_Prot_Sim5_noSim3_MulFil_merged.mat Domain_mapping_Prot_Sim5_noSim3_MulFil_merged
t=Domain_mapping_Prot_Sim5_noSim3_MulFil_merged';
fid=fopen('Domain_mapping_Prot_Sim5_noSim3_MulFil_merged_finalized.txt','w');
fprintf(fid,'%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n',t{:});
fclose(fid);
length(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged)
length(unique(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged(:,1)))
length(unique(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged(:,2)))
% (final domain association stats: 27,032 mappings between 8,165 compounds and 250 InterPro domains)



% Associating compounds with domain pairs (where it performs better
% compared to mappings to single domains in the corresponding pairs):

% (do not use DAs for this, selection of only 1 DA out of many results in omission of some significant domain pairs, instead use bag of domains obtained directly from UniProt, variable: 'InterPro_entriesmatched_onlytarget_Org_onlydom')

% (generating an array for InterPro entry pairs in the same hierarchy, so that these pairs can be subtracted from all domain pairs since these pairs are not 2 domains hits but 2 related entries hit to the same position)
load InterPro_v72_Dom_groups.mat
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
save InterPro_domain_relation_pairs.mat InterPro_domain_relation_pairs

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
save InterPro_entriesmatched_onlytarget_Org_bagofdom_onlydompairs.mat InterPro_entriesmatched_onlytarget_Org_bagofdom_onlydompairs

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
save InterPro_entriesmatched_allhuman_Org_bagofdom_onlydompairs.mat InterPro_entriesmatched_allhuman_Org_bagofdom_onlydompairs


% % Generation and organisation of the domain architectures for all ChEMBL and PubChem targets (change current directory to: /Documents/EBI Post-doctoral Research Study/):
% 
% % (generation of the all-targets protein file)
% load Target_UniProtacc_Combined_Org_act_inact.mat
% [InterPro_entriesmatched_onlyChEMBL_UniProtAcc,InterPro_entriesmatched_onlyChEMBL_IPRid]=textread('InterPro_UniProtAcc_IPRid_entriesmatched_onlyChEMBL.tab', '%s %s', 'delimiter', '\t', 'headerlines', 1);
% All_Targets_UniProtAcc=unique([InterPro_entriesmatched_onlyChEMBL_UniProtAcc;Target_UniProtacc_Combined_Org_act_inact]);
% save All_Targets_UniProtAcc.mat All_Targets_UniProtAcc
% 
% % (organisation of the InterPro_v72 annotation file)
% eval(['mkdir interpro_hit_file_v72'])
% unix(['split -l 5000 -a 4 protein2ipr.dat interpro_hit_file_v72/prot2ipr_'])
% no_of_files=length(dir('interpro_hit_file_v72'))-2;
% Interpro_dom_ACC=textread('InterPro_v72_Domain_entrylist_onlyid.txt', '%s', 'delimiter', '\t');
% save Interpro_dom_ACC_v72.mat Interpro_dom_ACC -v7
% entry_hierarchy_gen('InterPro_v72_ParentChildTreeFile.txt');
% Ipr_dom_hierarchy_groups=hierarchy_grouping('InterPro_v72_ParentChildTreeFile.txt','Interpro_dom_ACC_v72.mat');
% save Ipr_dom_v72_hierarchy_groups.mat Ipr_dom_hierarchy_groups
% 
% % (organisation of the input protein file using all ChEMBL and PubChem targets):
% [ChEMBL_Targets_UniProt_Acc,ChEMBL_Targets_UniProt_Size]=textread('ChEMBL_Targets_UniProt_Acc_Size_2019_01.tsv', '%s %d', 'delimiter', '\t', 'headerlines', 1);
% [ChEMBL_Targets_UniProt_Acc,idx]=sort(ChEMBL_Targets_UniProt_Acc);
% ChEMBL_Targets_UniProt_Size=ChEMBL_Targets_UniProt_Size(idx,1);
% save ChEMBL_Targets_UniProt_Acc.mat ChEMBL_Targets_UniProt_Acc
% save ChEMBL_Targets_UniProt_Size.mat ChEMBL_Targets_UniProt_Size
% fid=fopen('ChEMBL_Targets_UniProt_Acc.txt', 'w');
% fprintf(fid, '%s\n', ChEMBL_Targets_UniProt_Acc{:});
% fclose(fid);
% fid=fopen('ChEMBL_Targets_UniProt_Size.txt', 'w');
% fprintf(fid, '%d\n', ChEMBL_Targets_UniProt_Size);
% fclose(fid);
% 
% % (organisation of the input protein file using all human proteins):
% [HumanProt_UniProt_Acc,HumanProt_UniProt_Size]=textread('HumanProt_UniProt_Acc_Size_2019_01.tsv', '%s %d', 'delimiter', '\t', 'headerlines', 1);
% [HumanProt_UniProt_Acc,idx]=sort(HumanProt_UniProt_Acc);
% HumanProt_UniProt_Size=HumanProt_UniProt_Size(idx,1);
% save HumanProt_UniProt_Acc.mat HumanProt_UniProt_Acc
% save HumanProt_UniProt_Size.mat HumanProt_UniProt_Size
% fid=fopen('HumanProt_UniProt_Acc.txt', 'w');
% fprintf(fid, '%s\n', HumanProt_UniProt_Acc{:});
% fclose(fid);
% fid=fopen('HumanProt_UniProt_Size.txt', 'w');
% fprintf(fid, '%d\n', HumanProt_UniProt_Size);
% fclose(fid);
% 
% % (generation of the DAs)
% no_of_files=length(dir('interpro_hit_file_v72'))-2;
% Ipr_DA_generation_run('ChEMBL_Targets_UniProt_Acc.txt','ChEMBL_Targets_UniProt_Size.txt','Interpro_dom_ACC_v72.mat','interpro_hit_file_v72',no_of_files,'entry_hierarchy_v72.csv','Workfiles/ChEMBL_Targets_UniProt_2019_01_DA_generation',0,0,1,[],[],1)
% Ipr_DA_generation_run('HumanProt_UniProt_Acc.txt','HumanProt_UniProt_Size.txt','Interpro_dom_ACC_v72.mat','interpro_hit_file_v72',no_of_files,'entry_hierarchy_v72.csv','Workfiles/HumanProt_UniProt_2019_01_DA_generation',0,0,1,[],[],1)
% 
% % (organising the DA output files)
% [ChEMBL_Targets_UniProt_DA_assoc_proteins]=textread('ChEMBL_Targets_UniProt_DA_assoc_proteins.txt', '%s', 'delimiter', '\n', 'bufsize', 500000);
% [ChEMBL_Targets_UniProt_DA_output]=textread('ChEMBL_Targets_UniProt_DA_output.txt', '%s', 'delimiter', '\n', 'bufsize', 500000);
% ChEMBL_Targets_UniProt_DA_Org=cell(length(ChEMBL_Targets_UniProt_Acc),1);
% ChEMBL_Targets_UniProt_DomNum_Matrix=zeros(length(ChEMBL_Targets_UniProt_Acc),10);
% ChEMBL_Targets_UniProt_DomNum_ind=zeros(length(ChEMBL_Targets_UniProt_Acc),1);
% for i=1:length(ChEMBL_Targets_UniProt_DA_output)
%     disp(['Organising Domain number: ' num2str(i), ' / ' num2str(length(ChEMBL_Targets_UniProt_DA_output))])
%     DA_temp=strsplit(ChEMBL_Targets_UniProt_DA_output{i,1},'\t')';
%     DA_temp(1,:)=[];
%     DomNum_temp=length(DA_temp);
%     prot_temp=strsplit(ChEMBL_Targets_UniProt_DA_assoc_proteins{i,1},'\t')';
%     prot_temp(1,:)=[];
%     [Lia,Locb]=ismember(prot_temp,ChEMBL_Targets_UniProt_Acc);
%     for j=1:length(Locb)
%         ChEMBL_Targets_UniProt_DomNum_ind(Locb(j,1),1)=ChEMBL_Targets_UniProt_DomNum_ind(Locb(j,1),1)+1;
%         ChEMBL_Targets_UniProt_DomNum_Matrix(Locb(j,1),ChEMBL_Targets_UniProt_DomNum_ind(Locb(j,1),1))=DomNum_temp;
%         ChEMBL_Targets_UniProt_DA_Org{Locb(j,1),1}=[ChEMBL_Targets_UniProt_DA_Org{Locb(j,1),1};{DA_temp}];
%     end
% end
% 
% % (selecting one DA for each protein, which should be the DA with minimum number of domains among the DAs with multiple domains)
% ChEMBL_Targets_UniProt_DA_Org_sel=cell(length(ChEMBL_Targets_UniProt_DA_Org),1);
% for i=1:length(ChEMBL_Targets_UniProt_DA_Org)
%     if size(ChEMBL_Targets_UniProt_DA_Org{i,1},1)>1
%         ind_temp=find(ChEMBL_Targets_UniProt_DomNum_Matrix(i,:)>1);
%         if isempty(ind_temp)==0
%             Mat_temp=ChEMBL_Targets_UniProt_DomNum_Matrix(i,ind_temp);
%             ind=find(ChEMBL_Targets_UniProt_DomNum_Matrix(i,:)==min(Mat_temp));
%             ChEMBL_Targets_UniProt_DA_Org_sel{i,1}=ChEMBL_Targets_UniProt_DA_Org{i,1}{ind(1,1),1};
%         else
%             ChEMBL_Targets_UniProt_DA_Org_sel{i,1}=ChEMBL_Targets_UniProt_DA_Org{i,1}{1,1};
%         end
%     elseif size(ChEMBL_Targets_UniProt_DA_Org{i,1},1)==1
%         ChEMBL_Targets_UniProt_DA_Org_sel{i,1}=ChEMBL_Targets_UniProt_DA_Org{i,1}{1,1};
%     end
% end
% 
% % (preparing the domain pairs for each protein as only pairs will be considered for the mapping to compounds)
% ChEMBL_Targets_UniProt_DA_Org_sel_DomPairs=cell(length(ChEMBL_Targets_UniProt_DA_Org),1);
% for i=1:length(ChEMBL_Targets_UniProt_DA_Org_sel)
%     disp(['Organising the DA of the protein number: ' num2str(i), ' / ' num2str(length(ChEMBL_Targets_UniProt_DA_Org_sel))])
%     if size(ChEMBL_Targets_UniProt_DA_Org_sel{i,1},1)>1 && size(ChEMBL_Targets_UniProt_DA_Org_sel{i,1},1)<6
%         clear sel_temp
%         for j=1:size(ChEMBL_Targets_UniProt_DA_Org_sel{i,1},1)
%             sel_temp(j,1:length(strsplit(ChEMBL_Targets_UniProt_DA_Org_sel{i,1}{j,1},'/')))=strsplit(ChEMBL_Targets_UniProt_DA_Org_sel{i,1}{j,1},'/');
%         end
%         if size(sel_temp,2)>1
%             for j=2:size(sel_temp,2)
%                 indd=find(cellfun(@isempty,sel_temp(:,j)));
%                 if isempty(indd)==0
%                     sel_temp(indd,j)=sel_temp(indd,(j-1));
%                 end
%             end
%         end
%         for j=1:size(sel_temp,2)
%             ChEMBL_Targets_UniProt_DA_Org_sel_DomPairs{i,1}=[ChEMBL_Targets_UniProt_DA_Org_sel_DomPairs{i,1};nchoosek(sel_temp(:,j),2)];
%         end
%         ChEMBL_Targets_UniProt_DA_Org_sel_DomPairs{i,1}=uniqueRowsCA(ChEMBL_Targets_UniProt_DA_Org_sel_DomPairs{i,1});
%     end
% end
% save ChEMBL_Targets_UniProt_DA_var.mat ChEMBL_Targets_UniProt_DA_Org ChEMBL_Targets_UniProt_DomNum_Matrix ChEMBL_Targets_UniProt_DomNum_ind ChEMBL_Targets_UniProt_DA_Org_sel ChEMBL_Targets_UniProt_DA_output ChEMBL_Targets_UniProt_DA_assoc_proteins ChEMBL_Targets_UniProt_Acc ChEMBL_Targets_UniProt_DA_Org_sel_DomPairs
% 
% % (organising the protein domain pair annotation file by adding the domain pairs in inverse order and by concatanating the domain ids)
% InterPro_entriesmatched_onlytarget_Org_onlydompairs=cell(0,2);
% to=1;
% for i=1:length(ChEMBL_Targets_UniProt_Acc)
%     disp(['Line number: ' num2str(i), ' / ' num2str(length(ChEMBL_Targets_UniProt_Acc))])
%     if isempty(ChEMBL_Targets_UniProt_DA_Org_sel_DomPairs{i,1})==0
%         x=ChEMBL_Targets_UniProt_DA_Org_sel_DomPairs{i,1};
%         x=uniqueRowsCA([x;x(:,[2 1])]);
%         x=strcat(x(:,1), {'-'}, x(:,2));
%         InterPro_entriesmatched_onlytarget_Org_onlydompairs(to:(to+(size(x,1))-1),1)=repmat(ChEMBL_Targets_UniProt_Acc(i,1),size(x,1),1);
%         InterPro_entriesmatched_onlytarget_Org_onlydompairs(to:(to+(size(x,1))-1),2)=x;
%         to=to+(size(x,1));
%     end
% end
% save InterPro_entriesmatched_onlytarget_Org_onlydompairs.mat InterPro_entriesmatched_onlytarget_Org_onlydompairs



% Calculating the compound cluster (similarity based) domain pair association scores for all compounds:

load InterPro_entriesmatched_onlytarget_Org_bagofdom_onlydompairs.mat
load Combined_Org_unique_act_inact_Prot_5_var.mat
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
save DomainPair_mapping_Prot_5.mat DomainPair_mapping_Prot_5

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
del_ind=find(cellfun(@isempty,DomainPair_mapping_Prot_5_all_Org(:,1)));
DomainPair_mapping_Prot_5_all_Org(del_ind,:)=[];
length(unique(DomainPair_mapping_Prot_5_all_Org(:,1)))
length(unique(DomainPair_mapping_Prot_5_all_Org(:,2)))
t=DomainPair_mapping_Prot_5_all_Org';
fid=fopen('DomainPair_mapping_Prot_5_all_Org.txt','w');
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
save DomainPair_mapping_Prot_5_all_Org_part1.mat DomainPair_mapping_Prot_5_all_Org_part1
save DomainPair_mapping_Prot_5_all_Org_part2.mat DomainPair_mapping_Prot_5_all_Org_part2
save DomainPair_mapping_Prot_5_all_Org_part3.mat DomainPair_mapping_Prot_5_all_Org_part3
save DomainPair_mapping_Prot_5_all_Org_part4.mat DomainPair_mapping_Prot_5_all_Org_part4
save DomainPair_mapping_Prot_5_all_Org_part5.mat DomainPair_mapping_Prot_5_all_Org_part5
save DomainPair_mapping_Prot_5_all_Org_part6.mat DomainPair_mapping_Prot_5_all_Org_part6
save DomainPair_mapping_Prot_5_all_Org_part7.mat DomainPair_mapping_Prot_5_all_Org_part7
save DomainPair_mapping_Prot_5_all_Org_part8.mat DomainPair_mapping_Prot_5_all_Org_part8

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
save DomainPair_mapping_Prot_5_all_MulFil.mat DomainPair_mapping_Prot_5_all_MulFil
length(unique(DomainPair_mapping_Prot_5_all_MulFil(:,1)))
length(unique(DomainPair_mapping_Prot_5_all_MulFil(:,2)))



% Calculating the single compound (non-similarity based) domain pair assocation scores for all compounds:

load Combined_Org_unique_act_inact_Prot_noSim_3_var.mat
load InterPro_entriesmatched_onlytarget_Org_bagofdom_onlydompairs.mat
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
save DomainPair_mapping_Prot_noSim_3.mat DomainPair_mapping_Prot_noSim_3

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
del_ind=find(cellfun(@isempty,DomainPair_mapping_Prot_noSim_3_all_Org(:,1)));
DomainPair_mapping_Prot_noSim_3_all_Org(del_ind,:)=[];
length(unique(DomainPair_mapping_Prot_noSim_3_all_Org(:,1)))
length(unique(DomainPair_mapping_Prot_noSim_3_all_Org(:,2)))
t=DomainPair_mapping_Prot_noSim_3_all_Org';
fid=fopen('DomainPair_mapping_Prot_noSim_3_all_Org.txt','w');
fprintf(fid,'%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n',t{:});
fclose(fid);
save DomainPair_mapping_Prot_noSim_3_all_Org.mat DomainPair_mapping_Prot_noSim_3_all_Org

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
save DomainPair_mapping_Prot_noSim_3_all_MulFil.mat DomainPair_mapping_Prot_noSim_3_all_MulFil
length(unique(DomainPair_mapping_Prot_noSim_3_all_MulFil(:,1)))
length(unique(DomainPair_mapping_Prot_noSim_3_all_MulFil(:,2)))



% Merging compound similarity based domain pair associations with noSim associations:

load DomainPair_mapping_Prot_5_all_MulFil.mat
load DomainPair_mapping_Prot_noSim_3_all_MulFil.mat
% (extracting mutual associations and selecting maximum performance according to MCC for merging)
DomainPair_mapping_Prot_5_all_MulFil_23=strcat(DomainPair_mapping_Prot_5_all_MulFil(:,1), {' '}, DomainPair_mapping_Prot_5_all_MulFil(:,2));
DomainPair_mapping_Prot_noSim_3_all_MulFil_23=strcat(DomainPair_mapping_Prot_noSim_3_all_MulFil(:,1), {' '}, DomainPair_mapping_Prot_noSim_3_all_MulFil(:,2));
[IntSec,I1,I2]=intersect(DomainPair_mapping_Prot_5_all_MulFil_23,DomainPair_mapping_Prot_noSim_3_all_MulFil_23);

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
   [Lia,Locb]=ismember(dom_temp_rev,dom_temp);
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
save DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged.mat DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged
t=DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged';
fid=fopen('DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged_finalized.txt','w');
fprintf(fid,'%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n',t{:});
fclose(fid);
length(DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged)
length(unique(DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged(:,1)))
length(unique(DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged(:,2)))
% (final domain pair association stats: 3,721 mappings between 1,456 compounds and 270 InterPro domain pairs)



% Filtering single domain and domain pair mappings by comparing and removing the corresponding single or pair mapping to the same compound when the performance is inferior:

load Domain_mapping_Prot_Sim5_noSim3_MulFil_merged.mat
load DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged.mat
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
        [Lia,Locb]=ismember(Dom_temp_on,DomPair_temp_on);
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
save Domain_mapping_Prot_Sim5_noSim3_MulFil_merged_final.mat Domain_mapping_Prot_Sim5_noSim3_MulFil_merged_final
length(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged_final)
length(unique(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged_final(:,1)))
length(unique(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged_final(:,2)))
% (final single domain association stats: 27,021 mappings between 8,164 compounds and 250 InterPro domains)

DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged_final=DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged;
DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged_final(del_realind2,:)=[];
save DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged_final.mat DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged_final
length(DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged_final)
length(unique(DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged_final(:,1)))
length(unique(DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged_final(:,2)))
% (final domain pair association stats: 22 mappings between 10 compounds and 12 InterPro domain pairs)


% following inspection of the cases where domain pair mappings score better 
% compared to single domain mappings, in 41 mapping, the performance of
% domain pair mapping is higher compared to the one domain mapping. 
% Just an example case can be given as: 
% compound: CHEMBL3621867 mapped to IPR003349 alone yields:
% TP:3, FN:0, FP:1, TN:2,
% Rec:1, Pre:0.75, Acc:0.83, F1:0.86, MCC:0.71; 
% compound: CHEMBL3621867 mapped to IPR001965 alone yields:
% TP:3, FN:0, FP:1, TN:2,
% Rec:1, Pre:0.75, Acc:0.83, F1:0.86, MCC:0.71;
% and mapped to IPR003349-IPR001965 together yields:
% TP:3, FN:0, FP:0, TN:3,
% Rec:1, Pre:1, Acc:1, F1:1, MCC:1
% also, these domains 
% alone presented in IPR003349: 47 human protein entries including 
% isoforms (10 reviewed) and IPR001965: 387 human protein entries
% including isoforms (88 reviewed), whereas together 
% IPR003349-IPR001965: 22 human protein entries (7 reviewed)
save Examples_DomPair_better.mat Examples_DomPair_better



% Predicting DTIs based on the compound-domain associations:

% (merging the compound-domain association file with the protein domain annotation file by eliminating duplicates and by selecting the pairs with the maximum performance)
load InterPro_entriesmatched_allhuman_Org_onlydom.mat
load Domain_mapping_Prot_Sim5_noSim3_MulFil_merged_final.mat
h=cell(0,3);
Dom_map_unique=unique(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged_final(:,2));
for i=1:length(Dom_map_unique)
    disp(['Domain association number: ' num2str(i), ' / ' num2str(length(Dom_map_unique))])
    [Lia,Locb]=ismember(InterPro_entriesmatched_allhuman_Org_onlydom(:,2),Dom_map_unique(i,1));
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
save DomMap_DTI_Predictions_finalized.mat DomMap_DTI_Predictions_finalized
length(unique(DomMap_DTI_Predictions_finalized(:,1)))
length(unique(DomMap_DTI_Predictions_finalized(:,2)))
t=DomMap_DTI_Predictions_finalized';
fid=fopen('DomMap_DTI_Predictions_finalized.txt','w');
fprintf(fid,'%s\t%s\t%.4f\n',t{:});
fclose(fid);

% (extracting novel DTIs from the finalized prediction set)
load Combined_Org_act_var.mat
load DomMap_DTI_Predictions_finalized.mat
Combined_Org_act_23=strcat(Combined_Org_act(:,2), {' '}, Combined_Org_act(:,1));
DomMap_DTI_Predictions_finalized_23=strcat(DomMap_DTI_Predictions_finalized(:,1), {' '}, DomMap_DTI_Predictions_finalized(:,2));
[Lia,Locb]=ismember(DomMap_DTI_Predictions_finalized_23,Combined_Org_act_23);
sum(Lia)
% (40,442 predictions were shared with the actives in the training set)
load Combined_Org_inact_withact_var.mat
Combined_Org_inact_withact_23=strcat(Combined_Org_inact_withact(:,2), {' '}, Combined_Org_inact_withact(:,1));
[Lia,Locb]=ismember(DomMap_DTI_Predictions_finalized_23,Combined_Org_inact_withact_23);
sum(Lia)
% (6,887 predictions were shared and contradicted with the inactives in the training set)

[Lia,Locb]=ismember(DomMap_DTI_Predictions_finalized_23,Combined_Org_act_23);
DomMap_DTI_Predictions_finalized_novelpred=DomMap_DTI_Predictions_finalized(Lia==0,:);
save DomMap_DTI_Predictions_finalized_novelpred.mat DomMap_DTI_Predictions_finalized_novelpred
t=DomMap_DTI_Predictions_finalized_novelpred';
fid=fopen('DomMap_DTI_Predictions_finalized_novelpred.txt','w');
fprintf(fid,'%s\t%s\t%.4f\n',t{:});
fclose(fid);
length(DomMap_DTI_Predictions_finalized_novelpred)
length(unique(DomMap_DTI_Predictions_finalized_novelpred(:,1)))
length(unique(DomMap_DTI_Predictions_finalized_novelpred(:,2)))
% (single domain mapping based final prediction stats: 3,672,076 novel DTIs between 8,158 compounds and 5,563 proteins)



% Predicting DTIs based on the compound-domain pair associations:

% (merging the compound-domain association file with the protein domain annotation file by eliminating duplicates and by selecting the pairs with the maximum performance)
load InterPro_DA_based_entriesmatched_allhuman_Org_onlydom.mat
load DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged_final.mat
h=cell(0,3);
DomPair_map_unique=unique(DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged_final(:,2));
for i=1:length(DomPair_map_unique)
    disp(['Domain association number: ' num2str(i), ' / ' num2str(length(DomPair_map_unique))])
    [Lia,Locb]=ismember(InterPro_entriesmatched_allhuman_Org_bagofdom_onlydompairs(:,2),DomPair_map_unique(i,1));
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
save DomMap_DomPair_DTI_Predictions_finalized.mat DomMap_DomPair_DTI_Predictions_finalized
length(unique(DomMap_DomPair_DTI_Predictions_finalized(:,1)))
length(unique(DomMap_DomPair_DTI_Predictions_finalized(:,2)))
t=DomMap_DomPair_DTI_Predictions_finalized';
fid=fopen('DomMap_DomPair_DTI_Predictions_finalized.txt','w');
fprintf(fid,'%s\t%s\t%.4f\n',t{:});
fclose(fid);

% (extracting novel DTIs from the finalized prediction set)
load Combined_Org_act_var.mat
load DomMap_DomPair_DTI_Predictions_finalized.mat
Combined_Org_act_23=strcat(Combined_Org_act(:,2), {' '}, Combined_Org_act(:,1));
DomMap_DomPair_DTI_Predictions_finalized_23=strcat(DomMap_DomPair_DTI_Predictions_finalized(:,1), {' '}, DomMap_DomPair_DTI_Predictions_finalized(:,2));
[Lia,Locb]=ismember(DomMap_DomPair_DTI_Predictions_finalized_23,Combined_Org_act_23);
sum(Lia)
% (27 predictions were shared with the actives in the training set)
load Combined_Org_inact_withact_var.mat
Combined_Org_inact_withact_23=strcat(Combined_Org_inact_withact(:,2), {' '}, Combined_Org_inact_withact(:,1));
[Lia,Locb]=ismember(DomMap_DomPair_DTI_Predictions_finalized_23,Combined_Org_inact_withact_23);
sum(Lia)
% (1 predictions were shared and contradicted with the inactives in the training set)

[Lia,Locb]=ismember(DomMap_DomPair_DTI_Predictions_finalized_23,Combined_Org_act_23);
DomMap_DomPair_DTI_Predictions_finalized_novelpred=DomMap_DomPair_DTI_Predictions_finalized(Lia==0,:);
save DomMap_DomPair_DTI_Predictions_finalized_novelpred.mat DomMap_DomPair_DTI_Predictions_finalized_novelpred
t=DomMap_DomPair_DTI_Predictions_finalized_novelpred';
fid=fopen('DomMap_DomPair_DTI_Predictions_finalized_novelpred.txt','w');
fprintf(fid,'%s\t%s\t%.4f\n',t{:});
fclose(fid);
length(DomMap_DomPair_DTI_Predictions_finalized_novelpred)
length(unique(DomMap_DomPair_DTI_Predictions_finalized_novelpred(:,1)))
length(unique(DomMap_DomPair_DTI_Predictions_finalized_novelpred(:,2)))
% (domain pair mapping based final prediction stats: 631 novel DTIs between 9 compounds and 286 proteins)



% Merging the single domain and domain pair association based novel DTI prediction files:

load DomMap_DTI_Predictions_finalized_novelpred.mat
load DomMap_DomPair_DTI_Predictions_finalized_novelpred.mat
DomMap_DTI_Predictions_finalized_novelpred_23=strcat(DomMap_DTI_Predictions_finalized_novelpred(:,1), {' '}, DomMap_DTI_Predictions_finalized_novelpred(:,2));
DomMap_DomPair_DTI_Predictions_finalized_novelpred_23=strcat(DomMap_DomPair_DTI_Predictions_finalized_novelpred(:,1), {' '}, DomMap_DomPair_DTI_Predictions_finalized_novelpred(:,2));
[ins,ia,ib]=intersect(DomMap_DTI_Predictions_finalized_novelpred_23,DomMap_DomPair_DTI_Predictions_finalized_novelpred_23);
find(cell2mat(DomMap_DTI_Predictions_finalized_novelpred(ia,3))>cell2mat(DomMap_DomPair_DTI_Predictions_finalized_novelpred(ib,3)))
% (none of the single domain based preditions scored better compared to domain pair mappings -most of them had the exact same score-, so all intersections can be selected from the domain pair association based predictions)
% (there is one example where the predictions score was improved due to domain pair mapping: compound id: CHEMBL450519, protein acc: P06239, score of singe domain and domain pair mappings: 0.7746 & 0.8433)
DomMap_DTI_Predictions_finalized_novelpred(ia,:)=[];
DomMap_SingDom_DomPair_merged_DTI_Pred_final_novelpred=[DomMap_DTI_Predictions_finalized_novelpred;DomMap_DomPair_DTI_Predictions_finalized_novelpred];
DomMap_SingDom_DomPair_merged_DTI_Pred_final_novelpred=sortrows(DomMap_SingDom_DomPair_merged_DTI_Pred_final_novelpred);
save DomMap_SingDom_DomPair_merged_DTI_Pred_final_novelpred.mat DomMap_SingDom_DomPair_merged_DTI_Pred_final_novelpred
t=DomMap_SingDom_DomPair_merged_DTI_Pred_final_novelpred';
fid=fopen('DomMap_SingDom_DomPair_merged_DTI_Pred_final_novelpred.txt','w');
fprintf(fid,'%s\t%s\t%.4f\n',t{:});
fclose(fid);
length(DomMap_SingDom_DomPair_merged_DTI_Pred_final_novelpred)
length(unique(DomMap_SingDom_DomPair_merged_DTI_Pred_final_novelpred(:,1)))
length(unique(DomMap_SingDom_DomPair_merged_DTI_Pred_final_novelpred(:,2)))
% (finalized merged prediction stats: 3,672,220 novel DTIs between 8,163 compounds and 5,563 proteins)


% Propagating the predictions to other compounds based on molecular similarity:

% (similarity threshold is selected as 0.8)

load Combined_Org_Sim_07.mat
load Combined_Org_Comp1_07.mat
load Combined_Org_Comp2_07.mat
ind=find(Combined_Org_Sim>=0.8);
Combined_Org_Comp1_08=Combined_Org_Comp1(ind);
Combined_Org_Comp2_08=Combined_Org_Comp2(ind);
Combined_Org_Sim_08=Combined_Org_Sim(ind);
clear Combined_Org_Comp1 Combined_Org_Comp2 Combined_Org_Sim
save Combined_Org_Comp1_08.mat Combined_Org_Comp1_08
save Combined_Org_Comp2_08.mat Combined_Org_Comp2_08
save Combined_Org_Sim_08.mat Combined_Org_Sim_08


load DomMap_SingDom_DomPair_merged_DTI_Pred_final_novelpred.mat
pred_scores_cell=DomMap_SingDom_DomPair_merged_DTI_Pred_final_novelpred(:,3);
load Combined_Org_Sim_08.mat
load Combined_Org_Comp1_08.mat
load Combined_Org_Comp2_08.mat
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
[Lia,Locb]=ismember(h_23,h_23_unique);
h_23_unique_maxscr=zeros(length(h_23_unique),1);
h_23_unique_sep=split(h_23_unique,' ');
scr_all=cell2mat(h(:,3));
for i=3776363:5082192
    disp(['Unique pair number: ' num2str(i), ' / ' num2str(length(h_23_unique_maxscr))])
    h_23_unique_maxscr(i,1)=max(scr_all(Locb==i,1));
end
DomMap_SingDom_DomPair_merged_DTI_Pred_propagat=[h_23_unique_sep num2cell(h_23_unique_maxscr)];

load Combined_Org_act_var.mat
Combined_Org_act_23=strcat(Combined_Org_act(:,2), {' '}, Combined_Org_act(:,1));
DomMap_SingDom_DomPair_merged_DTI_Pred_propagat_23=strcat(DomMap_SingDom_DomPair_merged_DTI_Pred_propagat(:,1), {' '}, DomMap_SingDom_DomPair_merged_DTI_Pred_propagat(:,2));
[Lia,Locb]=ismember(DomMap_SingDom_DomPair_merged_DTI_Pred_propagat_23,Combined_Org_act_23);
sum(Lia)
% (26,533 predictions were shared with the actives in the training set)
DomMap_SingDom_DomPair_merged_DTI_Pred_propagat=DomMap_SingDom_DomPair_merged_DTI_Pred_propagat(Lia==0,:);

load Combined_Org_inact_withact_var.mat
Combined_Org_inact_withact_23=strcat(Combined_Org_inact_withact(:,2), {' '}, Combined_Org_inact_withact(:,1));
DomMap_SingDom_DomPair_merged_DTI_Pred_propagat_23=strcat(DomMap_SingDom_DomPair_merged_DTI_Pred_propagat(:,1), {' '}, DomMap_SingDom_DomPair_merged_DTI_Pred_propagat(:,2));
[Lia,Locb]=ismember(DomMap_SingDom_DomPair_merged_DTI_Pred_propagat_23,Combined_Org_inact_withact_23);
sum(Lia)
% (4,818 predictions were shared and contradicted with the inactives in the training set)
DomMap_SingDom_DomPair_merged_DTI_Pred_propagat_novelpred_final=DomMap_SingDom_DomPair_merged_DTI_Pred_propagat(Lia==0,:);

save DomMap_SingDom_DomPair_merged_DTI_Pred_propagat_novelpred_final.mat DomMap_SingDom_DomPair_merged_DTI_Pred_propagat_novelpred_final
length(unique(DomMap_SingDom_DomPair_merged_DTI_Pred_propagat_novelpred_final(:,1)))
length(unique(DomMap_SingDom_DomPair_merged_DTI_Pred_propagat_novelpred_final(:,2)))
t=DomMap_SingDom_DomPair_merged_DTI_Pred_propagat_novelpred_final';
fid=fopen('DomMap_SingDom_DomPair_merged_DTI_Pred_propagat_novelpred_final.txt','w');
fprintf(fid,'%s\t%s\t%.4f\n',t{:});
fclose(fid);
% (finalized propagated merged prediction stats: 5,050,841 novel DTIs between 10,944 compounds and 5,461 proteins)



% Inspecting the previous LIMK predictions in the new associations:

load Domain_mapping_Prot_5_all_Org_part1.mat
load Domain_mapping_Prot_5_all_Org_part2.mat
load Domain_mapping_Prot_5_all_Org_part3.mat
Domain_mapping_Prot_5_all_Org=[Domain_mapping_Prot_5_all_Org_part1;Domain_mapping_Prot_5_all_Org_part2;Domain_mapping_Prot_5_all_Org_part3];
clear Domain_mapping_Prot_5_all_Org_part1 Domain_mapping_Prot_5_all_Org_part2 Domain_mapping_Prot_5_all_Org_part3
LIMK_kinasedom=Domain_mapping_Prot_5_all_Org(ismember(Domain_mapping_Prot_5_all_Org(:,2),'IPR000719')==1,:);
LIMK_limdom=Domain_mapping_Prot_5_all_Org(ismember(Domain_mapping_Prot_5_all_Org(:,2),'IPR001781')==1,:);
LIMK_pdzdom=Domain_mapping_Prot_5_all_Org(ismember(Domain_mapping_Prot_5_all_Org(:,2),'IPR001478')==1,:);
LIMK_ligand296=Domain_mapping_Prot_5_all_Org(ismember(Domain_mapping_Prot_5_all_Org(:,1),'CHEMBL1512352')==1,:);
LIMK_ligand314=Domain_mapping_Prot_5_all_Org(ismember(Domain_mapping_Prot_5_all_Org(:,1),'CHEMBL1316589')==1,:);
LIMK_ligand379=Domain_mapping_Prot_5_all_Org(ismember(Domain_mapping_Prot_5_all_Org(:,1),'CHEMBL518653')==1,:);
LIMK_ligand395=Domain_mapping_Prot_5_all_Org(ismember(Domain_mapping_Prot_5_all_Org(:,1),'CHEMBL516650')==1,:);
load Domain_mapping_Prot_noSim_3_all_Org.mat
LIMK_kinasedom_noSim=Domain_mapping_Prot_noSim_3_all_Org(ismember(Domain_mapping_Prot_noSim_3_all_Org(:,2),'IPR000719')==1,:);
LIMK_limdom_noSim=Domain_mapping_Prot_noSim_3_all_Org(ismember(Domain_mapping_Prot_noSim_3_all_Org(:,2),'IPR001781')==1,:);
LIMK_pdzdom_noSim=Domain_mapping_Prot_noSim_3_all_Org(ismember(Domain_mapping_Prot_noSim_3_all_Org(:,2),'IPR001478')==1,:);
LIMK_ligand296_noSim=Domain_mapping_Prot_noSim_3_all_Org(ismember(Domain_mapping_Prot_noSim_3_all_Org(:,1),'CHEMBL1512352')==1,:);
LIMK_ligand314_noSim=Domain_mapping_Prot_noSim_3_all_Org(ismember(Domain_mapping_Prot_noSim_3_all_Org(:,1),'CHEMBL1316589')==1,:);
LIMK_ligand379_noSim=Domain_mapping_Prot_noSim_3_all_Org(ismember(Domain_mapping_Prot_noSim_3_all_Org(:,1),'CHEMBL518653')==1,:);
LIMK_ligand395_noSim=Domain_mapping_Prot_noSim_3_all_Org(ismember(Domain_mapping_Prot_noSim_3_all_Org(:,1),'CHEMBL516650')==1,:);
load Domain_mapping_Prot_Sim5_noSim3_MulFil_merged.mat
LIMK_kinasedom_MulFil_merged=Domain_mapping_Prot_Sim5_noSim3_MulFil_merged(ismember(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged(:,2),'IPR000719')==1,:);
LIMK_limdom_MulFil_merged=Domain_mapping_Prot_Sim5_noSim3_MulFil_merged(ismember(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged(:,2),'IPR001781')==1,:);
LIMK_pdzdom_MulFil_merged=Domain_mapping_Prot_Sim5_noSim3_MulFil_merged(ismember(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged(:,2),'IPR001478')==1,:);
LIMK_ligand296_MulFil_merged=Domain_mapping_Prot_Sim5_noSim3_MulFil_merged(ismember(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged(:,1),'CHEMBL1512352')==1,:);
LIMK_ligand314_MulFil_merged=Domain_mapping_Prot_Sim5_noSim3_MulFil_merged(ismember(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged(:,1),'CHEMBL1316589')==1,:);
LIMK_ligand379_MulFil_merged=Domain_mapping_Prot_Sim5_noSim3_MulFil_merged(ismember(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged(:,1),'CHEMBL518653')==1,:);
LIMK_ligand395_MulFil_merged=Domain_mapping_Prot_Sim5_noSim3_MulFil_merged(ismember(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged(:,1),'CHEMBL516650')==1,:);
save LIMK_Domain_mapping_var.mat LIMK_kinasedom LIMK_limdom LIMK_pdzdom LIMK_ligand296 LIMK_ligand314 LIMK_ligand379 LIMK_ligand395 LIMK_kinasedom_noSim LIMK_limdom_noSim LIMK_pdzdom_noSim LIMK_ligand296_noSim LIMK_ligand314_noSim LIMK_ligand379_noSim LIMK_ligand395_noSim LIMK_kinasedom_MulFil_merged LIMK_limdom_MulFil_merged LIMK_pdzdom_MulFil_merged LIMK_ligand296_MulFil_merged LIMK_ligand314_MulFil_merged LIMK_ligand379_MulFil_merged LIMK_ligand395_MulFil_merged

load DomainPair_mapping_Prot_5_all_Org_part1.mat DomainPair_mapping_Prot_5_all_Org_part1
load DomainPair_mapping_Prot_5_all_Org_part2.mat DomainPair_mapping_Prot_5_all_Org_part2
load DomainPair_mapping_Prot_5_all_Org_part3.mat DomainPair_mapping_Prot_5_all_Org_part3
load DomainPair_mapping_Prot_5_all_Org_part4.mat DomainPair_mapping_Prot_5_all_Org_part4
load DomainPair_mapping_Prot_5_all_Org_part5.mat DomainPair_mapping_Prot_5_all_Org_part5
load DomainPair_mapping_Prot_5_all_Org_part6.mat DomainPair_mapping_Prot_5_all_Org_part6
load DomainPair_mapping_Prot_5_all_Org_part7.mat DomainPair_mapping_Prot_5_all_Org_part7
load DomainPair_mapping_Prot_5_all_Org_part8.mat DomainPair_mapping_Prot_5_all_Org_part8
DomainPair_mapping_Prot_5_all_Org=[DomainPair_mapping_Prot_5_all_Org_part1;DomainPair_mapping_Prot_5_all_Org_part2;DomainPair_mapping_Prot_5_all_Org_part3;DomainPair_mapping_Prot_5_all_Org_part4;DomainPair_mapping_Prot_5_all_Org_part5;DomainPair_mapping_Prot_5_all_Org_part6;DomainPair_mapping_Prot_5_all_Org_part7;DomainPair_mapping_Prot_5_all_Org_part8];
clear DomainPair_mapping_Prot_5_all_Org_part1 DomainPair_mapping_Prot_5_all_Org_part2 DomainPair_mapping_Prot_5_all_Org_part3 DomainPair_mapping_Prot_5_all_Org_part4 DomainPair_mapping_Prot_5_all_Org_part5 DomainPair_mapping_Prot_5_all_Org_part6 DomainPair_mapping_Prot_5_all_Org_part7 DomainPair_mapping_Prot_5_all_Org_part8
LIMK_kinasedom_limdom_DomPair=[DomainPair_mapping_Prot_5_all_Org(ismember(DomainPair_mapping_Prot_5_all_Org(:,2),'IPR000719-IPR001781')==1,:);DomainPair_mapping_Prot_5_all_Org(ismember(DomainPair_mapping_Prot_5_all_Org(:,2),'IPR001781-IPR000719')==1,:)];
LIMK_limdom_pdzdom_DomPair=[DomainPair_mapping_Prot_5_all_Org(ismember(DomainPair_mapping_Prot_5_all_Org(:,2),'IPR001781-IPR001478')==1,:);DomainPair_mapping_Prot_5_all_Org(ismember(DomainPair_mapping_Prot_5_all_Org(:,2),'IPR001478-IPR001781')==1,:)];
LIMK_pdzdom_kinasedom_DomPair=[DomainPair_mapping_Prot_5_all_Org(ismember(DomainPair_mapping_Prot_5_all_Org(:,2),'IPR001478-IPR000719')==1,:);DomainPair_mapping_Prot_5_all_Org(ismember(DomainPair_mapping_Prot_5_all_Org(:,2),'IPR000719-IPR001478')==1,:)];
LIMK_ligand296_DomPair=DomainPair_mapping_Prot_5_all_Org(ismember(DomainPair_mapping_Prot_5_all_Org(:,1),'CHEMBL1512352')==1,:);
LIMK_ligand314_DomPair=DomainPair_mapping_Prot_5_all_Org(ismember(DomainPair_mapping_Prot_5_all_Org(:,1),'CHEMBL1316589')==1,:);
LIMK_ligand379_DomPair=DomainPair_mapping_Prot_5_all_Org(ismember(DomainPair_mapping_Prot_5_all_Org(:,1),'CHEMBL518653')==1,:);
LIMK_ligand395_DomPair=DomainPair_mapping_Prot_5_all_Org(ismember(DomainPair_mapping_Prot_5_all_Org(:,1),'CHEMBL516650')==1,:);
load DomainPair_mapping_Prot_noSim_3_all_Org.mat
LIMK_kinasedom_limdom_noSim_DomPair=[DomainPair_mapping_Prot_noSim_3_all_Org(ismember(DomainPair_mapping_Prot_noSim_3_all_Org(:,2),'IPR000719-IPR001781')==1,:);DomainPair_mapping_Prot_noSim_3_all_Org(ismember(DomainPair_mapping_Prot_noSim_3_all_Org(:,2),'IPR001781-IPR000719')==1,:)];
LIMK_limdom_pdzdom_noSim_DomPair=[DomainPair_mapping_Prot_noSim_3_all_Org(ismember(DomainPair_mapping_Prot_noSim_3_all_Org(:,2),'IPR001781-IPR001478')==1,:);DomainPair_mapping_Prot_noSim_3_all_Org(ismember(DomainPair_mapping_Prot_noSim_3_all_Org(:,2),'IPR001478-IPR001781')==1,:)];
LIMK_pdzdom_kinasedom_noSim_DomPair=[DomainPair_mapping_Prot_noSim_3_all_Org(ismember(DomainPair_mapping_Prot_noSim_3_all_Org(:,2),'IPR001478-IPR000719')==1,:);DomainPair_mapping_Prot_noSim_3_all_Org(ismember(DomainPair_mapping_Prot_noSim_3_all_Org(:,2),'IPR000719-IPR001478')==1,:)];
LIMK_ligand296_noSim_DomPair=DomainPair_mapping_Prot_noSim_3_all_Org(ismember(DomainPair_mapping_Prot_noSim_3_all_Org(:,1),'CHEMBL1512352')==1,:);
LIMK_ligand314_noSim_DomPair=DomainPair_mapping_Prot_noSim_3_all_Org(ismember(DomainPair_mapping_Prot_noSim_3_all_Org(:,1),'CHEMBL1316589')==1,:);
LIMK_ligand379_noSim_DomPair=DomainPair_mapping_Prot_noSim_3_all_Org(ismember(DomainPair_mapping_Prot_noSim_3_all_Org(:,1),'CHEMBL518653')==1,:);
LIMK_ligand395_noSim_DomPair=DomainPair_mapping_Prot_noSim_3_all_Org(ismember(DomainPair_mapping_Prot_noSim_3_all_Org(:,1),'CHEMBL516650')==1,:);
load DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged.mat
LIMK_kinasedom_limdom_MulFil_merged_DomPair=[DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged(ismember(DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged(:,2),'IPR000719-IPR001781')==1,:);DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged(ismember(DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged(:,2),'IPR001781-IPR000719')==1,:)];
LIMK_limdom_pdzdom_MulFil_merged_DomPair=[DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged(ismember(DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged(:,2),'IPR001781-IPR001478')==1,:);DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged(ismember(DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged(:,2),'IPR001478-IPR001781')==1,:)];
LIMK_pdzdom_kinasedom_MulFil_merged_DomPair=[DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged(ismember(DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged(:,2),'IPR001478-IPR000719')==1,:);DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged(ismember(DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged(:,2),'IPR000719-IPR001478')==1,:)];
LIMK_ligand296_MulFil_merged_DomPair=DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged(ismember(DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged(:,1),'CHEMBL1512352')==1,:);
LIMK_ligand314_MulFil_merged_DomPair=DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged(ismember(DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged(:,1),'CHEMBL1316589')==1,:);
LIMK_ligand379_MulFil_merged_DomPair=DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged(ismember(DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged(:,1),'CHEMBL518653')==1,:);
LIMK_ligand395_MulFil_merged_DomPair=DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged(ismember(DomainPair_mapping_Prot_Sim5_noSim3_MulFil_merged(:,1),'CHEMBL516650')==1,:);
save LIMK_Domain_mapping_DomPair_var.mat LIMK_kinasedom_limdom_DomPair LIMK_limdom_pdzdom_DomPair LIMK_pdzdom_kinasedom_DomPair LIMK_ligand296_DomPair LIMK_ligand314_DomPair LIMK_ligand379_DomPair LIMK_ligand395_DomPair LIMK_kinasedom_limdom_noSim_DomPair LIMK_limdom_pdzdom_noSim_DomPair LIMK_pdzdom_kinasedom_noSim_DomPair LIMK_ligand296_noSim_DomPair LIMK_ligand314_noSim_DomPair LIMK_ligand379_noSim_DomPair LIMK_ligand395_noSim_DomPair LIMK_kinasedom_limdom_MulFil_merged_DomPair LIMK_limdom_pdzdom_MulFil_merged_DomPair LIMK_pdzdom_kinasedom_MulFil_merged_DomPair LIMK_ligand296_MulFil_merged_DomPair LIMK_ligand314_MulFil_merged_DomPair LIMK_ligand379_MulFil_merged_DomPair LIMK_ligand395_MulFil_merged_DomPair
% (predicted compounds is not related to LIMK protein domains at all, compound similarities will be checked to observe indirect similarities)

[Combined_Org_Comp1,Combined_Org_Comp2,Combined_Org_Sim]=textread('Combined_SimS_Org.txt', '%s %s %.4f', 'delimiter', '\t');
SimComp05_LIMK_ligand296=unique([Combined_Org_Comp2(ismember(Combined_Org_Comp1,'CHEMBL1512352')==1,1);Combined_Org_Comp1(ismember(Combined_Org_Comp2,'CHEMBL1512352')==1,1)]);
SimComp05_LIMK_ligand314=unique([Combined_Org_Comp2(ismember(Combined_Org_Comp1,'CHEMBL1316589')==1,1);Combined_Org_Comp1(ismember(Combined_Org_Comp2,'CHEMBL1316589')==1,1)]);
SimComp05_LIMK_ligand379=unique([Combined_Org_Comp2(ismember(Combined_Org_Comp1,'CHEMBL518653')==1,1);Combined_Org_Comp1(ismember(Combined_Org_Comp2,'CHEMBL518653')==1,1)]);
SimComp05_LIMK_ligand395=unique([Combined_Org_Comp2(ismember(Combined_Org_Comp1,'CHEMBL516650')==1,1);Combined_Org_Comp1(ismember(Combined_Org_Comp2,'CHEMBL516650')==1,1)]);
save SimComp05_LIMK_ligand395_var.mat SimComp05_LIMK_ligand296 SimComp05_LIMK_ligand314 SimComp05_LIMK_ligand379 SimComp05_LIMK_ligand395

load LIMK_Domain_mapping_var.mat
LIMK_kinasedom_all=[LIMK_kinasedom;LIMK_kinasedom_noSim];
LIMK_kinasedom_ligand296_SimComp05=LIMK_kinasedom_all(ismember(LIMK_kinasedom_all(:,1),SimComp05_LIMK_ligand296)==1,:);
LIMK_kinasedom_ligand314_SimComp05=LIMK_kinasedom_all(ismember(LIMK_kinasedom_all(:,1),SimComp05_LIMK_ligand314)==1,:);
LIMK_kinasedom_ligand379_SimComp05=LIMK_kinasedom_all(ismember(LIMK_kinasedom_all(:,1),SimComp05_LIMK_ligand379)==1,:);
LIMK_kinasedom_ligand395_SimComp05=LIMK_kinasedom_all(ismember(LIMK_kinasedom_all(:,1),SimComp05_LIMK_ligand395)==1,:);
LIMK_limdom_all=[LIMK_limdom;LIMK_limdom_noSim];
LIMK_limdom_ligand296_SimComp05=LIMK_limdom_all(ismember(LIMK_limdom_all(:,1),SimComp05_LIMK_ligand296)==1,:);
LIMK_limdom_ligand314_SimComp05=LIMK_limdom_all(ismember(LIMK_limdom_all(:,1),SimComp05_LIMK_ligand314)==1,:);
LIMK_limdom_ligand379_SimComp05=LIMK_limdom_all(ismember(LIMK_limdom_all(:,1),SimComp05_LIMK_ligand379)==1,:);
LIMK_limdom_ligand395_SimComp05=LIMK_limdom_all(ismember(LIMK_limdom_all(:,1),SimComp05_LIMK_ligand395)==1,:);
LIMK_pdzdom_all=[LIMK_pdzdom;LIMK_pdzdom_noSim];
LIMK_pdzdom_ligand296_SimComp05=LIMK_pdzdom_all(ismember(LIMK_pdzdom_all(:,1),SimComp05_LIMK_ligand296)==1,:);
LIMK_pdzdom_ligand314_SimComp05=LIMK_pdzdom_all(ismember(LIMK_pdzdom_all(:,1),SimComp05_LIMK_ligand314)==1,:);
LIMK_pdzdom_ligand379_SimComp05=LIMK_pdzdom_all(ismember(LIMK_pdzdom_all(:,1),SimComp05_LIMK_ligand379)==1,:);
LIMK_pdzdom_ligand395_SimComp05=LIMK_pdzdom_all(ismember(LIMK_pdzdom_all(:,1),SimComp05_LIMK_ligand395)==1,:);
save LIMK_SimComp05_SingDom_Domain_mapping_var.mat LIMK_kinasedom_all LIMK_kinasedom_ligand296_SimComp05 LIMK_kinasedom_ligand314_SimComp05 LIMK_kinasedom_ligand379_SimComp05 LIMK_kinasedom_ligand395_SimComp05 LIMK_limdom_all LIMK_limdom_ligand296_SimComp05 LIMK_limdom_ligand314_SimComp05 LIMK_limdom_ligand379_SimComp05 LIMK_limdom_ligand395_SimComp05 LIMK_pdzdom_all LIMK_pdzdom_ligand296_SimComp05 LIMK_pdzdom_ligand314_SimComp05 LIMK_pdzdom_ligand379_SimComp05 LIMK_pdzdom_ligand395_SimComp05
% (kinase domain IPR000719 - ligand296 CHEMBL1512352 relation: LIMK protein predictions are due to the similarity of ligand to compounds 5307150 and CHEMBL1443093, which are directly mapped to IPR000719 with MCC:0.12)
% (kinase domain IPR000719 - ligand314 CHEMBL1316589 relation: LIMK protein predictions are due to the similarity of ligand to compounds 648334, CHEMBL1413201 and CHEMBL1460534, which are directly mapped to LIMK's kinase domain with MCC:0.22-0.28)
% (kinase domain IPR000719 - ligand379 CHEMBL518653 relation: LIMK protein predictions are due to the similarity of ligand to compounds CHEMBL1597288, which are directly mapped to LIMK's kinase domain with MCC:0.37)
% (kinase domain IPR000719 - ligand395 CHEMBL516650 relation: LIMK protein predictions are due to the similarity of ligand to compounds CHEMBL1326188, which are directly mapped to LIMK's kinase domain with MCC:0.23)
% (no results for lim and pdz domains)

load LIMK_Domain_mapping_DomPair_var.mat
LIMK_kinasedom_limdom_DomPair_all=[LIMK_kinasedom_limdom_DomPair;LIMK_kinasedom_limdom_noSim_DomPair];
LIMK_kinasedom_limdom_DomPair_ligand296_SimComp05=LIMK_kinasedom_limdom_DomPair_all(ismember(LIMK_kinasedom_limdom_DomPair_all(:,1),SimComp05_LIMK_ligand296)==1,:);
LIMK_kinasedom_limdom_DomPair_ligand314_SimComp05=LIMK_kinasedom_limdom_DomPair_all(ismember(LIMK_kinasedom_limdom_DomPair_all(:,1),SimComp05_LIMK_ligand314)==1,:);
LIMK_kinasedom_limdom_DomPair_ligand379_SimComp05=LIMK_kinasedom_limdom_DomPair_all(ismember(LIMK_kinasedom_limdom_DomPair_all(:,1),SimComp05_LIMK_ligand379)==1,:);
LIMK_kinasedom_limdom_DomPair_ligand395_SimComp05=LIMK_kinasedom_limdom_DomPair_all(ismember(LIMK_kinasedom_limdom_DomPair_all(:,1),SimComp05_LIMK_ligand395)==1,:);
LIMK_limdom_pdzdom_DomPair_all=[LIMK_limdom_pdzdom_DomPair;LIMK_limdom_pdzdom_noSim_DomPair];
LIMK_limdom_pdzdom_DomPair_ligand296_SimComp05=LIMK_limdom_pdzdom_DomPair_all(ismember(LIMK_limdom_pdzdom_DomPair_all(:,1),SimComp05_LIMK_ligand296)==1,:);
LIMK_limdom_pdzdom_DomPair_ligand314_SimComp05=LIMK_limdom_pdzdom_DomPair_all(ismember(LIMK_limdom_pdzdom_DomPair_all(:,1),SimComp05_LIMK_ligand314)==1,:);
LIMK_limdom_pdzdom_DomPair_ligand379_SimComp05=LIMK_limdom_pdzdom_DomPair_all(ismember(LIMK_limdom_pdzdom_DomPair_all(:,1),SimComp05_LIMK_ligand379)==1,:);
LIMK_limdom_pdzdom_DomPair_ligand395_SimComp05=LIMK_limdom_pdzdom_DomPair_all(ismember(LIMK_limdom_pdzdom_DomPair_all(:,1),SimComp05_LIMK_ligand395)==1,:);
LIMK_pdzdom_kinasedom_DomPair_all=[LIMK_pdzdom_kinasedom_DomPair;LIMK_pdzdom_kinasedom_noSim_DomPair];
LIMK_pdzdom_kinasedom_DomPair_ligand296_SimComp05=LIMK_pdzdom_kinasedom_DomPair_all(ismember(LIMK_pdzdom_kinasedom_DomPair_all(:,1),SimComp05_LIMK_ligand296)==1,:);
LIMK_pdzdom_kinasedom_DomPair_ligand314_SimComp05=LIMK_pdzdom_kinasedom_DomPair_all(ismember(LIMK_pdzdom_kinasedom_DomPair_all(:,1),SimComp05_LIMK_ligand314)==1,:);
LIMK_pdzdom_kinasedom_DomPair_ligand379_SimComp05=LIMK_pdzdom_kinasedom_DomPair_all(ismember(LIMK_pdzdom_kinasedom_DomPair_all(:,1),SimComp05_LIMK_ligand379)==1,:);
LIMK_pdzdom_kinasedom_DomPair_ligand395_SimComp05=LIMK_pdzdom_kinasedom_DomPair_all(ismember(LIMK_pdzdom_kinasedom_DomPair_all(:,1),SimComp05_LIMK_ligand395)==1,:);
% (no results from the domain pair mappings)



% Domain mapping performance analysis using structurally known binding regions
% of proteins from PDB using the mappings of the InteracDome study (folder: /Performance_analysis):

[MappingPI_Pfam,MappingPI_InterPro,~]=textread('Performance_analysis/pfam_v32_interpro_mapping.txt', '%s %s %s', 'delimiter', '\t', 'bufsize', 250000);
[MappingCP_ChEMBLid,MappingCP_PDBligandid]=textread('Performance_analysis/ChEMBL_PDBligand_mapping.txt', '%s %s', 'delimiter', '\t', 'headerlines', 1);
[MappingPP_PDBligandid,MappingPP_PubChemid]=textread('Performance_analysis/PubChem_PDBligand_mapping.txt', '%s %s', 'delimiter', '\t', 'headerlines', 1);

[InteracDome_repNR_Pfamid,~,InteracDome_repNR_PDBligid,~,~,~,~,~,~,~,~,~]=textread('Performance_analysis/InteracDome_v0.3-representableNR.tsv', '%s %s %s %s %s %s %s %s %s %s %s %s', 'delimiter', '\t', 'headerlines', 5, 'bufsize', 250000);
InteracDome_repNR_Pfamid=cellfun(@(x) x(1:7), InteracDome_repNR_Pfamid, 'un', 0);
InteracDome_repNR_pairs=unique(strcat(InteracDome_repNR_Pfamid, {' '}, InteracDome_repNR_PDBligid));
[InteracDome_rep_Pfamid,~,InteracDome_rep_PDBligid,~,~,~,~,~,~,~,~,~]=textread('Performance_analysis/InteracDome_v0.3-representable.tsv', '%s %s %s %s %s %s %s %s %s %s %s %s', 'delimiter', '\t', 'headerlines', 5, 'bufsize', 250000);
InteracDome_rep_Pfamid=cellfun(@(x) x(1:7), InteracDome_rep_Pfamid, 'un', 0);
InteracDome_rep_pairs=unique(strcat(InteracDome_rep_Pfamid, {' '}, InteracDome_rep_PDBligid));
InteracDome_rep_pairs_diff=setdiff(InteracDome_rep_pairs,InteracDome_repNR_pairs);

[InterPro_v72_Domain_id,~,~]=textread('Performance_analysis/InterPro_v72_Domain_entrylist.txt', '%s %s %s', 'delimiter', '\t');
[PDB_druglike_compounds]=textread('Performance_analysis/PDB_druglike_compounds_ids.txt', '%s', 'delimiter', '\t');

% (calculating the maximum possible coverage in terms of domains and ligands)
InterPro_v72_Domain_id=unique(InterPro_v72_Domain_id);
[Lia,Locb]=ismember(MappingPI_InterPro,InterPro_v72_Domain_id);
MappingPI_InterPro_domain=MappingPI_InterPro(Lia==1);
MappingPI_Pfam_domain=MappingPI_Pfam(Lia==1);
InteracDome_repNR_Pfamid_all_unique=unique(InteracDome_repNR_Pfamid);
[Lia,Locb]=ismember(InteracDome_repNR_Pfamid_all_unique,MappingPI_Pfam_domain);
InteracDome_repNR_Pfamid_InterProDomainMap_unique=InteracDome_repNR_Pfamid_all_unique(Lia==1);
[Lia,Locb]=ismember(InteracDome_repNR_Pfamid,InteracDome_repNR_Pfamid_InterProDomainMap_unique);
InteracDome_repNR_Pfamid_InterProDomainMap=InteracDome_repNR_Pfamid(Lia==1);
InteracDome_repNR_PDBligid_InterProDomainMap=InteracDome_repNR_PDBligid(Lia==1);
MappingPP_CP_PDBligandid=unique([MappingPP_PDBligandid;MappingCP_PDBligandid]);
[Lia,Locb]=ismember(InteracDome_repNR_PDBligid_InterProDomainMap,MappingPP_CP_PDBligandid);
InteracDome_repNR_PDBligid_InterProDomainMap_sharedligand=InteracDome_repNR_PDBligid_InterProDomainMap(Lia==1);
InteracDome_repNR_Pfamid_InterProDomainMap_sharedligand=InteracDome_repNR_Pfamid_InterProDomainMap(Lia==1);
Max_Cov_dom=length(unique(InteracDome_repNR_Pfamid_InterProDomainMap_sharedligand));
Max_Cov_lig=length(unique(InteracDome_repNR_PDBligid_InterProDomainMap_sharedligand));

% (calculating the performance for the recovery of InteracDome pairs, together with coverage)
load Domain_mapping_Prot_5.mat
load Domain_mapping_Prot_noSim_3.mat

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
        [Lia,Locb]=ismember(Mapping_InteracDome_repNR_Pfamid,InteracDome_repNR_Pfamid);
        Mapping_InteracDome_repNR_final_Pfamid=Mapping_InteracDome_repNR_Pfamid(Lia==1);
        Mapping_InteracDome_repNR_final_PDBligandid=Mapping_InteracDome_repNR_PDBligandid(Lia==1);
        %(above 2 variables are the Pfam-PDBligand mappings, where Pfams and PDBligands intersect with the InteracDome Pfams and PDBligands)
        
        [Lia,Locb]=ismember(InteracDome_repNR_Pfamid,Mapping_InteracDome_repNR_final_Pfamid);
        Shared_InteracDome_repNR_Pfamid=InteracDome_repNR_Pfamid(Lia==1);
        Shared_InteracDome_repNR_PDBligid=InteracDome_repNR_PDBligid(Lia==1);
        [Lia,Locb]=ismember(Shared_InteracDome_repNR_PDBligid,Mapping_InteracDome_repNR_final_PDBligandid);
        Shared2_InteracDome_repNR_Pfamid=Shared_InteracDome_repNR_Pfamid(Lia==1);
        Shared2_InteracDome_repNR_PDBligid=Shared_InteracDome_repNR_PDBligid(Lia==1);
        %(above 2 variables are the InteracDome mappings, where Pfams and PDBligands intersect with the Pfams and PDBligands in the mappings created in this study)
        
        %(numbers, coverage and extra coverage calculations)
        Num_map=length(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged);
        Num_dom=length(unique(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged(:,2)));
        Num_lig=length(unique(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged(:,1)));
        [Lia,~]=ismember(unique(InteracDome_repNR_Pfamid_InterProDomainMap_unique),Domain_mapping_PDB_Pfam_Pfamid);
        Cov_dom=length(find(Lia==1))/Max_Cov_dom;
        [Lia,Locb]=ismember(unique(InteracDome_repNR_PDBligid_InterProDomainMap_sharedligand),Domain_mapping_PDB_PDBligandid);
        Cov_lig=length(find(Lia==1))/Max_Cov_lig;
        ExCov_dom=(Num_dom/Max_Cov_dom)-Cov_dom;
        ExCov_lig=(Num_lig/Max_Cov_lig)-Cov_lig;
        %(performance calculation)
        InteracDome_perf_cal_pairs=unique(strcat(Shared2_InteracDome_repNR_Pfamid, {' '}, Shared2_InteracDome_repNR_PDBligid));
        DomMap_perf_cal_pairs=unique(strcat(Mapping_InteracDome_repNR_final_Pfamid, {' '}, Mapping_InteracDome_repNR_final_PDBligandid));
        TP=length(intersect(InteracDome_perf_cal_pairs,DomMap_perf_cal_pairs));
        [Lia,Locb]=ismember(DomMap_perf_cal_pairs,InteracDome_rep_pairs_diff);
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
save Performance_analysis/Perf_results_DomainMapping.mat Perf_results_DomainMapping


% Calculation of the statistics for InteracDome:

length(InteracDome_repNR_Pfamid)
length(unique(InteracDome_repNR_Pfamid))
length(unique(InteracDome_repNR_PDBligid))

[Lia,Locb]=ismember(InteracDome_repNR_PDBligid,PDB_druglike_compounds);
InteracDome_repNR_PDBligid_druglikesmallmol=InteracDome_repNR_PDBligid(Lia==1);
InteracDome_repNR_Pfamid_druglikesmallmol=InteracDome_repNR_Pfamid(Lia==1);
length(InteracDome_repNR_Pfamid_druglikesmallmol)
length(unique(InteracDome_repNR_Pfamid_druglikesmallmol))
length(unique(InteracDome_repNR_PDBligid_druglikesmallmol))

save Performance_analysis/Performance_analysis_all_variables.mat


% Calculation of domain mapping statistics:

subs = findgroups(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged(:, 2));
dom_uniq = unique(Domain_mapping_Prot_Sim5_noSim3_MulFil_merged(:, 2));
dom_freq_hist=hist(subs,0.5:1:250.5);
dom_freq_hist_sort=sort(dom_freq_hist,'descend');
sum(dom_freq_hist_sort(1,1:10))/sum(dom_freq_hist_sort)
dom_most_freq_10=dom_uniq(dom_freq_hist>1080);


% Calculation of protein & family statistics in the source bioactivty datasetand the output predictions dataset:

load Combined_Org_act_inact.mat
load Combined_Org_act_var.mat
load Combined_Org_inact_withact_var.mat
load DomMap_SingDom_DomPair_merged_DTI_Pred_final_novelpred.mat
Pred_prot_acc=DomMap_SingDom_DomPair_merged_DTI_Pred_final_novelpred(:,2);

length(unique(Combined_Org_act_inact(:,1)))
length(unique(Combined_Org_act(:,1)))
length(unique(Combined_Org_inact_withact(:,1)))
length(unique(Pred_prot_acc))

[Enzyme_UniProt_acc,~,~,~]=textread('Protein_family_datasets/Enzyme_uniprot_reviewed_human.tab', '%s %s %s %s', 'headerlines', 1, 'delimiter', '\t');
[GPCR_UniProt_acc,~,~]=textread('Protein_family_datasets/GPCR_uniprot_reviewed_human.tab', '%s %s %s', 'headerlines', 1, 'delimiter', '\t');
[IonChannel_UniProt_acc,~,~]=textread('Protein_family_datasets/IonChannel_uniprot_reviewed_human.tab', '%s %s %s', 'headerlines', 1, 'delimiter', '\t');
[TranscriptionFactor_UniProt_acc,~,~]=textread('Protein_family_datasets/TranscriptionFactor_uniprot_reviewed_human.tab', '%s %s %s', 'headerlines', 1, 'delimiter', '\t');
[OtherFamilies_UniProt_acc,~,~]=textread('Protein_family_datasets/OtherFamilies_uniprot_reviewed_human.tab', '%s %s %s', 'headerlines', 1, 'delimiter', '\t');

Freq_act=zeros(5,1);
Freq_uniq_prot_act=zeros(5,1);
[Lia,Locb]=ismember(Combined_Org_act(:,1),Enzyme_UniProt_acc);
Freq_act(1,1)=sum(Lia);
Freq_uniq_prot_act(1,1)=length(unique(Combined_Org_act(Lia,1)));
[Lia,Locb]=ismember(Combined_Org_act(:,1),GPCR_UniProt_acc);
Freq_act(2,1)=sum(Lia);
Freq_uniq_prot_act(2,1)=length(unique(Combined_Org_act(Lia,1)));
[Lia,Locb]=ismember(Combined_Org_act(:,1),IonChannel_UniProt_acc);
Freq_act(3,1)=sum(Lia);
Freq_uniq_prot_act(3,1)=length(unique(Combined_Org_act(Lia,1)));
[Lia,Locb]=ismember(Combined_Org_act(:,1),TranscriptionFactor_UniProt_acc);
Freq_act(4,1)=sum(Lia);
Freq_uniq_prot_act(4,1)=length(unique(Combined_Org_act(Lia,1)));
[Lia,Locb]=ismember(Combined_Org_act(:,1),OtherFamilies_UniProt_acc);
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
[Lia,Locb]=ismember(Combined_Org_inact_withact(:,1),Enzyme_UniProt_acc);
Freq_inact(1,1)=sum(Lia);
Freq_uniq_prot_inact(1,1)=length(unique(Combined_Org_inact_withact(Lia,1)));
[Lia,Locb]=ismember(Combined_Org_inact_withact(:,1),GPCR_UniProt_acc);
Freq_inact(2,1)=sum(Lia);
Freq_uniq_prot_inact(2,1)=length(unique(Combined_Org_inact_withact(Lia,1)));
[Lia,Locb]=ismember(Combined_Org_inact_withact(:,1),IonChannel_UniProt_acc);
Freq_inact(3,1)=sum(Lia);
Freq_uniq_prot_inact(3,1)=length(unique(Combined_Org_inact_withact(Lia,1)));
[Lia,Locb]=ismember(Combined_Org_inact_withact(:,1),TranscriptionFactor_UniProt_acc);
Freq_inact(4,1)=sum(Lia);
Freq_uniq_prot_inact(4,1)=length(unique(Combined_Org_inact_withact(Lia,1)));
[Lia,Locb]=ismember(Combined_Org_inact_withact(:,1),OtherFamilies_UniProt_acc);
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
[Lia,Locb]=ismember(Pred_prot_acc,Enzyme_UniProt_acc);
Freq_pred_act(1,1)=sum(Lia);
Freq_pred_uniq_prot_act(1,1)=length(unique(Pred_prot_acc(Lia,1)));
[Lia,Locb]=ismember(Pred_prot_acc,GPCR_UniProt_acc);
Freq_pred_act(2,1)=sum(Lia);
Freq_pred_uniq_prot_act(2,1)=length(unique(Pred_prot_acc(Lia,1)));
[Lia,Locb]=ismember(Pred_prot_acc,IonChannel_UniProt_acc);
Freq_pred_act(3,1)=sum(Lia);
Freq_pred_uniq_prot_act(3,1)=length(unique(Pred_prot_acc(Lia,1)));
[Lia,Locb]=ismember(Pred_prot_acc,TranscriptionFactor_UniProt_acc);
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


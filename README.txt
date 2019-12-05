# makeOrgPackageFromEmapper

===================================================================================
###### make_orgdb_2.r ###########
===================================================================================

# if you use make_orgdb_2.r, you should open it and copy the command follow that steps

inputfiles will used in make_orgdb_2.r

(1) eggnog-mapper annotation (should remove the line with "#") result from eggnog-mapper software
(2) gene_list or gene id you want to annotation (the )
(3) cog_funclass.tab for COG annotatinon
(4) kegg_info for KEGG analysis(dowload it by parse_kegg_json.r)


=========================================================================================================
######## core_specialist_sequence.py and extract_id.py#####
========================================================================================================

The two script were used to get gene it from orthofinder results

such as：

python3 core_specialist_sequence.py -i orthogroup.txt -f Orthogroups_gene_id.txt -s genome.fasta -o gene.fasta
python3 extract_id.py -i gene.fasta -o gene_id.txt

#orthogroup.txt: gene families name
OG0000002
OG0000007
OG0000008
OG0000009
OG0000011

#Orthogroups_gene_id.txt:(the file was get from orthofinder)

OG0000000: Af_XP_747013.2 Af_XP_749620.1 Af_XP_752835.1 
OG0000001: Af_XP_747104.1 Af_XP_751088.1 Bc_XP_001546796.1 Bc_XP_001553843.1
OG0000002: Bc_XP_024553062.1 Bg_Bg1174 Bg_Bg1579 Bg_Bg4144 Cf_PHH51114.1 
OG0000003: Hc_XP_001536408.1 Hc_XP_001536571.1 Hc_XP_001536622.1 Hc_XP_001536745.1 Hc_XP_001536810.1 Hc_XP_001536869.1 

#genome.fasta:

>Af_XP_001481391.1 
MRVEQRHTYYLYEMGHSPDRVTLGHLVHGNYAMPTKSLHYSAPRMRTLRRLITLSREELESWAEIRDVRGEYGQTDKTNYGIEVNVFGAGNFGVGFVTESGKAIIAEKGRRIELLEFFETKILTDEKAMEKLRTWL
>Af_XP_001481392.1 
MPSKNTTHYILLDSVNTLSPADLERLLRAIQQRDQGMPSRVRVLTIDITKHNKPDIKAFIVEELKRKGIFQGADEDSQRRKTMVEERLWVRSNNSYLTIQQDLRKVEEIITSGGTEEELHRVLHESSTDPKELVRAEIEGLEAMLKAREIEEINELLVWVIASDKSMLLVELEAALFLRFKTVSLQPLDKKITGKYSKIFTLIYGKYLALKDHVRDCVVAQRDRPRQSADDPKITATISITNGDLKLVQRFFWDLNHYSFLGGFAFEPTSDQ
>Af_XP_001481393.1 
MPLTYKNTIPCAMIKAEQMKGPVLTEDSALEFLAINGLPGPYVYVTETLGLYAVQEFYSALEITDCVSFWLCAFSSRPGVEPVLFQGRVDSQTVTPRGSNGFAFDPIFEVQGRIRRDAQTKACQARSKNKVIYR
>Af_XP_001481394.1
MLLNISAFALSTTGVLCMLLCLVVKLFETSPVPKNIPWVLSKKGFFGRAAERFLGVNPIGFIQEGYQQYSKNQQPFVMPSVNEGDEVILPKEHSNAVLFAKESEYSFKAHISDFFQLKYTSWPLAFAERYDFFVKL




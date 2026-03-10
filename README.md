# Data Jamboree

## Scientific Questions
- How does the individual omics level change over different locations (D01-D12),
concentration levels (Control, 1x, 10x), and their interactions?
- How do the interactions between the multiple omics levels change over different
locations, concentration levels, and their interactions?
- How does the change relate to the distribution of detected individual organic
chemical compounds? Which chemical may cause the most significant adverse
outcome on the Daphnia?
- What biological/environmental insights can you obtain from your data analysis
findings?

## File Descriptions

```bash
datajamboree/
├── data/
│   │ # 
│   ├── polar_neg_pqn_imputed.csv               #
│   ├── polar_neg_pqn_imputed_glog.csv
│   ├── polar_pos_pqn_imputed.csv
│   ├── polar_pos_pqn_imputed_glog.csv
│   └── 
├── annotations/
│   ├── polar_neg_pkl_to_kegg_annotations.tsv   #
│   ├── polar_pos_pkl_to_kegg_annotations.tsv   #
│   ├── rna_dma_to_hsa_gn_mappings.tsv          #
│   └── rna_dma_to_hsa_mappings.tsv             #
├── scripts/
│   ├── fileio.R
│   ├── id_mapping.R
│   ├── irf_demo.R
│   └── wgnca_demo.R
└── README.md
```
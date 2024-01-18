def get_markers_tcell(adata=None):
    return check_markers(
        {
            "general": [
                # "CD3D",
                "CD3E",
                "CD4",
                "CD8A",
                # "CD8B",
            ],
            "Trm": [
                "CD69",
                "CXCR6",
                # "ITGAE",
                "RGS1",
            ],
            "eff.": [
                "CCR6",
                "CXCR3",
                "IFNG",
                "TNF",
            ],
            # "cytokines": [
            #     "IFNG",
            #     "TNF",
            #     # "IL2",
            #     # "IL4",
            #     # "IL17A",
            #     # "IL17F",
            #     # "IL21",
            #     # "IL22",
            # ],
            "naive/Tcm": [
                "LEF1",
                # "TCF7",
                # "LTB",
                "CCR7",
                "SELL",
                # "KLF2",
            ],
            # "naive": [
            #     "LEF1",
            #     "TCF7",
            #     "LTB",
            # ],
            # "Tcm": [
            #     "CCR7",
            #     "SELL",
            #     "KLF2",
            #     "S1PR1",
            # ],
            "Treg": [
                "FOXP3",
                "IL2RA",
                # "CTLA4",
                # "IKZF2",
                # "TIGIT",
            ],
            "gdT": [
                "TRDV2",
                "TRGV9",
                # "TRGC1",
                # "TRGC2",
                # "TRDC",
            ],
            "MAIT": [
                "TRAV1-2",
                # "KLRB1",
            ],
            "cytotoxic": [
                # "GZMK",
                # "GZMA",
                "GZMB",
                # "GZMH",
                # "GNLY",
                "PRF1",
            ],
            "NK/NKT": [
                "KLRB1",
                "KLRK1",
                "NCAM1",
                "FCGR3A",
                # "B3GAT1",
                # "CD1D",
            ],
            "prolif.": [
                "STMN1",
                "MKI67",
                # "TOP2A",
            ],
            # "stress": [
            #     "LYZ",
            #     "CD14",
            #     "CD33",
            #     "ITGAM",
            #     "FCGR3A",
            #     "FCGR3B",
            #     "CD68",
            #     "CD64",
            #     "S100A8",
            #     "S100A9",
            #     "ICAM1",
            # ],
        },
        adata=adata,
    )


def get_markers_cd4(adata=None):
    return check_markers(
        {
            "general": [
                "CD4",
                "CD8A",
                # "CD8B",
            ],
            "Th17": [
                "CCR6",
                "RORC",
                "IL17A",
                "IL17F",
                "IL23R",
                # "IL22",
            ],
            "Th1": [
                "CXCR3",
                "IFNG",
                "TNF",
                "TBX21",
            ],
            "Tfh": [
                "CXCR5",
                "CD200",
                "IL21",
                "BCL6",
                # "POU2AF1",
                # "ASCL2",
                # "ID3",
                # "ICOS",
                # "ICOSLG",
            ],
        },
        adata=adata,
    )


def get_markers_cd8(adata=None):
    return check_markers(
        {
            "general": [
                "CD4",
                "CD8A",
            ],
            "Tc1": [
                "CXCR3",
                "IFNG",
                "TNF",
                "TBX21",
            ],
            "Tc17": [
                "CCR6",
                "RORC",
                "IL17A",
                "IL17F",
                "IL23R",
            ],
            "naive/Tcm": [
                "LEF1",
                "CCR7",
                "SELL",
            ],
            "NKT": [
                "KLRB1",
                # "KLRK1",
                "NCAM1",
                "FCGR3A",
            ],
        },
        adata=adata,
    )


def check_markers(marker_dict, adata=None):
    if adata is None:
        return marker_dict
    else:
        return {
            k: [x for x in v if x in adata.var_names] for k, v in marker_dict.items()
        }

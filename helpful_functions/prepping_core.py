import pandas as pd
import anndata as ad


def annotate_adata_with_tissue_core(
        adata: ad.AnnData,
        obs_key: str,
        barcode_csv_path: str,
        drop_nan: bool = False
) -> ad.AnnData:
    """
    Annotate adata[obs_key].obs with tissue core metadata using a barcode CSV.

    Parameters
    ----------
    adata : AnnData
        The input AnnData object.
    obs_key : str
        The key under `adata` where the .obs table lives (e.g. 'square_008um').
    barcode_csv_path : str
        Path to CSV file with 'Barcode' and TMA annotations (e.g., 'TMA2140_Row4_VC').
    drop_nan : bool, default False
        Whether to drop cells with no tissue annotation.

    Returns
    -------
    AnnData
        Annotated AnnData object with 'Row_num' and 'Tissue_core' in .obs.
    """
    # Load tissue core info
    barcode_df = pd.read_csv(barcode_csv_path)

    # Identify the annotation column (assumes second column is TMA annotation)
    tma_col = barcode_df.columns[1]

    # Clean whitespace
    barcode_df['Barcode'] = barcode_df['Barcode'].astype(str).str.strip()
    barcode_df[tma_col] = barcode_df[tma_col].astype(str).str.strip()

    # Extract Row_num and Tissue_core from values like "TMA2140_Row4_VC"
    extracted = barcode_df[tma_col].str.extract(r'^(TMA\d+_Row\d+)_(\w+)$')
    extracted.columns = ['Row_num', 'Tissue_core']
    barcode_df = pd.concat([barcode_df, extracted], axis=1)

    # Drop rows with missing annotation if requested
    if drop_nan:
        barcode_df = barcode_df.dropna(subset=['Row_num', 'Tissue_core'])

    # Prepare adata.obs
    obs_df = adata[obs_key].obs.copy()
    obs_df['Barcode'] = obs_df.index

    # Merge annotation

    merged = pd.merge(obs_df, barcode_df[['Barcode', 'Row_num', 'Tissue_core']],
                      how='left', left_on='Barcode', right_on='Barcode')
    # Write to adata.obs
    adata[obs_key].obs['Row_num'] = merged.set_index('Barcode').loc[adata[obs_key].obs.index, 'Row_num']
    adata[obs_key].obs['Tissue_core'] = merged.set_index('Barcode').loc[adata[obs_key].obs.index, 'Tissue_core']

    # Drop unannotated rows if needed
    if drop_nan:
        keep_mask = ~adata[obs_key].obs['Tissue_core'].isna()
        adata[obs_key] = adata[obs_key][keep_mask]

    return adata


def add_patient_condition_from_metadata(
    adata,
    obs_key,
    metadata_path,
    rownum_col='Row_num'
):
    """
    Merge patient and condition metadata into adata.obs using Row_num as key.

    Parameters
    ----------
    adata : AnnData
        Your annotated AnnData object.
    obs_key : str
        The key for the sub-AnnData in adata.
    metadata_path : str
        Path to Excel or CSV file with columns: Row_num, Patient, Condition.
    rownum_col : str
        Column in both adata.obs and metadata to merge on (default = 'Row_num').

    Returns
    -------
    AnnData
        The AnnData object with 'Patient' and 'Condition' added to .obs.
    """
    import pandas as pd

    # Load metadata
    if metadata_path.endswith('.xlsx') or metadata_path.endswith('.xls'):
        meta_df = pd.read_excel(metadata_path)
    else:
        meta_df = pd.read_csv(metadata_path)

    # Strip whitespace from merge key
    meta_df[rownum_col] = meta_df[rownum_col].astype(str).str.strip()
    adata_rownums = adata[obs_key].obs[rownum_col].astype(str).str.strip()

    # Create mapping
    row_to_patient = meta_df.set_index(rownum_col)['Patient'].to_dict()
    row_to_condition = meta_df.set_index(rownum_col)['Condition'].to_dict()

    # Map into .obs
    adata[obs_key].obs['Patient'] = adata_rownums.map(row_to_patient)
    adata[obs_key].obs['Condition'] = adata_rownums.map(row_to_condition)

    return adata

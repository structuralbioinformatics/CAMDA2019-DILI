import os, sys, re
import numpy as np
import pandas as pd

def main():
    """
    python /home/quim/PHD/Projects/camda/scripts/get_cell_painting_data.py
    """
    precalculated_features_file = '/home/quim/PHD/Projects/camda/camda_data/cell_painting_data/CAMDA_l1000_826compounds-CellPainting_images.npz'
    data = np.load(precalculated_features_file)
    interest_features = ["Area", "Compactness", "FormFactor", "MajorAxisLength", "MeanRadius", "MinorAxisLength", "Perimeter", "Solidity"]
    
    # Find the names of the features of interest
    selected_features = []
    selected_feature_indexes = []
    for index in xrange(len(data["colnames"])):
        feature = data["colnames"][index]
        fields = feature.split('_')
        place = fields[0]
        type_feature = fields[1]
        name = fields[2]
        #name = '_'.join(fields[2:])
        if place != "Cytoplasm" and name in interest_features:
            selected_features.append(feature)
            selected_feature_indexes.append(index)
    #print(selected_features)
    #print(selected_feature_indexes)

    # Create dataframe with all features
    output_features_file = '/home/quim/PHD/Projects/camda/camda_data/cell_painting_data/features_values.tsv'
    features = data["X"]
    columns=["broad_ids", "sample_ids"]+list(data["colnames"])
    features_df = pd.DataFrame(columns=columns)
    for x in xrange(len(features)):
        #if x == 10:
        #    break
        broad_id = data["broad_ids"][x]
        sample_id = data["sample_ids"][x]
        features_values = features[x]
        row=[broad_id, sample_id]+list(features_values)
        df2 = pd.DataFrame([row], columns=columns, index=[x])
        features_df = features_df.append(df2)

    # Calculate additional features
    features_df['Cells_AreaShape_AspectRatio'] = features_df['Cells_AreaShape_MajorAxisLength']/features_df['Cells_AreaShape_MinorAxisLength']
    features_df['Cytoplasm_AreaShape_AspectRatio'] = features_df['Cytoplasm_AreaShape_MajorAxisLength']/features_df['Cytoplasm_AreaShape_MinorAxisLength']
    features_df['Nuclei_AreaShape_AspectRatio'] = features_df['Nuclei_AreaShape_MajorAxisLength']/features_df['Nuclei_AreaShape_MinorAxisLength']
    features_df['Cells_AreaShape_Perimeter/Area'] = features_df['Cells_AreaShape_Perimeter']/features_df['Cells_AreaShape_Area']
    features_df['Cytoplasm_AreaShape_Perimeter/Area'] = features_df['Cytoplasm_AreaShape_Perimeter']/features_df['Cytoplasm_AreaShape_Area']
    features_df['Nuclei_AreaShape_Perimeter/Area'] = features_df['Nuclei_AreaShape_Perimeter']/features_df['Nuclei_AreaShape_Area']

    # Write dataframe
    features_df.to_csv(output_features_file, sep="\t", index=False)

    # Filter by selected features
    output_features_file_filtered = '/home/quim/PHD/Projects/camda/camda_data/cell_painting_data/features_values_filtered.tsv'
    new_columns = ['broad_ids', 'sample_ids'] + selected_features + ['Cells_AreaShape_AspectRatio', 'Nuclei_AreaShape_AspectRatio', 'Cells_AreaShape_Perimeter/Area', 'Nuclei_AreaShape_Perimeter/Area']
    features_df_filtered = features_df.loc[:,new_columns]

    # Write filtered dataframe
    features_df_filtered.to_csv(output_features_file_filtered, sep="\t", index=False)

    return


def obtain_header_fields(fields):
    """ 
    Obtain a dictionary: "field_name" => "position" 
    """
    fields_dict = {}

    for x in xrange(0, len(fields)):
        fields_dict[fields[x]] = x

    return fields_dict

def fileExist(file):
    """
    Checks if a file exists AND is a file
    """
    return os.path.exists(file) and os.path.isfile(file)

def create_directory(directory):
    """
    Checks if a directory exists and if not, creates it
    """
    try:
        os.stat(directory)
    except:
        os.mkdir(directory)
    return


if  __name__ == "__main__":
    main()

# Function Parser Refactoring Summary

## Overview
Refactored the functional classification system to work more like the location parser, using a YAML-based configuration file for better maintainability and extensibility.

## Changes Made

### 1. Created YAML Configuration File
**File:** `data/functional_classification_terms.yaml`

- Defines functional groups and their matching patterns
- Includes 25+ functional categories (Ion Channel, GPCR, Kinase, etc.)
- Includes membrane protein indicators
- Easy to update and extend without modifying code
- Well-documented with clear structure

### 2. Updated `llps_functions.py`

#### New Functions:
- `_load_functional_terms(yaml_path)`: Loads and caches YAML terms
- `parse_function_categories(function_str, protein_name_str, yaml_path)`: 
  - Similar to `parse_location()` 
  - Cleans input text (removes brackets, parentheses, etc.)
  - Matches against YAML patterns
  - Returns list of functional categories
  - Searches both Function [CC] and Protein names fields

#### Updated Functions:
- `is_membrane_protein()`: Now uses YAML patterns
- `classify_protein_function()`: Marked as deprecated, wraps new function for backward compatibility
- `add_functional_categories()`: 
  - Now uses `parse_function_categories()`
  - Adds 'Functional_Categories' column (list of categories)
  - Creates binary columns for each category (Is_Ion_Channel, etc.)
  - Provides detailed statistics on classification results

#### Added Import:
- `import yaml`: Required for loading YAML configuration

### 3. Refactored Notebook 05

**File:** `05_interactive_functional_group_networks.ipynb`

#### Added New Cell (after imports):
- Loads pLLPS data from results
- Applies functional classification using new YAML-based parser
- Creates both single-category and exploded dataframes
- Generates pLLPS lookup dictionary
- Shows detailed classification statistics

#### Updated Configuration Cell:
- Dynamically loads available groups from classified data
- Displays all groups with protein counts
- Validates selected TARGET_GROUP exists

#### Updated Data Loading Cell:
- Uses `df_exploded` (one row per category per protein)
- Handles proteins with multiple categories correctly
- Uses `.unique()` to avoid duplicate protein IDs

## Benefits

1. **Maintainability**: Classification terms in YAML are easier to update
2. **Consistency**: Works like the location parser (familiar pattern)
3. **Flexibility**: Easy to add new categories without code changes
4. **Transparency**: All matching patterns visible in one place
5. **Extensibility**: Can add metadata to categories (descriptions, colors, etc.)
6. **Multi-field Search**: Searches both function and protein name fields
7. **Backward Compatibility**: Old function still works via wrapper

## Usage Examples

### Basic Classification
```python
import llps_functions as lf

# Classify a protein
categories = lf.parse_function_categories(
    "Voltage-gated potassium channel",
    "KCNA1_HUMAN Potassium voltage-gated channel"
)
print(categories)  # ['Ion Channel']
```

### Batch Classification
```python
# Add categories to entire dataframe
df = lf.load_llps_data('data/protein_data.xlsx')
df_classified = lf.add_functional_categories(df)

# Check results
print(df_classified['Functional_Categories'].head())
print(df_classified.columns)  # Will show Is_Ion_Channel, Is_GPCR, etc.
```

### Custom YAML Path
```python
# Use a custom classification file
categories = lf.parse_function_categories(
    function_str,
    protein_name_str,
    yaml_path='custom_terms.yaml'
)
```

## Testing Recommendations

1. Run notebook 05 to verify classification works
2. Check that all functional groups are detected correctly
3. Verify protein counts match expectations
4. Compare with previous classification results if available
5. Test with different TARGET_GROUP values

## Future Enhancements

Possible additions to the YAML file:
- Category descriptions
- Category colors for visualization
- Category hierarchies (parent-child relationships)
- Alternative names/synonyms
- References/citations for patterns
- Confidence scores for patterns

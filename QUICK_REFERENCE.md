# Quick Reference - Function Parser Refactoring

## What Changed?

### Before
```python
# Classification was hardcoded in Python
def classify_protein_function(function_str, protein_name_str):
    functional_groups = {
        'Ion Channel': [r'\bion channel\b', ...],
        'GPCR': [r'\bGPCR\b', ...],
        # ... many patterns hardcoded here
    }
```

### After
```python
# Classification rules are in YAML file
def parse_function_categories(function_str, protein_name_str, yaml_path=None):
    terms = _load_functional_terms(yaml_path)
    # Uses patterns from data/functional_classification_terms.yaml
```

## Key Improvements

| Aspect | Before | After |
|--------|--------|-------|
| Configuration | Hardcoded in Python | YAML file |
| Adding Categories | Modify code | Edit YAML |
| Pattern Visibility | Scattered in code | All in one place |
| Maintainability | Low | High |
| Extensibility | Requires coding | Edit YAML file |
| Text Fields | Function only | Function + Name |

## Files You Need to Know About

1. **`data/functional_classification_terms.yaml`**
   - Contains all classification patterns
   - Edit this to modify categories
   - Add new entries here

2. **`llps_functions.py`**
   - Main functions module
   - Contains `parse_function_categories()`
   - Contains `is_membrane_protein()`
   - Contains `add_functional_categories()`

3. **`05_interactive_functional_group_networks.ipynb`**
   - Example notebook using new parser
   - Data loading cell uses new classification

## Common Tasks

### Add a New Functional Category

Edit `data/functional_classification_terms.yaml`:
```yaml
functional_groups:
  # ... existing categories ...
  
  My_New_Category:
    patterns:
      - '\bpattern1\b'
      - '\bpattern2\b'
      - '\bpattern3\b'
```

### Add a Pattern to Existing Category

Edit `data/functional_classification_terms.yaml`:
```yaml
Ion Channel:
  patterns:
    # ... existing patterns ...
    - '\bnew_pattern\b'  # Add here
```

### Test Classification Manually

```python
import llps_functions as lf

# Test single protein
categories = lf.parse_function_categories(
    "Voltage-gated potassium channel",
    "KCNA1_HUMAN"
)
print(categories)  # ['Ion Channel']
```

### Classify All Proteins in DataFrame

```python
import llps_functions as lf

df = lf.load_llps_data('data.xlsx')
df = lf.add_functional_categories(df)

# Now df has:
# - 'Functional_Categories' column (list)
# - 'Is_Ion_Channel', 'Is_GPCR', etc. columns (boolean)
```

## Understanding the Regex Patterns

Patterns in YAML use Python regex:

```yaml
'\bion channel\b'          # Matches "ion channel" as whole word
'\bsodium channel\b'       # Matches "sodium channel"
'\bvoltage-gated\b'        # Matches "voltage-gated"
'\bligand-gated\b'         # Matches "ligand-gated"
```

- `\b` = word boundary (prevents matching substrings)
- Letters match as written
- `|` = OR operator (e.g., `\b(cat|dog)\b`)

## Text Preprocessing

The parser automatically:
1. Removes `{...}` curly brackets (citations)
2. Removes `[...]` square brackets (notes)
3. Removes `(...)` parentheses
4. Converts to lowercase
5. Combines Function and Name fields

## Troubleshooting

### Pattern Not Matching?
1. Check regex syntax in YAML
2. Check word boundaries `\b` are correct
3. Test pattern manually: `re.search(pattern, text.lower())`

### Category Not Found?
1. Verify category is in YAML
2. Check pattern spelling
3. Check text contains pattern

### Performance Slow?
1. Reduce number of patterns
2. Cache results: patterns are cached after first load
3. Use specific patterns instead of general ones

## Statistics Available

After classification:
```python
df_classified = lf.add_functional_categories(df)

# See what was classified
print(df_classified['Functional_Categories'].value_counts())

# Check binary columns
print(df_classified['Is_Ion_Channel'].sum())  # Count Ion Channels
```

## YAML Structure

```yaml
functional_groups:
  Category_Name:
    patterns:
      - '\bpattern1\b'
      - '\bpattern2\b'

membrane_indicators:
  patterns:
    - '\bmembrane\b'
    - '\btransmembrane\b'
```

## Backward Compatibility

Old code still works:
```python
# Old function name still available
categories = lf.classify_protein_function(func_str, name_str)
# Returns same result as: lf.parse_function_categories(func_str, name_str)
```

## Performance Notes

- YAML is cached after first load (fast subsequent calls)
- Classification is ~O(n) per protein (n = number of patterns)
- 20,000+ proteins classified in seconds
- DataFrame operation is vectorized (fast)

## Next Steps

1. ✅ Review YAML file structure
2. ✅ Test classification on sample proteins
3. ✅ Run notebook 05 with new classification
4. ✅ Add/modify patterns as needed
5. ✅ Update other notebooks to use new parser

---

**Reference Version:** 1.0  
**Last Updated:** January 15, 2026

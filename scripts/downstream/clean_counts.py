import pandas as pd
import re
import argparse
import sys

def clean_columns(columns, fraction_type):
    cleaned_cols = []
    for col in columns:
        # featureCounts sometimes replaces '/' and '-' with '.' in column names
        # We use a flexible regex to catch the sample number
        if fraction_type == "large":
            # Target: 24248R-01-16 or 24248R.01.16 -> extract 16
            match = re.search(r'24248R[-.]\d+[-.](\d+)', col)
        elif fraction_type == "small":
            # Target: 24248FL-07-01-02 or 24248FL.07.01.02 -> extract 02
            match = re.search(r'24248FL[-.]\d+[-.]\d+[-.](\d+)', col)
        else:
            match = None
            
        if match:
            sample_num = int(match.group(1)) # int() removes leading zeros (e.g., '02' -> 2)
            cleaned_cols.append(f"NT.{sample_num}")
        else:
            # If it doesn't match the pattern (e.g., Geneid, Length), keep original
            cleaned_cols.append(col)
            
    return cleaned_cols

def main():
    parser = argparse.ArgumentParser(description="Clean RNA-seq count matrix column names.")
    parser.add_argument("--input", required=True, help="Path to raw counts file")
    parser.add_argument("--output", required=True, help="Path to save clean counts file")
    parser.add_argument("--fraction", required=True, choices=["large", "small"], help="Fraction type")
    args = parser.parse_args()

    # Read the featureCounts table (skip the first comment line usually present)
    try:
        df = pd.read_csv(args.input, sep='\t', comment='#', index_col=0)
    except Exception as e:
        print(f"Error reading {args.input}: {e}")
        sys.exit(1)

    # Clean the column names
    df.columns = clean_columns(df.columns, args.fraction)
    
    # Drop unnecessary featureCounts columns if they exist (Chr, Start, End, Strand, Length)
    cols_to_drop = ['Chr', 'Start', 'End', 'Strand', 'Length']
    df = df.drop(columns=[c for c in cols_to_drop if c in df.columns], errors='ignore')

    # Save the cleaned matrix
    df.to_csv(args.output, sep='\t')
    print(f"Successfully cleaned {args.fraction} fraction counts. Saved to {args.output}")

if __name__ == "__main__":
    main()

from pathlib import Path
from collections import defaultdict
import pandas as pd
import re

ID_NORMALIZER_PATTERNS = [
    (re.compile(r'\s+.*$'), ''),           # drop everything after first whitespace
    (re.compile(r'\|.*$'), ''),            # drop pipe-suffixed annotations (e.g., |m.123)
    (re.compile(r'\.p\d+$'), ''),          # drop protein isoform suffix like .p1, .p2
]

def normalize_transcript_id(tid: str) -> str:
    tid = tid.strip()
    for pat, repl in ID_NORMALIZER_PATTERNS:
        tid = pat.sub(repl, tid)
    return tid

def select_folders(folders):
    """Let user select folders from a numbered list"""
    if not folders:
        print("No folders found!")
        return []
    
    print("\n📂 Available Folders:")
    for i, folder in enumerate(folders, 1):
        print(f"{i:2d}. {folder}")
    
    print("\n🔹 Select folders (comma separated numbers, 'all', or 'none'):")
    print("   Example: 1,3,5 or 'all' or 'none'")
    
    while True:
        choice = input("Your selection: ").strip().lower()
        
        if choice == 'all':
            return folders.copy()
        elif choice == 'none':
            return []
        
        try:
            selected_indices = [int(x.strip()) for x in choice.split(',')]
            selected_folders = []
            
            for idx in selected_indices:
                if 1 <= idx <= len(folders):
                    selected_folders.append(folders[idx-1])
                else:
                    print(f"⚠️  Invalid number: {idx}")
            
            if selected_folders:
                return selected_folders
            else:
                print("❌ No valid selections. Please try again.")
                
        except ValueError:
            print("❌ Please enter numbers separated by commas, or 'all'/'none'")

def is_kaas_txt_file(filename):
    """Check if a .txt file contains 'KAAS' in the name (case insensitive)"""
    return filename.lower().endswith('.txt') and 'kaas' in filename.lower()

def scan_folder_for_files(folder_path, target_extensions):
    """Scan a folder for files with specific extensions and special conditions"""
    found_files = defaultdict(list)
    
    for file_path in folder_path.iterdir():
        if file_path.is_file():
            file_extension = file_path.suffix.lower()
            filename = file_path.name
            
            # Check for standard extensions
            if file_extension in target_extensions:
                found_files[file_extension].append(filename)
            
            # Special handling for KAAS .txt files
            elif file_extension == '.txt' and is_kaas_txt_file(filename):
                found_files['.txt (KAAS)'].append(filename)
            
            # Special handling for annotations files (Eggnog)
            elif file_extension == '.annotations':
                found_files['.annotations (Eggnog)'].append(filename)
    
    return found_files

def scan_selected_folders(selected_folders, base_dir, target_extensions):
    """Scan all selected folders for target file types"""
    results = {}
    
    for folder_name in selected_folders:
        folder_path = base_dir / folder_name
        print(f"\n🔍 Scanning folder: {folder_name}")
        
        if folder_path.exists() and folder_path.is_dir():
            found_files = scan_folder_for_files(folder_path, target_extensions)
            results[folder_name] = found_files
            
            # Print results for this folder
            if found_files:
                for ext, files in found_files.items():
                    print(f"   ✅ {ext}: {len(files)} files")
            else:
                print("   ❌ No target files found")
        else:
            print(f"   ⚠️  Folder not found or inaccessible: {folder_name}")
    
    return results

def print_detailed_results(results):
    """Print detailed results of the file scan"""
    print("\n" + "="*60)
    print("📊 DETAILED SCAN RESULTS")
    print("="*60)
    
    for folder_name, file_dict in results.items():
        print(f"\n📁 {folder_name}:")
        
        if not file_dict:
            print("   No target files found")
            continue
            
        total_files = sum(len(files) for files in file_dict.values())
        print(f"   Total target files: {total_files}")
        
        for ext, files in file_dict.items():
            print(f"   {ext}: {len(files)} files")
            for file in sorted(files):
                print(f"      - {file}")

def print_summary_table(results):
    """Print a summary table showing which files are found in each folder"""
    print("\n" + "="*80)
    print("📋 SUMMARY TABLE")
    print("="*80)
    
    # Get all unique file types found
    all_file_types = set()
    for file_dict in results.values():
        all_file_types.update(file_dict.keys())
    
    all_file_types = sorted(all_file_types)
    
    # Print header
    header = "Folder".ljust(20)
    for file_type in all_file_types:
        header += file_type.ljust(15)
    print(header)
    print("-" * (20 + 15 * len(all_file_types)))
    
    # Print rows
    for folder_name, file_dict in results.items():
        row = folder_name.ljust(20)
        for file_type in all_file_types:
            count = len(file_dict.get(file_type, []))
            row += ("✅" if count > 0 else "❌").ljust(15)
        print(row)

def parse_kaas_file(file_path):
    """
    Parse KAAS file and return:
      - transcript_to_knumber: {transcript: knumber}
      - knumber_to_transcripts: {knumber: [transcripts]}
    """
    transcript_to_knumber = {}
    knumber_to_transcripts = defaultdict(list)
    
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        transcript_raw = parts[0].strip()
                        knumber = parts[1].strip()
                        if knumber and knumber not in {'N/A', '-'}:
                            transcript = normalize_transcript_id(transcript_raw)
                            transcript_to_knumber[transcript] = knumber
                            knumber_to_transcripts[knumber].append(transcript)
    except Exception as e:
        print(f"Error reading KAAS file {file_path}: {e}")
    
    return transcript_to_knumber, knumber_to_transcripts

def parse_sf_file(file_path):
    """Parse salmon quant (.sf) files and extract TPM and NumReads separately"""
    tpm_data = {}
    reads_data = {}
    
    try:
        with open(file_path, 'r') as f:
            header = next(f).strip().split('\t')
            
            # Find the column indices for TPM and NumReads
            tpm_idx = header.index('TPM') if 'TPM' in header else 3
            reads_idx = header.index('NumReads') if 'NumReads' in header else 4
            
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= max(tpm_idx, reads_idx) + 1:
                    transcript = normalize_transcript_id(parts[0])
                    
                    # Extract TPM and NumReads values
                    tpm_value = parts[tpm_idx] if len(parts) > tpm_idx else '0'
                    reads_value = parts[reads_idx] if len(parts) > reads_idx else '0'
                    
                    # Store in separate dictionaries
                    tpm_data[transcript] = tpm_value
                    reads_data[transcript] = reads_value
    except Exception as e:
        print(f"Error reading .sf file {file_path}: {e}")
    
    return {'tpm': tpm_data, 'reads': reads_data}


def parse_other_file_type(file_path, file_type):
    """Parse other file types to extract transcript-indexed data"""
    transcripts_data = {}

    try:
        if file_type == '.sf':  # salmon quant files
            return parse_sf_file(file_path)
            
        elif file_type == '.pep':  # protein sequences
            current_transcript = None
            sequence = []
            with open(file_path, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        # flush previous
                        if current_transcript and sequence:
                            seq = ''.join(sequence)
                            transcripts_data[current_transcript] = {'sequence': seq, 'length': len(seq)}
                        # start new
                        current_transcript = normalize_transcript_id(line[1:].split()[0])
                        sequence = []
                    else:
                        sequence.append(line.strip())
            if current_transcript and sequence:
                seq = ''.join(sequence)
                transcripts_data[current_transcript] = {'sequence': seq, 'length': len(seq)}

        elif file_type == '.fasta':  # nucleotide sequences
            current_transcript = None
            full_entry = []  # Store the entire FASTA entry including header
            with open(file_path, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        # Save previous entry if exists
                        if current_transcript and full_entry:
                            # Join all lines of the FASTA entry
                            full_fasta = ''.join(full_entry)
                            transcripts_data[current_transcript] = {
                                'full_fasta': full_fasta,
                                'length': len(''.join(full_entry[1:]).replace('\n', ''))  # Calculate sequence length
                            }
                        # Start new entry
                        current_transcript = normalize_transcript_id(line[1:].split()[0])
                        full_entry = [line]  # Start with the header line
                    else:
                        full_entry.append(line)
            # Save the last entry
            if current_transcript and full_entry:
                full_fasta = ''.join(full_entry)
                transcripts_data[current_transcript] = {
                    'full_fasta': full_fasta,
                    'length': len(''.join(full_entry[1:]).replace('\n', ''))
                }

        elif file_type == '.annotations' or file_type == '.annotations (Eggnog)':
            # Handle EggNOG annotation files
            transcripts_data = parse_eggnog_annotations(file_path)

        else:
            # generic fallback: first token per non-comment line treated as transcript ID
            with open(file_path, 'r') as f:
                for i, line in enumerate(f):
                    if line.strip() and not line.startswith('#'):
                        first = line.split()[0].strip()
                        transcript = normalize_transcript_id(first)
                        transcripts_data[transcript] = {
                            'line': i + 1,
                            'content': line.strip()[:100] + ('...' if len(line.strip()) > 100 else '')
                        }

    except Exception as e:
        print(f"Error reading {file_type} file {file_path}: {e}")

    return transcripts_data



def parse_eggnog_annotations(file_path):
    """
    Parse EggNOG annotation files and extract structured data
    Returns: dictionary with {transcript_id: annotation_data}
    """
    annotations = {}
    
    try:
        with open(file_path, 'r') as f:
            # Skip header lines until we find the column headers
            line = f.readline()
            while line and not line.startswith('#query'):
                line = f.readline()
            
            # If we found the header line, process the columns
            if line.startswith('#query'):
                headers = line.strip().split('\t')
                
                # Process each data line
                for line in f:
                    if line.startswith('#'):
                        continue  # Skip comment lines
                    
                    parts = line.strip().split('\t')
                    if len(parts) < len(headers):
                        continue  # Skip malformed lines
                    
                    # Extract transcript ID (normalize it)
                    transcript_id = normalize_transcript_id(parts[0])
                    
                    # Create a structured annotation dictionary
                    annotation_data = {}
                    
                    # Map each column to its header
                    for i, header in enumerate(headers):
                        if i < len(parts):
                            annotation_data[header.replace('#', '')] = parts[i]
                    
                    # Parse specific fields for easier access
                    if 'eggNOG_OGs' in annotation_data:
                        # Parse taxonomic levels
                        tax_levels = {}
                        if annotation_data['eggNOG_OGs']:
                            og_parts = annotation_data['eggNOG_OGs'].split(',')
                            for og_part in og_parts:
                                if '@' in og_part:
                                    og, tax = og_part.split('@', 1)
                                    tax_id, tax_name = tax.split('|', 1) if '|' in tax else (tax, tax)
                                    tax_levels[tax_name] = og
                        annotation_data['parsed_tax_levels'] = tax_levels
                    
                    if 'GOs' in annotation_data and annotation_data['GOs']:
                        # Parse GO terms
                        go_terms = annotation_data['GOs'].split(',')
                        annotation_data['parsed_GOs'] = [go for go in go_terms if go]
                    
                    if 'KEGG_ko' in annotation_data and annotation_data['KEGG_ko']:
                        # Parse KEGG orthology
                        kegg_ko = annotation_data['KEGG_ko'].split(',')
                        annotation_data['parsed_KEGG_ko'] = [ko for ko in kegg_ko if ko]
                    
                    if 'KEGG_Pathway' in annotation_data and annotation_data['KEGG_Pathway']:
                        # Parse KEGG pathways
                        kegg_pathways = annotation_data['KEGG_Pathway'].split(',')
                        annotation_data['parsed_KEGG_Pathway'] = [pathway for pathway in kegg_pathways if pathway]
                    
                    if 'COG_category' in annotation_data and annotation_data['COG_category']:
                        # Parse COG categories
                        cog_categories = list(annotation_data['COG_category'])
                        annotation_data['parsed_COG_categories'] = cog_categories
                    
                    if 'PFAMs' in annotation_data and annotation_data['PFAMs']:
                        # Parse PFAM domains
                        pfam_domains = annotation_data['PFAMs'].split(',')
                        annotation_data['parsed_PFAMs'] = [pfam for pfam in pfam_domains if pfam]
                    
                    # Store the annotation data
                    annotations[transcript_id] = annotation_data
                    
    except Exception as e:
        print(f"Error reading EggNOG annotation file {file_path}: {e}")
    
    return annotations

def parse_other_file_type(file_path, file_type):
    """Parse other file types to extract transcript-indexed data"""
    transcripts_data = {}

    try:
        if file_type == '.sf':  # salmon quant files
            return parse_sf_file(file_path)
            
        elif file_type == '.pep':  # protein sequences
            current_transcript = None
            sequence = []
            with open(file_path, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        # flush previous
                        if current_transcript and sequence:
                            seq = ''.join(sequence)
                            transcripts_data[current_transcript] = {'sequence': seq, 'length': len(seq)}
                        # start new
                        current_transcript = normalize_transcript_id(line[1:].split()[0])
                        sequence = []
                    else:
                        sequence.append(line.strip())
            if current_transcript and sequence:
                seq = ''.join(sequence)
                transcripts_data[current_transcript] = {'sequence': seq, 'length': len(seq)}

        elif file_type == '.fasta':  # nucleotide sequences
            current_transcript = None
            full_entry = []  # Store the entire FASTA entry including header
            with open(file_path, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        # Save previous entry if exists
                        if current_transcript and full_entry:
                            # Join all lines of the FASTA entry
                            full_fasta = ''.join(full_entry)
                            transcripts_data[current_transcript] = {
                                'full_fasta': full_fasta,
                                'length': len(''.join(full_entry[1:]).replace('\n', ''))  # Calculate sequence length
                            }
                        # Start new entry
                        current_transcript = normalize_transcript_id(line[1:].split()[0])
                        full_entry = [line]  # Start with the header line
                    else:
                        full_entry.append(line)
            # Save the last entry
            if current_transcript and full_entry:
                full_fasta = ''.join(full_entry)
                transcripts_data[current_transcript] = {
                    'full_fasta': full_fasta,
                    'length': len(''.join(full_entry[1:]).replace('\n', ''))
                }

        elif file_type == '.annotations' or file_type == '.annotations (Eggnog)':
            # Handle EggNOG annotation files
            transcripts_data = parse_eggnog_annotations(file_path)

        else:
            # generic fallback: first token per non-comment line treated as transcript ID
            with open(file_path, 'r') as f:
                for i, line in enumerate(f):
                    if line.strip() and not line.startswith('#'):
                        first = line.split()[0].strip()
                        transcript = normalize_transcript_id(first)
                        transcripts_data[transcript] = {
                            'line': i + 1,
                            'content': line.strip()[:100] + ('...' if len(line.strip()) > 100 else '')
                        }

    except Exception as e:
        print(f"Error reading {file_type} file {file_path}: {e}")

    return transcripts_data


def get_kaas_files_from_results(results, base_dir):
    """Extract KAAS file paths from scan results"""
    kaas_files = {}
    
    for folder_name, file_dict in results.items():
        if '.txt (KAAS)' in file_dict:
            for kaas_filename in file_dict['.txt (KAAS)']:
                kaas_file_path = base_dir / folder_name / kaas_filename
                kaas_files[folder_name] = kaas_file_path
                break  # Take first KAAS file found in each folder
    
    return kaas_files


def create_knumber_matrix(selected_folders, kaas_files, results, base_dir, additional_file_types=None):
    """Create matrix with Knumbers as rows and Folders as columns"""
    if additional_file_types is None:
        additional_file_types = []
    
    # First, collect all knumber mappings
    knumber_data = defaultdict(dict)
    all_knums = set()
    folder_knums = {}  # Store which knumbers are in which folders
    
    for folder_name in selected_folders:
        if folder_name in kaas_files:
            transcript_to_knumber, _ = parse_kaas_file(kaas_files[folder_name])
            
            for transcript, knumber in transcript_to_knumber.items():
                knumber_data[knumber][folder_name] = transcript
                all_knums.add(knumber)
                
                # Store folder-knumber relationship
                if folder_name not in folder_knums:
                    folder_knums[folder_name] = set()
                folder_knums[folder_name].add(knumber)
    
    # Add additional file type data
    additional_data = defaultdict(lambda: defaultdict(dict))
    
    for file_type in additional_file_types:
        print(f"  Processing {file_type} files...")
        for folder_name in selected_folders:
            if folder_name in results and file_type in results[folder_name]:
                for filename in results[folder_name][file_type]:
                    file_path = base_dir / folder_name / filename
                    file_data = parse_other_file_type(file_path, file_type)
                    
                    # Link to knumbers via transcripts
                    if folder_name in kaas_files:
                        transcript_to_knumber, _ = parse_kaas_file(kaas_files[folder_name])
                        
                        # Special handling for .sf files
                        if file_type == '.sf':
                            # Process TPM data
                            for transcript, tpm_value in file_data.get('tpm', {}).items():
                                if transcript in transcript_to_knumber:
                                    knumber = transcript_to_knumber[transcript]
                                    if 'tpm' not in additional_data[file_type][knumber]:
                                        additional_data[file_type][knumber]['tpm'] = {}
                                    additional_data[file_type][knumber]['tpm'][folder_name] = tpm_value
                            
                            # Process Reads data
                            for transcript, reads_value in file_data.get('reads', {}).items():
                                if transcript in transcript_to_knumber:
                                    knumber = transcript_to_knumber[transcript]
                                    if 'reads' not in additional_data[file_type][knumber]:
                                        additional_data[file_type][knumber]['reads'] = {}
                                    additional_data[file_type][knumber]['reads'][folder_name] = reads_value
                        else:
                            # For other file types
                            for transcript, data_value in file_data.items():
                                if transcript in transcript_to_knumber:
                                    knumber = transcript_to_knumber[transcript]
                                    additional_data[file_type][knumber][folder_name] = data_value
    
    return knumber_data, additional_data, sorted(all_knums), folder_knums



def export_matrix_to_file(knumber_data, additional_data, all_knums, selected_folders, output_file):
    """Export the matrix to an Excel file with multiple sheets"""
    print(f"📊 Starting export process...")
    
    # Create main dataframe for transcript matrix
    print("  Creating main transcript matrix...")
    try:
        rows = []
        for i, knumber in enumerate(all_knums):
            if i % 500 == 0:
                print(f"    Processing knumber {i+1}/{len(all_knums)}")
            
            row = {'Knumber': knumber}
            for folder in selected_folders:
                row[folder] = knumber_data.get(knumber, {}).get(folder, '')
            rows.append(row)
        
        df_main = pd.DataFrame(rows)
        print(f"  Main matrix shape: {df_main.shape}")
        
    except Exception as e:
        print(f"❌ Error creating main DataFrame: {e}")
        import traceback
        traceback.print_exc()
        return
    
    # Create Excel file with multiple sheets
    print("  Writing to Excel file...")
    try:
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            # Write main transcript matrix
            df_main.to_excel(writer, sheet_name='Transcripts', index=False)
            print(f"    ✓ Main transcript matrix written")
            
            # Write additional data sheets
            for file_type, data_dict in additional_data.items():
                print(f"    Processing {file_type} data...")
                
                if not data_dict:
                    print(f"      No data found for {file_type}, skipping")
                    continue
                
                # Special handling for .sf files - create separate TPM and Reads sheets
                if file_type == '.sf':
                    # Create TPM sheet
                    tpm_rows = []
                    for i, knumber in enumerate(all_knums):
                        if i % 500 == 0:
                            print(f"      Processing knumber {i+1}/{len(all_knums)} for TPM")
                        
                        row = {'Knumber': knumber}
                        for folder in selected_folders:
                            # Get TPM value for this knumber and folder
                            tpm_value = data_dict.get(knumber, {}).get('tpm', {}).get(folder, '')
                            row[folder] = tpm_value
                        tpm_rows.append(row)
                    
                    if tpm_rows:
                        df_tpm = pd.DataFrame(tpm_rows)
                        df_tpm.to_excel(writer, sheet_name='TPM', index=False)
                        print(f"      ✓ TPM data written to sheet")
                    
                    # Create Reads sheet
                    reads_rows = []
                    for i, knumber in enumerate(all_knums):
                        if i % 500 == 0:
                            print(f"      Processing knumber {i+1}/{len(all_knums)} for Reads")
                        
                        row = {'Knumber': knumber}
                        for folder in selected_folders:
                            # Get Reads value for this knumber and folder
                            reads_value = data_dict.get(knumber, {}).get('reads', {}).get(folder, '')
                            row[folder] = reads_value
                        reads_rows.append(row)
                    
                    if reads_rows:
                        df_reads = pd.DataFrame(reads_rows)
                        df_reads.to_excel(writer, sheet_name='Reads', index=False)
                        print(f"      ✓ Reads data written to sheet")
                    
                    continue
                
                # Special handling for EggNOG annotations
                if file_type == '.annotations' or file_type == '.annotations (Eggnog)':
                    # Create multiple sheets for different aspects of EggNOG data
                    create_eggnog_sheets(writer, data_dict, all_knums, selected_folders, knumber_data)
                    continue
                
                # For all other file types
                add_rows = []
                for i, knumber in enumerate(all_knums):
                    if i % 500 == 0:
                        print(f"      Processing knumber {i+1}/{len(all_knums)} for {file_type}")
                    
                    row = {'Knumber': knumber}
                    for folder in selected_folders:
                        folder_data = data_dict.get(knumber, {}).get(folder, {})
                        
                        # Convert the data to a string representation
                        if isinstance(folder_data, dict) and folder_data:
                            # Special handling for FASTA files
                            if file_type == '.fasta' and 'full_fasta' in folder_data:
                                data_str = folder_data['full_fasta']
                            elif 'sequence' in folder_data and 'length' in folder_data:
                                data_str = f"Len:{folder_data['length']}"
                            else:
                                data_str = "; ".join([f"{k}:{v}" for k, v in folder_data.items()])
                        else:
                            data_str = str(folder_data) if folder_data else ''
                        
                        row[folder] = data_str
                    add_rows.append(row)
                
                # Create DataFrame and write to sheet
                if add_rows:
                    df_add = pd.DataFrame(add_rows)
                    # Create a safe sheet name
                    sheet_name = file_type.replace('.', '').replace('(', '').replace(')', '')
                    sheet_name = sheet_name[:30]  # Excel sheet name limit
                    
                    # Check if sheet name is valid
                    invalid_chars = ['\\', '/', '*', '[', ']', ':', '?']
                    for char in invalid_chars:
                        sheet_name = sheet_name.replace(char, '_')
                    
                    try:
                        df_add.to_excel(writer, sheet_name=sheet_name, index=False)
                        print(f"      ✓ {file_type} data written to sheet '{sheet_name}'")
                    except Exception as e:
                        print(f"      ❌ Error writing {file_type} to sheet: {e}")
                        # Try with a simpler sheet name
                        simple_name = f"Sheet_{len(writer.sheets)}"
                        df_add.to_excel(writer, sheet_name=simple_name, index=False)
                        print(f"      ✓ {file_type} data written to sheet '{simple_name}' as fallback")
        
        print(f"✅ Matrix successfully exported to: {output_file}")
        
    except Exception as e:
        print(f"❌ Error exporting to Excel: {e}")
        import traceback
        traceback.print_exc()
        
        # Fallback: Export to CSV files
        import os
        output_dir = output_file.replace('.xlsx', '_files')
        os.makedirs(output_dir, exist_ok=True)
        
        # Export main matrix
        main_csv = os.path.join(output_dir, 'transcript_matrix.csv')
        df_main.to_csv(main_csv, index=False)
        print(f"✓ Main matrix exported to: {main_csv}")
        
        # Export additional data
        for file_type, data_dict in additional_data.items():
            if not data_dict:
                continue
                
            # Special handling for .sf files in CSV export
            if file_type == '.sf':
                # Export TPM data
                tpm_rows = []
                for knumber in all_knums:
                    row = {'Knumber': knumber}
                    for folder in selected_folders:
                        folder_data = data_dict.get(knumber, {}).get(folder, {})
                        tpm_value = folder_data.get('tpm', '') if isinstance(folder_data, dict) else ''
                        row[folder] = tpm_value
                    tpm_rows.append(row)
                
                if tpm_rows:
                    df_tpm = pd.DataFrame(tpm_rows)
                    tpm_csv = os.path.join(output_dir, 'tpm.csv')
                    df_tpm.to_csv(tpm_csv, index=False)
                    print(f"✓ TPM data exported to: {tpm_csv}")
                
                # Export Reads data
                reads_rows = []
                for knumber in all_knums:
                    row = {'Knumber': knumber}
                    for folder in selected_folders:
                        folder_data = data_dict.get(knumber, {}).get(folder, {})
                        reads_value = folder_data.get('reads', '') if isinstance(folder_data, dict) else ''
                        row[folder] = reads_value
                    reads_rows.append(row)
                
                if reads_rows:
                    df_reads = pd.DataFrame(reads_rows)
                    reads_csv = os.path.join(output_dir, 'reads.csv')
                    df_reads.to_csv(reads_csv, index=False)
                    print(f"✓ Reads data exported to: {reads_csv}")
                
                continue
                
            # For other file types
            add_rows = []
            for knumber in all_knums:
                row = {'Knumber': knumber}
                for folder in selected_folders:
                    folder_data = data_dict.get(knumber, {}).get(folder, {})
                    
                    if isinstance(folder_data, dict) and folder_data:
                        if file_type == '.fasta' and 'full_fasta' in folder_data:
                            data_str = folder_data['full_fasta']
                        elif 'sequence' in folder_data and 'length' in folder_data:
                            data_str = f"Len:{folder_data['length']}"
                        else:
                            data_str = "; ".join([f"{k}:{v}" for k, v in folder_data.items()])
                    else:
                        data_str = str(folder_data) if folder_data else ''
                    
                    row[folder] = data_str
                add_rows.append(row)
            
            if add_rows:
                df_add = pd.DataFrame(add_rows)
                csv_file = os.path.join(output_dir, f"{file_type.replace('.', '').replace('(', '').replace(')', '')}.csv")
                df_add.to_csv(csv_file, index=False)
                print(f"✓ {file_type} data exported to: {csv_file}")
        
        print(f"✅ Data exported to CSV files in directory: {output_dir}")


def create_eggnog_sheets(writer, eggnog_data, all_knums, selected_folders, knumber_data):
    """Create multiple sheets for different aspects of EggNOG annotation data"""
    # Sheet 1: Basic annotation summary
    basic_rows = []
    for knumber in all_knums:
        row = {'Knumber': knumber}
        for folder in selected_folders:
            annotation = eggnog_data.get(knumber, {}).get(folder, {})
            if annotation:
                # Get basic information
                desc = annotation.get('Description', '')
                evalue = annotation.get('evalue', '')
                score = annotation.get('score', '')
                row[folder] = f"{desc} (evalue: {evalue}, score: {score})"
            else:
                row[folder] = ''
        basic_rows.append(row)
    
    if basic_rows:
        df_basic = pd.DataFrame(basic_rows)
        df_basic.to_excel(writer, sheet_name='EggNOG_Basic', index=False)
        print("      ✓ EggNOG basic annotations written")
    
    # Sheet 2: GO terms
    go_rows = []
    for knumber in all_knums:
        row = {'Knumber': knumber}
        for folder in selected_folders:
            annotation = eggnog_data.get(knumber, {}).get(folder, {})
            if annotation and 'parsed_GOs' in annotation:
                row[folder] = ", ".join(annotation['parsed_GOs'])
            else:
                row[folder] = ''
        go_rows.append(row)
    
    if go_rows:
        df_go = pd.DataFrame(go_rows)
        df_go.to_excel(writer, sheet_name='EggNOG_GOs', index=False)
        print("      ✓ EggNOG GO terms written")
    
    # Sheet 3: KEGG pathways
    kegg_rows = []
    for knumber in all_knums:
        row = {'Knumber': knumber}
        for folder in selected_folders:
            annotation = eggnog_data.get(knumber, {}).get(folder, {})
            if annotation and 'parsed_KEGG_Pathway' in annotation:
                row[folder] = ", ".join(annotation['parsed_KEGG_Pathway'])
            else:
                row[folder] = ''
        kegg_rows.append(row)
    
    if kegg_rows:
        df_kegg = pd.DataFrame(kegg_rows)
        df_kegg.to_excel(writer, sheet_name='EggNOG_KEGG', index=False)
        print("      ✓ EggNOG KEGG pathways written")
    
    # Sheet 4: Taxonomic distribution
    tax_rows = []
    for knumber in all_knums:
        row = {'Knumber': knumber}
        for folder in selected_folders:
            annotation = eggnog_data.get(knumber, {}).get(folder, {})
            if annotation and 'parsed_tax_levels' in annotation:
                tax_info = []
                for tax_name, og in annotation['parsed_tax_levels'].items():
                    tax_info.append(f"{tax_name}: {og}")
                row[folder] = "; ".join(tax_info)
            else:
                row[folder] = ''
        tax_rows.append(row)
    
    if tax_rows:
        df_tax = pd.DataFrame(tax_rows)
        df_tax.to_excel(writer, sheet_name='EggNOG_Taxonomy', index=False)
        print("      ✓ EggNOG taxonomy written")
    
    # Sheet 5: PFAM domains
    pfam_rows = []
    for knumber in all_knums:
        row = {'Knumber': knumber}
        for folder in selected_folders:
            annotation = eggnog_data.get(knumber, {}).get(folder, {})
            if annotation and 'parsed_PFAMs' in annotation:
                row[folder] = ", ".join(annotation['parsed_PFAMs'])
            else:
                row[folder] = ''
        pfam_rows.append(row)
    
    if pfam_rows:
        df_pfam = pd.DataFrame(pfam_rows)
        df_pfam.to_excel(writer, sheet_name='EggNOG_PFAM', index=False)
        print("      ✓ EggNOG PFAM domains written")


def second_selection_phase(results, base_dir):
    """Second phase: select folders for Knumber matrix creation"""
    available_folders = list(results.keys())
    
    print("\n" + "="*60)
    print("🎯 SECOND SELECTION PHASE - Knumber Matrix Creation")
    print("="*60)
    
    print("Folders with KAAS files available:")
    kaas_files = get_kaas_files_from_results(results, base_dir)
    for folder in available_folders:
        status = "✅" if folder in kaas_files else "❌"
        print(f" {status} {folder}")
    
    # Select folders for matrix
    selected_matrix_folders = select_folders(available_folders)
    
    if not selected_matrix_folders:
        print("❌ No folders selected for matrix creation.")
        return
    
    # Select additional file types to include
    print("\n📁 Available additional file types to include:")
    all_file_types = set()
    for folder_data in results.values():
        all_file_types.update(folder_data.keys())
    
    file_type_list = sorted([ft for ft in all_file_types if ft != '.txt (KAAS)'])
    
    for i, ft in enumerate(file_type_list, 1):
        print(f"{i:2d}. {ft}")
    
    print("\nSelect additional file types to include (comma separated, 'none' to skip):")
    selected_file_types = []
    
    while True:
        choice = input("Your selection: ").strip().lower()
        
        if choice == 'none':
            break
        
        try:
            selected_indices = [int(x.strip()) for x in choice.split(',')]
            for idx in selected_indices:
                if 1 <= idx <= len(file_type_list):
                    selected_file_types.append(file_type_list[idx-1])
            
            break
        except ValueError:
            print("❌ Please enter numbers separated by commas, or 'none'")
    
    # Create the matrix
    print(f"\n🔨 Creating Knumber matrix for {len(selected_matrix_folders)} folders...")
    
    knumber_data, additional_data, all_knums, folder_knums = create_knumber_matrix(
        selected_matrix_folders, kaas_files, results, base_dir, selected_file_types
    )
    
    # Print statistics
    print(f"📊 Found {len(all_knums)} unique Knumbers")
    for folder in selected_matrix_folders:
        count = len(folder_knums.get(folder, []))
        print(f"   {folder}: {count} Knumbers")
    
    # Show preview of additional data
    if selected_file_types and additional_data:
        print("\n📝 Additional data preview:")
        for file_type in selected_file_types:
            if file_type in additional_data:
                data_count = sum(len(folder_data) for knumber_data in additional_data[file_type].values() 
                               for folder_data in knumber_data.values())
                print(f"   {file_type}: {data_count} data points")
            else:
                print(f"   {file_type}: No data found")
    else:
        print("   No additional data selected or found")
    
    # Export to file
    output_filename = input("\nEnter output filename (e.g., knumber_matrix.xlsx): ").strip()
    if not output_filename:
        output_filename = "knumber_matrix.xlsx"
    
    if not output_filename.endswith('.xlsx'):
        output_filename += '.xlsx'
    
    export_matrix_to_file(knumber_data, additional_data, all_knums, selected_matrix_folders, output_filename)
    
    return output_filename

def main():
    """Main function to run the complete analysis"""
    # Define target file extensions (case insensitive)
    TARGET_EXTENSIONS = {
        '.sf',              # salmon quants
        '.fasta',           # Trinity Output
        '.bed', 
        '.cds', 
        '.pep',
        '.hits',
        '.seed_orthologs',
        '.domtblout',
        '.gene_trans_map'
        # Note: .txt (KAAS) and .annotations (Eggnog) are handled separately
    }
    
    base_dir = Path(__file__).parent.resolve()
    folders = [item.name for item in base_dir.iterdir() if item.is_dir()]
    sorted_folders = sorted(folders)
    
    print(f"Found {len(folders)} folders in: {base_dir}")
    print("Looking for these file types:")
    print("  .sf (salmon quants)")
    print("  .txt (with 'KAAS' in filename)")
    print("  .fasta (Trinity Output)")
    print("  .bed")
    print("  .cds")
    print("  .pep")
    print("  .annotations (Eggnog)")
    print("  .hits")
    print("  .seed_orthologs")
    print("  .domtblout")
    print("  .gene_trans_map")
    
    # Let user select folders
    selected_folders = select_folders(sorted_folders)
    
    if not selected_folders:
        print("❌ No folders selected. Exiting.")
        return
    
    print(f"\n✅ Selected {len(selected_folders)} folders for scanning:")
    for folder in selected_folders:
        print(f" - {folder}")
    
    # Scan the selected folders
    results = scan_selected_folders(selected_folders, base_dir, TARGET_EXTENSIONS)
    
    # Print summary
    print("\n" + "="*60)
    print("📈 SCAN SUMMARY")
    print("="*60)
    
    total_folders_with_files = 0
    total_files_found = 0
    extension_counts = defaultdict(int)
    
    for folder_name, file_dict in results.items():
        folder_total = sum(len(files) for files in file_dict.values())
        if folder_total > 0:
            total_folders_with_files += 1
            total_files_found += folder_total
            for ext, files in file_dict.items():
                extension_counts[ext] += len(files)
    
    print(f"Folders with target files: {total_folders_with_files}/{len(selected_folders)}")
    print(f"Total files found: {total_files_found}")
    for ext, count in sorted(extension_counts.items()):
        print(f"{ext}: {count} files")
    
    # Print summary table
    if total_files_found > 0:
        print_summary_table(results)
    
    # Ask if user wants detailed results
    if total_files_found > 0:
        show_details = input("\nShow detailed file listings? (y/n): ").strip().lower()
        if show_details == 'y':
            print_detailed_results(results)
    
    # After initial scan, offer second phase
    if total_files_found > 0:
        proceed = input("\nProceed to Knumber matrix creation? (y/n): ").strip().lower()
        if proceed == 'y':
            output_file = second_selection_phase(results, base_dir)
            print(f"\n🎉 Analysis complete! Output file: {output_file}")

if __name__ == "__main__":
    main()
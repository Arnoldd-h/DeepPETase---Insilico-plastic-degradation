#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
01_mining_ncbi.py - Genomic Data Mining Script from NCBI
================================================================================

Description:
    This script performs mining of hypothetical proteins from metagenomes
    related to plastic degradation. It connects to the NCBI Entrez API to
    search and download candidate PETase sequences.

    The search strategy focuses on:
    - Proteins of metagenomic origin (unculturable environments)
    - Labeled as "hypothetical" (not experimentally characterized)
    - Originating from environments with plastic exposure

Date: December 2025
Project: Deep-PETase-Mining - Phase 1

Usage:
    python 01_mining_ncbi.py

Dependencies:
    - Python 3.9+
    - Biopython (Bio.Entrez, Bio.SeqIO)
    - Pandas
    - tqdm (for progress bars)

Technical Notes:
    - The length filter 100-2000 aa removes fragmented peptides and
      giant proteins that would complicate downstream AI analyses.
    - It is recommended not to exceed 200-500 sequences per run to avoid
      IP blocking by NCBI.
================================================================================
"""

# ============================================================================
# IMPORTS
# ============================================================================
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from datetime import datetime
import logging
import time
import sys

# Biopython for NCBI connection
try:
    from Bio import Entrez, SeqIO
    from Bio.SeqRecord import SeqRecord
except ImportError:
    print("❌ Error: Biopython is not installed.")
    print("   Run: pip install biopython")
    sys.exit(1)

# Pandas for metadata handling
try:
    import pandas as pd
except ImportError:
    print("❌ Error: Pandas is not installed.")
    print("   Run: pip install pandas")
    sys.exit(1)

# tqdm for progress bars (optional but recommended)
try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False
    print("⚠️ Warning: tqdm is not installed. Progress bars will not be shown.")
    print("   To install: pip install tqdm")


# ============================================================================
# GLOBAL CONFIGURATION
# ============================================================================

# Project paths (using pathlib for cross-platform compatibility)
BASE_DIR: Path = Path(__file__).parent.resolve()
RAW_DATA_DIR: Path = BASE_DIR / "data" / "raw"
METADATA_DIR: Path = BASE_DIR / "data" / "metadata"
LOGS_DIR: Path = BASE_DIR / "results" / "logs"

# Output files
OUTPUT_FASTA: Path = RAW_DATA_DIR / "raw_candidates.fasta"
OUTPUT_CSV: Path = METADATA_DIR / "candidates_info.csv"

# NCBI search configuration
# This query searches for hypothetical proteins from metagenomes related to plastics
SEARCH_QUERY: str = (
    "(metagenome[Organism] OR landfill[All Fields] OR plastic degradation[All Fields]) "
    "AND hypothetical protein[Title] "
    "AND 100:2000[Sequence Length]"
)

# NCBI limits (respect to avoid blocking)
# Without API key: maximum 3 requests/second → delay 0.34s
# With API key: maximum 10 requests/second → delay 0.1s
BATCH_SIZE: int = 200  # Number of sequences per request (NCBI allows up to 500)
DELAY_BETWEEN_REQUESTS: float = 0.35  # Seconds between requests (without API key)
DELAY_WITH_API_KEY: float = 0.1  # Seconds between requests (with API key)


# ============================================================================
# LOGGING CONFIGURATION
# ============================================================================
def setup_logging() -> logging.Logger:
    """
    Configures the logging system with console and file output.
    
    Returns:
        logging.Logger: Logger configured for the script.
    
    The log file is saved in results/logs/ with a timestamp to
    maintain a history of all executions.
    """
    # Create logs directory if it doesn't exist
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    
    # Filename with timestamp
    log_filename = LOGS_DIR / f"mining_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    
    # Configure format
    log_format = '%(asctime)s - %(levelname)s - %(message)s'
    
    logging.basicConfig(
        level=logging.INFO,
        format=log_format,
        handlers=[
            logging.FileHandler(log_filename, encoding='utf-8'),
            logging.StreamHandler(sys.stdout)
        ]
    )
    
    logger = logging.getLogger(__name__)
    logger.info(f"📝 Log saved to: {log_filename}")
    
    return logger


# ============================================================================
# CONFIGURATION FUNCTIONS
# ============================================================================
def configure_entrez() -> str:
    """
    Configures the NCBI Entrez API by requesting the user's email.
    
    Returns:
        str: Email configured for requests.
    
    Note for students:
        NCBI requires a valid email to track API usage.
        This allows them to contact you if your script causes problems.
        Always use your real email when working with NCBI!
    """
    print("\n" + "="*60)
    print("⚙️  NCBI ENTREZ CONFIGURATION")
    print("="*60)
    print("\nNCBI requires an email to use its API.")
    print("This is mandatory and helps NCBI contact you if there are problems.\n")
    
    email: str = input("📧 Enter your email: ").strip()
    
    # Basic email validation
    while "@" not in email or "." not in email:
        print("❌ Invalid email. Must contain '@' and '.'")
        email = input("📧 Enter your email: ").strip()
    
    # Configure Entrez
    Entrez.email = email
    
    # Optional API key (increases request limit)
    print("\n💡 Tip: If you have an NCBI API key, you can increase the download speed.")
    api_key: str = input("🔑 API Key (press Enter to skip): ").strip()
    
    if api_key:
        Entrez.api_key = api_key
        logger.info("✅ API Key configured")
    
    logger.info(f"✅ Email configured: {email}")
    
    return email


# ============================================================================
# SEARCH FUNCTIONS
# ============================================================================
def search_candidates(query: str, max_results: int = 200) -> List[str]:
    """
    Searches for candidate proteins in NCBI using the specified query.
    
    Args:
        query (str): Search term for NCBI Protein.
        max_results (int): Maximum number of results to obtain.
    
    Returns:
        List[str]: List of NCBI IDs (accession numbers).
    
    Example:
        >>> ids = search_candidates("hypothetical protein", max_results=100)
        >>> print(f"Found: {len(ids)} IDs")
    
    Note for students:
        This function uses ESearch, which is the first step of any
        NCBI search. It only returns IDs, not the sequences.
    """
    logger.info(f"🔍 Starting search in NCBI Protein...")
    logger.info(f"📝 Query: {query}")
    logger.info(f"📊 Maximum results requested: {max_results}")
    
    try:
        # ESearch: Searches and returns list of IDs
        handle = Entrez.esearch(
            db="protein",           # Protein database
            term=query,             # Our complex query
            retmax=max_results,     # Results limit
            usehistory="y",         # Use history (more efficient for many results)
            idtype="acc"            # Return accession numbers
        )
        
        # Parse result
        search_results: Dict = Entrez.read(handle)
        handle.close()
        
        # Extract IDs
        id_list: List[str] = search_results["IdList"]
        total_count: int = int(search_results["Count"])
        
        logger.info(f"✅ Search completed")
        logger.info(f"📊 Total sequences found in NCBI: {total_count}")
        logger.info(f"📥 IDs retrieved: {len(id_list)}")
        
        if total_count > max_results:
            logger.warning(
                f"⚠️ There are {total_count} results but only {max_results} were requested. "
                f"Consider increasing max_results if you need more data."
            )
        
        return id_list
        
    except Exception as e:
        logger.error(f"❌ Error during search: {e}")
        raise


def fetch_details(id_list: List[str]) -> Tuple[List[SeqRecord], List[Dict]]:
    """
    Downloads the complete details of sequences from their IDs.
    
    Args:
        id_list (List[str]): List of NCBI IDs to download.
    
    Returns:
        Tuple[List[SeqRecord], List[Dict]]: 
            - List of SeqRecord objects with sequences
            - List of dictionaries with extracted metadata
    
    Note for students:
        This function uses EFetch to download sequences in GenBank
        format (gb), which contains rich annotations. It then extracts
        the FASTA sequence and relevant metadata.
        
        We download in batches to:
        1. Avoid timeouts on large downloads
        2. Respect NCBI limits
        3. Show progress to the user
    """
    logger.info(f"📥 Starting download of {len(id_list)} sequences...")
    
    records: List[SeqRecord] = []
    metadata: List[Dict] = []
    
    # Calculate the number of batches needed
    num_batches: int = (len(id_list) + BATCH_SIZE - 1) // BATCH_SIZE
    
    # Configure the progress bar
    if TQDM_AVAILABLE:
        batch_iterator = tqdm(range(0, len(id_list), BATCH_SIZE), 
                              desc="Downloading", 
                              unit="batch",
                              total=num_batches)
    else:
        batch_iterator = range(0, len(id_list), BATCH_SIZE)
        print("Downloading sequences (this may take a few minutes)...")
    
    for start in batch_iterator:
        end = min(start + BATCH_SIZE, len(id_list))
        batch_ids = id_list[start:end]
        
        try:
            # EFetch: Downloads records in GenBank format
            handle = Entrez.efetch(
                db="protein",
                id=batch_ids,
                rettype="gb",           # GenBank format (richer in annotations)
                retmode="text"
            )
            
            # Parse GenBank records
            batch_records = list(SeqIO.parse(handle, "genbank"))
            handle.close()
            
            # Process each record
            for record in batch_records:
                records.append(record)
                
                # Extract metadata
                meta: Dict = extract_metadata(record)
                metadata.append(meta)
            
            # Respect NCBI limits (shorter delay if API key is present)
            delay = DELAY_WITH_API_KEY if Entrez.api_key else DELAY_BETWEEN_REQUESTS
            time.sleep(delay)
            
        except Exception as e:
            logger.error(f"❌ Error downloading batch {start}-{end}: {e}")
            # Continue with next batch instead of aborting
            continue
    
    logger.info(f"✅ Download completed: {len(records)} sequences retrieved")
    
    return records, metadata


def extract_metadata(record: SeqRecord) -> Dict:
    """
    Extracts relevant metadata from a GenBank record.
    
    Args:
        record (SeqRecord): Biopython record with GenBank annotations.
    
    Returns:
        Dict: Dictionary with extracted metadata.
    
    Extracted metadata:
        - id: Accession number
        - description: Protein description
        - organism: Source organism
        - length: Length in amino acids
        - source: Record source
        - taxonomy: Taxonomic classification
        - date: Deposit date
    """
    # We try to extract each field safely
    meta: Dict = {
        "id": record.id,
        "description": record.description,
        "length": len(record.seq),
        "organism": "Unknown",
        "source": "Unknown",
        "taxonomy": "",
        "date": ""
    }
    
    # Extract annotations from the record
    annotations = record.annotations
    
    if "organism" in annotations:
        meta["organism"] = annotations["organism"]
    
    if "source" in annotations:
        meta["source"] = annotations["source"]
    
    if "taxonomy" in annotations:
        meta["taxonomy"] = "; ".join(annotations["taxonomy"])
    
    if "date" in annotations:
        meta["date"] = annotations["date"]
    
    # Try to extract information from features
    for feature in record.features:
        if feature.type == "source":
            qualifiers = feature.qualifiers
            
            # Look for isolation information
            if "isolation_source" in qualifiers:
                meta["isolation_source"] = qualifiers["isolation_source"][0]
            
            # Look for environmental information
            if "environmental_sample" in qualifiers:
                meta["environmental_sample"] = "Yes"
    
    return meta


# ============================================================================
# SAVE FUNCTIONS
# ============================================================================
def save_fasta(records: List[SeqRecord], output_path: Path) -> int:
    """
    Saves sequences in FASTA format.
    
    Args:
        records (List[SeqRecord]): List of sequences to save.
        output_path (Path): Output file path.
    
    Returns:
        int: Number of sequences saved.
    
    Note for students:
        FASTA format is the de facto standard for biological sequences.
        It's a simple format: header (>) followed by the sequence.
    """
    # Ensure the directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    try:
        count = SeqIO.write(records, output_path, "fasta")
        logger.info(f"💾 {count} sequences saved to: {output_path}")
        return count
    except Exception as e:
        logger.error(f"❌ Error saving FASTA: {e}")
        raise


def save_metadata(metadata: List[Dict], output_path: Path) -> None:
    """
    Saves metadata in CSV format using Pandas.
    
    Args:
        metadata (List[Dict]): List of dictionaries with metadata.
        output_path (Path): Output CSV file path.
    
    Note for students:
        CSV is ideal for metadata because:
        1. It can be opened in Excel
        2. It's easy to read and filter with Pandas
        3. It's a universal format
    """
    # Ensure the directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    try:
        # Create DataFrame
        df = pd.DataFrame(metadata)
        
        # Sort by length (longer proteins first)
        df = df.sort_values("length", ascending=False)
        
        # Save CSV
        df.to_csv(output_path, index=False, encoding='utf-8')
        
        logger.info(f"💾 Metadata saved to: {output_path}")
        logger.info(f"📊 Columns: {', '.join(df.columns)}")
        
    except Exception as e:
        logger.error(f"❌ Error saving CSV: {e}")
        raise


# ============================================================================
# SUMMARY FUNCTION
# ============================================================================
def print_summary(records: List[SeqRecord], metadata: List[Dict]) -> None:
    """
    Prints a statistical summary of the downloaded data.
    
    Args:
        records (List[SeqRecord]): Downloaded sequences.
        metadata (List[Dict]): Extracted metadata.
    """
    df = pd.DataFrame(metadata)
    
    print("\n" + "="*60)
    print("📊 MINING SUMMARY")
    print("="*60)
    
    print(f"\n📈 General statistics:")
    print(f"   • Total sequences: {len(records)}")
    print(f"   • Average length: {df['length'].mean():.1f} aa")
    print(f"   • Minimum length: {df['length'].min()} aa")
    print(f"   • Maximum length: {df['length'].max()} aa")
    print(f"   • Standard deviation: {df['length'].std():.1f} aa")
    
    print(f"\n🦠 Top 5 organisms/sources:")
    top_organisms = df['organism'].value_counts().head(5)
    for i, (org, count) in enumerate(top_organisms.items(), 1):
        print(f"   {i}. {org}: {count} sequences")
    
    print("\n" + "="*60)
    print(f"📁 Generated files:")
    print(f"   • FASTA: {OUTPUT_FASTA}")
    print(f"   • CSV:   {OUTPUT_CSV}")
    print("="*60)


# ============================================================================
# MAIN FUNCTION
# ============================================================================
def main() -> None:
    """
    Main function that orchestrates the entire mining process.
    
    Flow:
    1. Configure Entrez with user's email
    2. Search for candidate protein IDs
    3. Download sequences and metadata
    4. Save FASTA and CSV files
    5. Display statistical summary
    """
    print("\n" + "🧬"*30)
    print("   DEEP-PETASE-MINING - Phase 1: NCBI Mining")
    print("🧬"*30)
    
    # Verify that folders exist
    if not RAW_DATA_DIR.exists() or not METADATA_DIR.exists():
        logger.error("❌ Project folders do not exist.")
        logger.error("   Run first: python 00_setup_project.py")
        sys.exit(1)
    
    # Step 1: Configure Entrez
    configure_entrez()
    
    # Step 2: First search to know how many sequences are available
    print("\n" + "-"*60)
    print("📊 SEARCH CONFIGURATION")
    print("-"*60)
    
    print(f"\n🔍 Search query:")
    print(f"   {SEARCH_QUERY}")
    
    # First query how many sequences are available
    logger.info("🔎 Querying number of sequences available in NCBI...")
    try:
        handle = Entrez.esearch(db="protein", term=SEARCH_QUERY, retmax=0)
        search_info = Entrez.read(handle)
        handle.close()
        total_available: int = int(search_info["Count"])
        logger.info(f"📊 Sequences available in NCBI: {total_available:,}")
    except Exception as e:
        logger.error(f"❌ Error querying NCBI: {e}")
        total_available = 0
    
    print(f"\n📊 Sequences available matching the query: {total_available:,}")
    
    # Ask how many to download
    try:
        print(f"\nHow many sequences do you want to download?")
        print(f"   - Enter a number (e.g., 500)")
        print(f"   - Type 'all' to download all {total_available:,} available")
        print(f"   - Press Enter to use the default value [200]")
        
        max_results_input = input("\n👉 Your choice: ").strip().lower()
        
        if max_results_input in ['todas', 'all', 'todo', 'max']:
            max_results: int = total_available
            logger.info(f"📥 Downloading ALL available sequences: {max_results:,}")
        elif max_results_input == '':
            max_results = 200
        else:
            max_results = int(max_results_input)
        
        # Only warn if it's a very large number, but do NOT limit
        if max_results > 5000:
            logger.warning(f"⚠️ You are about to download {max_results:,} sequences. This may take a long time.")
            estimated_time = (max_results / BATCH_SIZE) * DELAY_BETWEEN_REQUESTS / 60
            print(f"   ⏱️ Estimated time: ~{estimated_time:.1f} minutes")
            confirm_large = input("   Continue? [y/N]: ").strip().lower()
            if confirm_large not in ['s', 'si', 'sí', 'y', 'yes']:
                max_results = 1000
                logger.info("📉 Reduced to 1000 sequences by user choice")
            
    except ValueError:
        logger.warning("⚠️ Invalid input, using default value: 200")
        max_results = 200
    
    # Confirmation
    print("\n" + "-"*60)
    confirm = input("Proceed with download? [y/N]: ").strip().lower()
    if confirm not in ['s', 'si', 'sí', 'y', 'yes']:
        print("❌ Operation cancelled by user.")
        sys.exit(0)
    
    # Step 3: Search IDs
    print("\n" + "-"*60)
    id_list = search_candidates(SEARCH_QUERY, max_results)
    
    if not id_list:
        logger.error("❌ No sequences found. Check the query.")
        sys.exit(1)
    
    # Step 4: Download details
    print("\n" + "-"*60)
    records, metadata = fetch_details(id_list)
    
    if not records:
        logger.error("❌ Could not download sequences.")
        sys.exit(1)
    
    # Step 5: Save files
    print("\n" + "-"*60)
    logger.info("💾 Saving files...")
    save_fasta(records, OUTPUT_FASTA)
    save_metadata(metadata, OUTPUT_CSV)
    
    # Step 6: Summary
    print_summary(records, metadata)
    
    # Final message
    print("\n✨ Mining completed successfully!")
    print("👉 Next step: python 02_exploratory_analysis.py")
    print("\n" + "🧬"*30 + "\n")


# ============================================================================
# ENTRY POINT
# ============================================================================
if __name__ == "__main__":
    # Configure logging
    logger = setup_logging()
    
    try:
        main()
    except KeyboardInterrupt:
        print("\n\n⚠️ Execution cancelled by user.")
        logger.info("Execution cancelled by user")
    except Exception as e:
        logger.error(f"❌ Fatal error: {e}")
        raise

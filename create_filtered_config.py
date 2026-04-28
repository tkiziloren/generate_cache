#!/usr/bin/env python3
"""
pdbbind_optimized.yml'den failed_cases.txt'deki proteinleri çıkararak yeni config oluştur.
"""

import re
from pathlib import Path

YAML_FILE = "/Users/tevfik/Sandbox/Tevfik/Projects/phd_examples/3dunet-apbs/3dunet/config/codon/pdbbind_optimized.yml"
FAILED_CASES = "/Users/tevfik/Sandbox/github/PHD/data/pdbbind/refined-set-only-used-in-codon-tests/box161/failed_cases.txt"
OUTPUT_FILE = "/Users/tevfik/Sandbox/Tevfik/Projects/phd_examples/generate_cache/pdbbind_optimized_filtered.yml"

def read_failed_cases(filepath):
    """Failed cases dosyasını oku."""
    with open(filepath, 'r') as f:
        failed = set(line.strip() for line in f if line.strip())
    print(f"Failed cases: {len(failed)} proteins")
    return failed

def filter_yaml_config(yaml_path, failed_set, output_path):
    """YAML'i oku, failed cases'leri çıkar, yeni YAML yaz."""
    
    with open(yaml_path, 'r') as f:
        content = f.read()
    
    # Train, validation, test bölümlerini bul
    def extract_and_filter_section(section_name, content, failed_set):
        """Bir section'daki protein listesini çıkar ve filtrele."""
        pattern = rf'{section_name}:\s*\n((?:\s+-\s+\w+\n)+)'
        match = re.search(pattern, content, re.MULTILINE)
        
        if not match:
            return [], content
        
        section_text = match.group(1)
        proteins = re.findall(r'-\s+(\w+)', section_text)
        
        # Failed olanları çıkar
        filtered = [p for p in proteins if p not in failed_set]
        
        print(f"\n{section_name.capitalize()}:")
        print(f"  Orijinal: {len(proteins)}")
        print(f"  Başarısız: {len(proteins) - len(filtered)}")
        print(f"  Kalan: {len(filtered)}")
        
        # Yeni listeyi oluştur
        new_section = f"{section_name}:\n"
        for p in filtered:
            new_section += f"    - {p}\n"
        
        # Content'te replace et
        new_content = content.replace(match.group(0), new_section)
        
        return filtered, new_content
    
    # Her section'ı filtrele
    train_filtered, content = extract_and_filter_section("train", content, failed_set)
    val_filtered, content = extract_and_filter_section("validation", content, failed_set)
    test_filtered, content = extract_and_filter_section("test", content, failed_set)
    
    # Yeni dosyayı yaz
    with open(output_path, 'w') as f:
        f.write(content)
    
    print(f"\n{'='*60}")
    print(f"Toplam özet:")
    print(f"  Train: {len(train_filtered)}")
    print(f"  Validation: {len(val_filtered)}")
    print(f"  Test: {len(test_filtered)}")
    print(f"  Toplam: {len(train_filtered) + len(val_filtered) + len(test_filtered)}")
    print(f"\nYeni config dosyası: {output_path}")
    print(f"{'='*60}")

if __name__ == "__main__":
    # Failed cases'leri oku
    failed = read_failed_cases(FAILED_CASES)
    
    # YAML'i filtrele ve yaz
    filter_yaml_config(YAML_FILE, failed, OUTPUT_FILE)

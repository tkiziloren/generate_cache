#!/usr/bin/env python3
"""
YAML'deki protein listelerini okuyup kaynak dizinden hedef dizine kopyalayan script.
"""

import os
import shutil
import re
from pathlib import Path

# Dizin ayarları
YAML_FILE = "/Users/tevfik/Sandbox/Tevfik/Projects/phd_examples/3dunet-apbs/3dunet/config/codon/pdbbind_optimized.yml"
SOURCE_DIR = "/Users/tevfik/Sandbox/github/PHD/data/pdbbind/refined-set"
TARGET_DIR = "/Users/tevfik/Sandbox/github/PHD/data/pdbbind/refined-set-only-used-in-codon-tests"

def read_yaml_proteins(yaml_path):
    """YAML dosyasından train, val, test protein listelerini oku."""
    with open(yaml_path, 'r') as f:
        content = f.read()
    
    # Basit regex ile protein listelerini çıkar
    train_proteins = re.findall(r'train:\s*\n((?:\s+-\s+\w+\n)+)', content, re.MULTILINE)
    val_proteins = re.findall(r'validation:\s*\n((?:\s+-\s+\w+\n)+)', content, re.MULTILINE)
    test_proteins = re.findall(r'test:\s*\n((?:\s+-\s+\w+\n)+)', content, re.MULTILINE)
    
    def extract_names(text):
        """Liste string'inden protein isimlerini çıkar."""
        return re.findall(r'-\s+(\w+)', text)
    
    train_list = extract_names(train_proteins[0]) if train_proteins else []
    val_list = extract_names(val_proteins[0]) if val_proteins else []
    test_list = extract_names(test_proteins[0]) if test_proteins else []
    
    all_proteins = set(train_list + val_list + test_list)
    
    print(f"Train: {len(train_list)} proteins")
    print(f"Validation: {len(val_list)} proteins")
    print(f"Test: {len(test_list)} proteins")
    print(f"Total unique: {len(all_proteins)} proteins")
    
    return all_proteins

def copy_protein_folders(proteins, source_dir, target_dir):
    """Protein klasörlerini kaynak dizinden hedef dizine kopyala."""
    
    # Hedef dizini oluştur
    os.makedirs(target_dir, exist_ok=True)
    
    copied = 0
    skipped = 0
    missing = 0
    
    for protein in sorted(proteins):
        source_path = os.path.join(source_dir, protein)
        target_path = os.path.join(target_dir, protein)
        
        # Kaynak klasör var mı kontrol et
        if not os.path.exists(source_path):
            print(f"⚠️  MISSING: {protein} (kaynak dizinde bulunamadı)")
            missing += 1
            continue
        
        # Hedef klasör zaten var mı?
        if os.path.exists(target_path):
            print(f"⏭️  SKIP: {protein} (zaten var)")
            skipped += 1
            continue
        
        # Klasörü kopyala
        try:
            shutil.copytree(source_path, target_path)
            print(f"✅ COPIED: {protein}")
            copied += 1
        except Exception as e:
            print(f"❌ ERROR copying {protein}: {e}")
    
    print(f"\n{'='*60}")
    print(f"Özet:")
    print(f"  Kopyalanan: {copied}")
    print(f"  Zaten var (atlandı): {skipped}")
    print(f"  Kaynak dizinde bulunamayan: {missing}")
    print(f"{'='*60}")

if __name__ == "__main__":
    print(f"YAML dosyası: {YAML_FILE}")
    print(f"Kaynak dizin: {SOURCE_DIR}")
    print(f"Hedef dizin: {TARGET_DIR}")
    print()
    
    # YAML'den protein listelerini oku
    proteins = read_yaml_proteins(YAML_FILE)
    print()
    
    # Kopyalama işlemini başlat
    copy_protein_folders(proteins, SOURCE_DIR, TARGET_DIR)

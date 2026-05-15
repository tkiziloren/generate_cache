# Gridfix Cache Pipeline

Bu pipeline mevcut `generate_cache_pdbbind_resolution_fix.py` dosyasini degistirmez. Yeni ciktilar ayri klasore yazilir:

```text
/Users/tevfik/Sandbox/github/PHD/data/pdbbind/refined-set_minimal_cache_gridfix_v1/box72/
```

## Neden Ayri Pipeline?

APBS icin iyi calisan kaynak grid `161 x 161 x 161` grid noktasi olarak korunur. Bu grid `160 A` fiziksel span tasir: `center - 80 A` ile `center + 80 A` arasi, 1 A grid araligi.

Eski `resolution_fix` yolunda `resolution = 161.0 / box_size` kullaniliyordu. Bu, 161 grid noktasi ile 161 A fiziksel uzunlugu karistiriyor. Yeni gridfix yolu hedef grid icin su formulu kullanir:

```text
resolution = 160.0 / (box_size - 1)
```

`box_size=72` icin:

```text
resolution = 2.25352113 A/voxel
span = 71 * 2.25352113 = 160 A
```

Boylece 72 grid de 161 grid ile ayni fiziksel alani kaplar: `center +/- 80 A`.

## Neler Duzeldi?

- APBS kaynak grid 161 nokta / 160 A olarak uretilir.
- APBS `.dx` header'i `counts`, `origin`, `delta` ile parse edilir.
- 72 grid'e resampling fiziksel koordinat uzerinden yapilir.
- Atomic feature'lar `tfbio.make_grid` ile otomatik box hesabi yaptirilmadan, ayni `grid_origin` ve `resolution` ile voxelize edilir.
- Tum feature ve label tensorleri kaydetmeden once ayni shape icin assert edilir.
- Kaynak PDBBind dosyalari degistirilmez; protein/ligand gecici klasore kopyalanir.

## Tek Protein Uretme

```bash
cd /Users/tevfik/Sandbox/Tevfik/Projects/phd_examples/generate_cache
PATH="$PWD/.conda/bin:$PATH" .conda/bin/python generate_cache_pdbbind_gridfix.py \
  --cases 1a4w \
  --box-size 72 \
  --nproc 1 \
  --pdb2pqr-bin "$PWD/.conda/bin/pdb2pqr30" \
  --fail-on-error
```

## Kucuk Local Smoke Set

```bash
cd /Users/tevfik/Sandbox/Tevfik/Projects/phd_examples/generate_cache
PATH="$PWD/.conda/bin:$PATH" .conda/bin/python generate_cache_pdbbind_gridfix.py \
  --case-list pdbbind_gridfix160_full_fits.txt \
  --limit 5 \
  --box-size 72 \
  --nproc 1 \
  --pdb2pqr-bin "$PWD/.conda/bin/pdb2pqr30" \
  --fail-on-error
```

## Tum Local Fits Listesini Uretme

Eski `pdbbind_60_fits.txt` ve `pdbbind_72_fits.txt` listeleri Angstrom cinsinden protein kutu boyuna gore uretilmisti. Gridfix `box72`, 72 A degil 160 A span kullandigi icin bu listeler artik fazla konservatif kalir.

Tum refined-set uzerinde yapilan gridfix coverage kontrolu:

```text
valid protein cases: 5316
full protein fits center +/-80 A: 5300 / 5316
pocket fits center +/-80 A:       5316 / 5316
```

Mevcut eski listeler icin sonuc:

```text
pdbbind_60_fits.txt:  full 1640/1640, pocket 1640/1640
pdbbind_72_fits.txt:  full 2611/2611, pocket 2611/2611
pdbbind_128_fits.txt: full 3739/3740, pocket 3740/3740
```

Yani `box72` gridfix icin sadece 72 A icine giren kucuk proteinlerle sinirlamak gerekmiyor. Tam protein gridde kalsin istiyorsak 16 buyuk/asimetrik case'i disarida birakmak yeterli; pocket/label kapsami acisindan tum valid case'ler grid icinde kaliyor.

Coverage'i tekrar uretmek icin:

```bash
cd /Users/tevfik/Sandbox/Tevfik/Projects/phd_examples/generate_cache
.conda/bin/python audit_gridfix_coverage.py \
  --box-size 72 \
  --nproc 8 \
  --write-prefix /Users/tevfik/Sandbox/github/PHD/data/pdbbind/refined-set_minimal_cache_gridfix_v1/coverage/gridfix160
```

Bu komut su dosyalari yazar:

```text
gridfix160_coverage.csv
gridfix160_full_fits.txt
gridfix160_full_unfits.txt
gridfix160_pocket_fits.txt
gridfix160_pocket_unfits.txt
```

```bash
cd /Users/tevfik/Sandbox/Tevfik/Projects/phd_examples/generate_cache
PATH="$PWD/.conda/bin:$PATH" .conda/bin/python generate_cache_pdbbind_gridfix.py \
  --case-list pdbbind_gridfix160_full_fits.txt \
  --box-size 72 \
  --nproc 4 \
  --pdb2pqr-bin "$PWD/.conda/bin/pdb2pqr30"
```

## Cache Dogrulama

```bash
cd /Users/tevfik/Sandbox/Tevfik/Projects/phd_examples/generate_cache
.conda/bin/python validate_cache_gridfix.py \
  --h5-dir /Users/tevfik/Sandbox/github/PHD/data/pdbbind/refined-set_minimal_cache_gridfix_v1/box72 \
  --limit 5 \
  --fail-on-error
```

Validator ayni klasore su CSV ozetini yazar:

```text
cache_validation_summary.csv
```

## Ilk Smoke Sonucu

`1a4w` icin local smoke uretimi yapildi:

```text
/Users/tevfik/Sandbox/github/PHD/data/pdbbind/refined-set_minimal_cache_gridfix_v1/box72/1a4w.h5
```

Validator sonucu:

```text
[OK] 1a4w
Checked files: 1
Failed files: 0
```

Onemli shape kontrolu:

```text
old resolution_fix atomic_C: (33, 33, 33)
new gridfix atomic_C:      (72, 72, 72)
new gridfix labels:        (72, 72, 72)
```

APBS degeri eski 72 cache ile pratik olarak ayni kaldi:

```text
old-vs-new APBS corr = 0.999999
old-vs-new APBS MAE  = 0.000068
```

Kontrol ettigi ana kurallar:

- Tum `features/*` tensorleri `(box_size, box_size, box_size)` olmali.
- Tum `label/*` tensorleri ayni shape'te olmali.
- Label tensorleri bos olmamali.
- Continuous feature'larda NaN/Inf olmamali.
- `resolution * (box_size - 1)` ile `physical_span_angstrom` uyumlu olmali.

## Egitimde Kullanma

`3dunet_configurable` config dosyasinda `h5_directory` su klasore alinmali:

```yaml
h5_directory: "/Users/tevfik/Sandbox/github/PHD/data/pdbbind/refined-set_minimal_cache_gridfix_v1/box72"
```

Bu cache uretildikten sonra `electrostatic_grid + shape`, `electrostatic + atomic` ve `full_context` ablation'lari ayni grid uzerinde denenebilir.

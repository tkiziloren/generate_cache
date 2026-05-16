# Deep-APBS H5 Öznitelik Şeması

Bu doküman, Deep-APBS cache dosyalarında kullanılan HDF5 şemasını açıklar. Amaç, üretilen veri setinin başka bir araştırmacı tarafından da anlaşılabilir, tekrar kullanılabilir ve model girdisi açısından güvenli olmasıdır.

## H5 Dosya Yapısı

| Grup | Anlamı | Model girdisi olarak kullanılmalı mı? |
|---|---|---|
| `features/` | Kör bağlanma bölgesi tahmini için kullanılabilecek protein tabanlı öznitelikler. | Evet, ablation/deney tasarımına göre seçilir. |
| `label/` | Eğitim ve değerlendirme için kullanılan bağlanma bölgesi etiketleri. | Hayır, sadece hedef/ground truth. |
| `auxiliary/` | Görselleştirme, hata analizi veya veri kontrolü için tutulan ligand kaynaklı yardımcı kanallar. | Hayır. Kör tahminde ligand bilinmediği için leakage oluşturur. |
| root attributes | Grid, APBS ve şema metadata bilgileri. | Model girdisi değildir. |

`features/electrostatic_grid` eski isimlendirmedir ve public şemada tutulmaz. Bu kanal, açık isimli `electrostatic_grid_v1_ligand_proximal_chains_7A_raw` üretildikten sonra silinir.

## Grid ve Fiziksel Ölçek

| Ayar | Anlamı | Yaklaşık çözünürlük |
|---|---|---|
| `box36_span70` | Kalasanty/PUResNet benzeri küçük kutu, 70 Å fiziksel alan. | `70 / 35 = 2.00 Å/voxel` |
| `box72_span120` | Daha geniş bağlam ve 36'ya göre daha iyi çözünürlük dengesi. | `120 / 71 = 1.69 Å/voxel` |
| `box161_span160` | APBS alanına en yakın tam çözünürlüklü kutu. | `160 / 160 = 1.00 Å/voxel` |

APBS hesaplaması varsayılan olarak 161 grid noktası ve 160 Å fiziksel span ile yapılır. Daha küçük kutularda APBS alanı hedef grid üzerine örneklenir.

## APBS Kaynak Tanımları

| Prefix | Açıklama | Tez açısından rolü |
|---|---|---|
| `electrostatic_grid_v1_ligand_proximal_chains_7A_*` | Eski pipeline ile uyumlu APBS. Ligand çevresine yakın protein chain/residue seçimi üzerinden hesaplanır. | Legacy sonuçlarla karşılaştırma ve önceki deneyleri tekrar üretme. |
| `electrostatic_grid_v2_full_protein_*` | Ligand çevresi üzerinden chain kırpmadan, tüm protein yapısı üzerinden hesaplanan APBS. | Kör tahmin senaryosuna daha temiz ve daha savunulabilir APBS temsili. |

## APBS Temsil Türleri

Aynı APBS potansiyel alanı farklı normalize edilmiş kanallar olarak saklanır. Bunun nedeni elektrostatik potansiyelin uç değerler içerebilmesi ve farklı modellerin farklı ölçeklere duyarlı olmasıdır.

| Suffix | Açıklama | Değer aralığı | Kullanım yorumu |
|---|---|---|---|
| `raw` | APBS potansiyelinin hedef grid üzerindeki ham hali. | Sınırsız | Analiz ve alternatif normalizasyonlar için saklanır. |
| `clip5_minmax` | Değerler `[-5, 5]` aralığına kırpılır, sonra `[0, 1]` aralığına ölçeklenir. | `[0, 1]` | Zayıf ve orta elektrostatik bölgeleri öne çıkarır. |
| `clip10_minmax` | Değerler `[-10, 10]` aralığına kırpılır, sonra `[0, 1]` aralığına ölçeklenir. | `[0, 1]` | Orta ölçekli elektrostatik sinyali korur. |
| `clip20_minmax` | Değerler `[-20, 20]` aralığına kırpılır, sonra `[0, 1]` aralığına ölçeklenir. | `[0, 1]` | Daha güçlü potansiyel farklarını da taşır. |
| `clip20_signed` | Değerler `[-20, 20]` aralığına kırpılır, sonra `[-1, 1]` aralığına ölçeklenir. | `[-1, 1]` | Pozitif/negatif işaret bilgisini korur. |
| `full_signed150` | Ham değer `150` değerine bölünür, kırpma yapılmaz. | Genellikle yaklaşık `[-1, 1]`, ama teorik olarak sınırsız | Uç değerleri tamamen silmeden imzalı alanı verir. |
| `clip150_signed` | Değerler `[-150, 150]` aralığına kırpılır, sonra `[-1, 1]` aralığına ölçeklenir. | `[-1, 1]` | Full signed fikrinin güvenli kırpılmış versiyonu. |
| `positive_clip20` | Pozitif APBS bileşeni `[0, 20]` aralığına kırpılır ve `[0, 1]` yapılır. | `[0, 1]` | Pozitif elektrostatik bölgeleri ayrı kanal olarak verir. |
| `negative_clip20` | Negatif APBS bileşeni mutlak büyüklük olarak `[0, 20]` aralığına kırpılır ve `[0, 1]` yapılır. | `[0, 1]` | Negatif elektrostatik bölgeleri ayrı kanal olarak verir. |
| `gradient_magnitude_robust` | APBS alanının uzaysal gradient büyüklüğü hesaplanır ve robust percentile ile `[0, 1]` aralığına normalize edilir. | `[0, 1]` | Elektrostatik potansiyelin hızlı değiştiği bölgeleri gösterir. |
| `clip20_signed_surface_weighted` | `clip20_signed * protein_proximity_exp3`. | `[-1, 1]` civarı | Elektrostatik sinyali protein/yüzey yakınında vurgular. |
| `full_signed150_surface_weighted` | `full_signed150 * protein_proximity_exp3`. | Genellikle `[-1, 1]` civarı | Kırpmasız signed APBS bilgisini protein/yüzey yakınında vurgular. |

## Public Model Feature Listesi

### Atomik Kanallar

| Feature | Açıklama |
|---|---|
| `atomic_B` | Voxelize edilmiş bor atomu varlığı. |
| `atomic_C` | Voxelize edilmiş karbon atomu varlığı. |
| `atomic_N` | Voxelize edilmiş azot atomu varlığı. |
| `atomic_O` | Voxelize edilmiş oksijen atomu varlığı. |
| `atomic_P` | Voxelize edilmiş fosfor atomu varlığı. |
| `atomic_S` | Voxelize edilmiş kükürt atomu varlığı. |
| `atomic_Se` | Voxelize edilmiş selenyum atomu varlığı. |
| `atomic_acceptor` | Hidrojen bağı alıcısı atom/property kanalı. |
| `atomic_aromatic` | Aromatik yapı veya aromatik atom katılımı. |
| `atomic_donor` | Hidrojen bağı vericisi atom/property kanalı. |
| `atomic_halogen` | Halojen atom/property kanalı. |
| `atomic_heavydegree` | Atomun ağır atomlarla bağ sayısını temsil eden kanal. |
| `atomic_heterodegree` | Atomun hetero atomlarla bağ sayısını temsil eden kanal. |
| `atomic_hyb` | Atom hibridizasyon bilgisinin voxelize edilmiş hali. |
| `atomic_hydrophobic` | Atom seviyesinde hidrofobiklik özelliği. |
| `atomic_metal` | Metal atom/property kanalı. |
| `atomic_molcode` | Atom featurizer tarafından üretilen molekül kodu kanalı. |
| `atomic_partialcharge` | Atom kısmi yük bilgisinin voxelize edilmiş hali. |
| `atomic_ring` | Atomun halka yapısına dahil olup olmadığını temsil eder. |

### Geometrik ve Fizikokimyasal Kanallar

| Feature | Açıklama |
|---|---|
| `shape` | Protein atom doluluk maskesi. Protein hacmini ve genel şekli temsil eder. |
| `hydrophobicity` | Amino asit hidrofobiklik değerlerinin grid üzerine rasterize edilmiş hali. |
| `dist_to_surface` | Her grid noktasının protein atomlarına/yüzey proxy'sine uzaklığı. Cep derinliği ve yüzey yakınlığı hakkında sinyal taşır. |
| `protein_proximity_exp3` | `exp(-dist_to_surface / 3A)` ile hesaplanan protein yakınlığı kanalı. Yüzeye/protein atomlarına yakın bölgeleri yumuşak biçimde vurgular. |
| `protein_near_shell_0_3A` | `dist_to_surface` değeri 0-3 Å arasında olan yakın protein çevresi shell maskesi. |
| `protein_near_shell_3_6A` | `dist_to_surface` değeri 3-6 Å arasında olan daha dış protein çevresi shell maskesi. |
| `hydrophobicity_surface_weighted` | `hydrophobicity * protein_proximity_exp3`. Hidrofobik sinyali protein/yüzey yakınında vurgular. |

### APBS v1 Kanalları

Aşağıdaki tüm kanallar ligand-proximal legacy APBS kaynağından türetilir:

- `electrostatic_grid_v1_ligand_proximal_chains_7A_raw`
- `electrostatic_grid_v1_ligand_proximal_chains_7A_clip5_minmax`
- `electrostatic_grid_v1_ligand_proximal_chains_7A_clip10_minmax`
- `electrostatic_grid_v1_ligand_proximal_chains_7A_clip20_minmax`
- `electrostatic_grid_v1_ligand_proximal_chains_7A_clip20_signed`
- `electrostatic_grid_v1_ligand_proximal_chains_7A_full_signed150`
- `electrostatic_grid_v1_ligand_proximal_chains_7A_clip150_signed`
- `electrostatic_grid_v1_ligand_proximal_chains_7A_positive_clip20`
- `electrostatic_grid_v1_ligand_proximal_chains_7A_negative_clip20`
- `electrostatic_grid_v1_ligand_proximal_chains_7A_gradient_magnitude_robust`
- `electrostatic_grid_v1_ligand_proximal_chains_7A_clip20_signed_surface_weighted`
- `electrostatic_grid_v1_ligand_proximal_chains_7A_full_signed150_surface_weighted`

### APBS v2 Kanalları

Aşağıdaki tüm kanallar tüm protein üzerinden hesaplanan APBS kaynağından türetilir:

- `electrostatic_grid_v2_full_protein_raw`
- `electrostatic_grid_v2_full_protein_clip5_minmax`
- `electrostatic_grid_v2_full_protein_clip10_minmax`
- `electrostatic_grid_v2_full_protein_clip20_minmax`
- `electrostatic_grid_v2_full_protein_clip20_signed`
- `electrostatic_grid_v2_full_protein_full_signed150`
- `electrostatic_grid_v2_full_protein_clip150_signed`
- `electrostatic_grid_v2_full_protein_positive_clip20`
- `electrostatic_grid_v2_full_protein_negative_clip20`
- `electrostatic_grid_v2_full_protein_gradient_magnitude_robust`
- `electrostatic_grid_v2_full_protein_clip20_signed_surface_weighted`
- `electrostatic_grid_v2_full_protein_full_signed150_surface_weighted`

## Etiketler

| Label | Açıklama | Kullanım |
|---|---|---|
| `binding_site_calculated` | Protein-ligand yakınlığına göre pipeline tarafından hesaplanan bağlanma bölgesi etiketi. | Eğitim veya değerlendirme. |
| `binding_site_in_dataset` | Kaynak veri setinin sağladığı bağlanma bölgesi etiketi. | Eğitim veya değerlendirme. |
| `binding_site_fpocket_selected` | BU48/COACH gibi dış test setlerinde ligand konumuna göre seçilen fpocket pocket etiketi. | Dış veri seti değerlendirme. |

## Yardımcı Kanallar

| Auxiliary | Açıklama | Model girdisi olmalı mı? |
|---|---|---|
| `ligand` | Voxelize ligand maskesi. | Hayır. Sadece görselleştirme/debug. |
| `dist_to_ligand` | Grid noktalarının ligand atomlarına uzaklığı. | Hayır. Ligand konumunu doğrudan sızdırır. |

## Önerilen Ablation Grupları

| Grup | Feature set |
|---|---|
| Shape only | `shape` |
| APBS only | Tek bir APBS temsil kanalı, tercihen güçlü aday olarak `electrostatic_grid_v2_full_protein_full_signed150` veya önceki sonuçlara göre en iyi v1/v2 temsil. |
| Shape + APBS | `shape` + seçilen APBS kanalı. |
| Selected chemistry | `atomic_C`, `atomic_N`, `atomic_O`, `atomic_S`, `atomic_acceptor`, `atomic_donor`, `atomic_aromatic`, `atomic_hydrophobic`, `atomic_partialcharge`, `hydrophobicity` |
| APBS + selected chemistry | Seçilen APBS kanalı + selected chemistry. |
| APBS + shape + selected chemistry | Seçilen APBS kanalı + `shape` + selected chemistry. |
| Surface/hydro add-on | Yukarıdaki güçlü gruplara `dist_to_surface` ve `hydrophobicity` eklenmiş varyantlar. |

## Güvenlik Notu

Blind prediction senaryosunda kullanıcı sadece protein yapısını verecektir. Bu yüzden `label/`, `auxiliary/ligand` ve `auxiliary/dist_to_ligand` model girdisi olamaz. Bu kanallar sadece doğrulama, görselleştirme ve hata analizi için tutulur.

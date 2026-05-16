# Deep-APBS Öznitelik Literatür Taraması

Bu doküman, protein-ligand bağlanma bölgesi tahmini literatüründe kullanılan öznitelik ailelerini, mevcut Deep-APBS şemamızla ilişkilerini ve bundan sonra eklenebilecek öznitelik adaylarını özetler. Odak noktası, kör protein-only binding-site prediction senaryosudur: model inference sırasında ligand konumu bilinmez.

## Kısa Sonuç

Mevcut Deep-APBS şeması literatürdeki ana temsil ailelerinin büyük bölümünü kapsıyor: atomik kanallar, protein şekli, hidrofobiklik, yüzey uzaklığı ve APBS elektrostatik potansiyeli. En özgün katkı APBS'nin sistematik biçimde farklı normalizasyonlarla, özellikle APBS-only ve APBS+kimyasal kombinasyonları içinde ölçülmesidir.

Literatüre göre en mantıklı yeni öznitelik adayları şunlardır:

1. **Yüzey/cep geometrisi:** SAS noktaları, protrusion, curvature, concavity, alpha-sphere/fpocket tabanlı yoğunluk.
2. **Enerji tabanlı gridler:** Lennard-Jones/van der Waals, hidrojen bağı potansiyeli, Coulomb/partial-charge alanı.
3. **Evrimsel ve embedding tabanlı bilgi:** conservation skorları, PSSM/HMM profilleri, ESM-2 veya ESM-IF benzeri residue embeddingleri.
4. **Molecular interaction field benzeri probe kanalları:** hidrofobik, donor, acceptor, pozitif ve negatif probe enerjileri.
5. **Yapısal bağlam:** solvent accessibility, residue depth, secondary structure, B-factor/flexibility.

10 günlük tez teslim takvimi açısından en doğru sıra: önce mevcut APBS v1/v2 şemasını sağlamlaştırmak, sonra `dist_to_surface` ve mevcut atomik/kimyasal kanallarla ablation yapmak, daha sonra yüzey/cep geometrisi ve basit enerji gridlerini eklemek. Evrimsel embedding ve MaSIF/dMaSIF benzeri surface deep learning yaklaşımları güçlü ama tez teslimi öncesi riskli ve ayrı çalışma olacak kadar büyük işlerdir.

## Mevcut Deep-APBS Öznitelik Aileleri

| Aile | Bizdeki karşılığı | Literatürdeki yeri | Değerlendirme |
|---|---|---|---|
| Atom türü ve atomik özellikler | `atomic_*` kanalları | PUResNet, Kalasanty, DeepSite, DeepPocket gibi 3B grid yaklaşımlarında yaygın. | Güçlü ve standart baseline. |
| Protein şekli | `shape` | LIGSITE, SURFNET, fpocket, Kalasanty/PUResNet tarzı grid yöntemlerinde temel sinyal. | En temel geometrik sinyal. |
| Hidrofobiklik | `hydrophobicity`, `atomic_hydrophobic` | Cep karakterizasyonunda ve surface yöntemlerde sık kullanılır. | Cep oluşumu için biyofiziksel olarak anlamlı. |
| Yüzey uzaklığı | `dist_to_surface` | Surface/depth/concavity fikrine yakın. | Leakage yok; ablation'a eklenmeli. |
| Elektrostatik potansiyel | APBS v1/v2 kanalları | Elektrostatik potansiyel literatürde bilinen bir biyofiziksel sinyal, ancak Kalasanty/PUResNet ekseninde sistematik APBS grid ablation katkısı daha özgün. | Tezin ana katkısı. |
| Ligand maskesi ve ligand uzaklığı | `auxiliary/ligand`, `auxiliary/dist_to_ligand` | Eğitim target/debug için faydalı olabilir, ama blind prediction girdisi olamaz. | Public şemada auxiliary; training feature değildir. |

## Literatürdeki Ana Çizgiler

### 1. Kalasanty: 3B U-Net ile segmentation yaklaşımı

Kalasanty, bağlanma bölgesi tahminini 3B image segmentation problemi olarak kurar. Protein grid'e çevrilir ve model her grid noktasının pocket olup olmadığını tahmin eder. U-Net mimarisi encoder-decoder ve skip connection yapısıyla lokal detayları korumaya çalışır. Makalede DCC ve DVO'nun birbirini tamamladığı vurgulanır: DCC cebin doğru lokasyonda olup olmadığını, DVO ise şekil örtüşmesini ölçer.

Deep-APBS açısından anlamı: Bizim `shape`, atomik kanallar ve voxel mask değerlendirmemiz Kalasanty çizgisine yakın. Farkımız, APBS elektrostatik potansiyel alanını sistematik bir feature ailesi olarak eklememiz ve APBS-only performansı ayrıca ölçmemizdir.

Kaynak: Stepniewska-Dziubinska ve ark., 2020, Scientific Reports, "Improving detection of protein-ligand binding sites with 3D segmentation"  
https://www.nature.com/articles/s41598-020-61860-z

### 2. PUResNet: 36 x 36 x 36 x 18 atomik grid ve residual U-Net benzeri yapı

PUResNet, proteinleri `36 x 36 x 36 x 18` grid olarak temsil eder. Makalede dokuz atomik özellik ve bunların one-hot/kanal temsilleri kullanılır: hybridization, heavy atoms, heteroatoms, hydrophobic, aromatic, partial charge, acceptor, donor ve ring. Mimari U-Net ve ResNet fikrini birleştirir; encoder-decoder yapı içinde residual/identity bağlantılar kullanır.

Deep-APBS açısından anlamı: Bizim atomik kanallarımız PUResNet'in kullandığı atomik feature ailesine çok yakın. Bu nedenle APBS katkısını ölçerken PUResNet benzeri atomik setleri güçlü bir referans olarak kullanabiliriz.

Kaynak: Kandel ve ark., 2021, Journal of Cheminformatics, "PUResNet: prediction of protein-ligand binding sites using deep residual neural network"  
https://jcheminf.biomedcentral.com/articles/10.1186/s13321-021-00547-7

### 3. PUResNetV2.0: sparse representation ve yeni dış veri setleri

PUResNetV2.0, PUResNet çizgisini genişletip sparse representation ile daha verimli ve güçlü hale getirmeyi hedefler. Holo801 ve Apoholo45 gibi farklı durumları içeren veri setlerinde performans raporlar. Bu çalışma, sadece scPDB içi başarı değil, farklı dataset ve apo/holo bağlamında genelleme konusunu öne çıkarır.

Deep-APBS açısından anlamı: Bizim de BU48, COACH ve PDBBind gibi dış değerlendirme setlerini hazırlamamız doğru yöndür. APBS'nin gerçek katkısı yalnızca local fold skorunda değil, dış veri setlerine genellenip genellenmediğinde daha güçlü savunulur.

Kaynak: Kandel ve ark., 2024, Journal of Cheminformatics, "PUResNetV2.0: a deep learning model leveraging sparse representation for improved ligand binding site prediction"  
https://jcheminf.biomedcentral.com/articles/10.1186/s13321-024-00865-6

### 4. DeepSite: 3B CNN ve atom-type grid yaklaşımı

DeepSite, protein çevresini 3B grid alt bölgelerine ayırıp 3B CNN ile bağlanma olasılığı üretir. Bu yaklaşım U-Net segmentation'dan daha çok local patch classification mantığına yakındır. Atom/property type kanalları temel girdidir.

Deep-APBS açısından anlamı: DeepSite, atomik grid + 3B CNN yaklaşımının erken ve güçlü bir örneğidir. Bizim segmentation tabanlı modelimiz bu çizgiyi U-Net/ResNet/UNet++ ile genişletir.

Kaynak: Jimenez ve ark., 2017, Bioinformatics, "DeepSite: protein-binding site predictor using 3D-convolutional neural networks"  
https://academic.oup.com/bioinformatics/article/33/19/3036/3859178

### 5. P2Rank/PrankWeb: solvent-accessible surface noktaları ve ligandability score

P2Rank grid yerine protein solvent-accessible surface üzerinde örneklenen SAS noktalarını sınıflandırır. Her SAS noktası local çevresinden hesaplanan fizikokimyasal, geometrik ve istatistiksel feature vektörü ile temsil edilir. Makalede 35 sayısal feature kullanıldığı ve en önemli feature'lardan birinin protrusion olduğu belirtilir.

Deep-APBS açısından anlamı: Bizim grid tabanlı yaklaşımımıza ek olarak yüzey noktası veya yüzeyden türetilmiş channel eklemek mantıklı. `dist_to_surface`, bu yönde ilk basit adımdır. Daha güçlü adaylar: SAS point density, protrusion, curvature ve local surface patch descriptors.

Kaynak: Krivak ve Hoksza, 2018, Journal of Cheminformatics, "P2Rank: machine learning based tool for rapid and accurate prediction of ligand binding sites from protein structure"  
https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0285-8

### 6. fpocket ve alpha-sphere tabanlı cep geometrisi

fpocket, Voronoi tessellation ve alpha-sphere fikrine dayanır. Alpha sphere, yüzeydeki boşlukları ve kavisli bölgeleri temsil eder. Orta büyüklükte alpha sphere kümeleri pocket adayları olarak yorumlanır; daha sonra pocket özellikleri skorlanır.

Deep-APBS açısından anlamı: fpocket tek başına model değildir, ama feature üretici olarak çok değerli olabilir. Grid'e alpha-sphere yoğunluğu, fpocket pocket score, pocket rank veya pocket mask olarak eklenebilir. Ancak dikkat: fpocket seçimi label üretmek için kullanılıyorsa aynı bilgiyi input olarak vermek değerlendirmeyi karıştırabilir. Bu nedenle fpocket feature kullanımı ayrı ve net bir deney olmalıdır.

Kaynak: Le Guilloux ve ark., 2009, BMC Bioinformatics, "Fpocket: An open source platform for ligand pocket detection"  
https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-168

### 7. DeepPocket: fpocket adaylarını 3B CNN ile yeniden skorlama ve segmentation

DeepPocket, fpocket'ın ürettiği cep adaylarını 3B CNN ile yeniden skorlar ve ayrıca segmentation modeliyle pocket şeklini iyileştirir. Bu nedenle saf end-to-end grid segmentation değil, aday cep üretimi + deep learning refinement pipeline'ıdır.

Deep-APBS açısından anlamı: Eğer top-k pocket metriklerini iyileştirmek istiyorsak, postprocess ve candidate-pocket refinement kritik olabilir. Bizim modelin ham voxel maskesini connected component, closing, clear border, min volume ve top-k seçimiyle Kalasanty/DeepPocket çizgisine yaklaştırmamız önemlidir.

Kaynak: Aggarwal ve ark., 2021/2022, Journal of Chemical Information and Modeling, "DeepPocket: Ligand Binding Site Detection and Segmentation using 3D Convolutional Neural Networks"  
https://pubs.acs.org/doi/10.1021/acs.jcim.1c00799

### 8. ConCavity: yapı + evrimsel korunum

ConCavity, protein 3B yapısını evrimsel sequence conservation ile birleştirir. Makalede yapı tabanlı yöntemlerin genellikle conservation-only yöntemlerden güçlü olduğu, ama ikisinin birleşiminin daha iyi sonuç verdiği gösterilir.

Deep-APBS açısından anlamı: Eğer ek feature arayacaksak conservation çok güçlü bir biyolojik sinyal olabilir. Ancak MSA/profil üretimi computational olarak pahalıdır ve veri pipeline'ını büyütür. Tez teslimi öncesi değil, makale genişletmesi için daha uygun olabilir.

Kaynak: Capra ve ark., 2009, PLOS Computational Biology, "Predicting Protein Ligand Binding Sites by Combining Evolutionary Sequence Conservation and 3D Structure"  
https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000585

### 9. DeepSurf ve surface tabanlı deep learning

DeepSurf, protein yüzeyindeki local voxel gridleri ve fizikokimyasal özellikleri kullanarak surface-based ligandability tahmini yapar. Bu aile, tüm 3B kutuyu doldurmak yerine yüzeye yakın bölgeleri modellemeye çalışır.

Deep-APBS açısından anlamı: Bağlanma cepleri çoğunlukla yüzey ve yüzeye yakın kavitelerde oluştuğu için surface-centered representation verimlidir. Bizim `dist_to_surface` channel'ı bu fikrin basit grid karşılığıdır. Daha ileri aşamada surface point sampling veya local surface patch özellikleri eklenebilir.

Kaynak: Mylonas ve Axenopoulos, 2021, Bioinformatics, "DeepSurf: a surface-based deep learning approach for the prediction of ligand binding sites on proteins"  
https://academic.oup.com/bioinformatics/article/37/12/1681/6104838

### 10. MaSIF/dMaSIF: moleküler yüzey üzerinde geometric deep learning

MaSIF, protein moleküler yüzeyini kimyasal ve geometrik surface fingerprint olarak temsil eder. Ligand-binding pocket, protein-protein interaction ve surface matching gibi görevlerde geometric deep learning kullanır.

Deep-APBS açısından anlamı: Bu yaklaşım güçlü ama grid tabanlı tez pipeline'ından farklı bir model ailesidir. Tez için ana çizgiyi dağıtmadan "gelecek çalışma" olarak konumlandırılmalı. Ancak yüzey curvature, shape index ve chemical surface descriptors bizim grid feature'larımıza dönüştürülebilir.

Kaynak: Gainza ve ark., 2020, Nature Methods, "Deciphering interaction fingerprints from protein molecular surfaces using geometric deep learning"  
https://www.nature.com/articles/s41592-019-0666-6

### 11. Enerji tabanlı grid descriptor çalışmaları

Bazı çalışmalar grid tabanlı protein descriptor'larını şekil, van der Waals/Lennard-Jones potansiyeli, hidrojen bağı potansiyeli ve Coulomb/partial charge gibi enerji kanallarıyla kurar. Bu fikir APBS'ye çok yakındır: protein etrafındaki fiziksel alanlar model input'u olur.

Deep-APBS açısından anlamı: APBS ana katkı olarak korunurken, yanına basit `vdw_potential_grid`, `hbond_donor_potential_grid`, `hbond_acceptor_potential_grid` gibi probe/energy kanalları eklemek makale için güçlü bir ikinci aşama olabilir.

Kaynak: Jiang ve ark., 2019, BMC Bioinformatics, "A novel protein descriptor for the prediction of drug binding sites"  
https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-3058-0

### 12. APBS ve elektrostatik potansiyel

APBS, Poisson-Boltzmann temelli elektrostatik potansiyel alanı üretir. Resmi dokümantasyon, APBS/PDB2PQR akışının PyMOL ile elektrostatik isosurface ve surface potential görselleştirmek için kullanılmasını tarif eder. Elektrostatik potansiyel özellikle uzak etkili, işaretli ve fiziksel anlam taşıyan bir protein alanıdır.

Deep-APBS açısından anlamı: APBS'nin sadece görselleştirme değil, model input'u olarak sistematik test edilmesi tezin ayrıştırıcı yönüdür. Özellikle pozitif/negatif ayrıştırma, clipped signed, full signed ve full-protein APBS v2 varyantları savunulabilir feature engineering katkılarıdır.

Kaynak: APBS documentation, "Visualization with PyMOL"  
https://apbs.readthedocs.io/en/stable/using/examples/visualization-pymol.html

### 13. Büyük karşılaştırmalı benchmark çalışmaları

2024 tarihli karşılaştırmalı çalışma, modern binding-site prediction araçlarının feature ailelerini ve postprocess stratejilerini yan yana verir. PUResNet'in atom+one-hot encoding ile 18 feature kullandığı, P2Rank'in 35 atom/residue feature ve SAS noktaları kullandığı, DeepPocket'ın fpocket + atom feature + 3B CNN çizgisinde olduğu belirtilir.

Deep-APBS açısından anlamı: Tezde sadece tek bir metrikle değil, DCC/DCA/DVO/F1/top-k gibi literature-style metriklerle raporlamak doğru. Ayrıca postprocess ve top-k protokolünün sonucu doğrudan etkilediği açıkça anlatılmalı.

Kaynak: "Comparative evaluation of methods for the prediction of protein-ligand binding sites", Journal of Cheminformatics, 2024  
https://jcheminf.biomedcentral.com/articles/10.1186/s13321-024-00923-z

## Eklenebilecek Öznitelik Adayları

### A. Hemen uygulanabilir adaylar

| Aday feature | Nasıl üretilebilir? | Neden faydalı olabilir? | Risk |
|---|---|---|---|
| `dist_to_surface` ablation | Zaten şemada var. | Cep derinliği/yüzeye yakınlık sinyali verir. | Düşük. |
| APBS v2 full protein | Kaynak PDB üzerinden APBS tekrar çalıştırılır. | Kör tahmin senaryosuna v1'den daha temiz. | Orta, APBS başarısız proteinler olabilir. |
| APBS gradient magnitude | `raw` APBS gridinden finite difference ile türetilir. | Elektrostatik alanın hızlı değiştiği bölgeleri yakalayabilir. | Düşük-orta, ekstra normalization gerekir. |
| APBS positive/negative split | Zaten üretilecek. | Pozitif ve negatif potansiyelin karışmasını engeller. | Düşük. |
| Surface shell mask | `dist_to_surface` üzerinden 0-3 Å, 3-6 Å gibi shell kanalları. | Modeli protein içi hacim yerine yüzeye odaklayabilir. | Düşük. |
| Hydrophobic surface shell | `hydrophobicity` ile surface shell çarpımı. | Pocket surface chemistry sinyalini güçlendirebilir. | Düşük. |

Bu listedeki düşük riskli H5-türevli adaylar cache düzeltme scriptine eklenmiştir: `protein_proximity_exp3`, `protein_near_shell_0_3A`, `protein_near_shell_3_6A`, `hydrophobicity_surface_weighted`, APBS `gradient_magnitude_robust` ve surface-weighted APBS kanalları.

### B. Orta vadeli güçlü adaylar

| Aday feature | Nasıl üretilebilir? | Neden faydalı olabilir? | Risk |
|---|---|---|---|
| fpocket alpha-sphere density | fpocket çıktısındaki alpha sphere merkezleri grid'e rasterize edilir. | Cep geometrisini doğrudan taşır. | Label fpocket ise leakage tartışması yaratabilir; ayrı deney olmalı. |
| fpocket pocket score/rank grid | En yakın pocket score'u grid'e yazılır. | Candidate pocket ranking bilgisini modele verir. | Benchmark ile karışmamalı. |
| SAS point density/protrusion | Protein yüzeyi örneklenir, P2Rank benzeri protrusion hesaplanır. | P2Rank'in güçlü sinyal ailesidir. | Implementasyon süresi orta. |
| Curvature/concavity | Moleküler yüzey veya alpha sphere üzerinden hesaplanır. | Cavity shape için doğrudan bilgi verir. | Surface mesh pipeline gerekir. |
| Residue depth | Her atom/residue için yüzeyden gömülme derinliği. | Binding site derinliği ile ilişkili olabilir. | Orta. |
| Secondary structure | DSSP benzeri araçla helix/sheet/coil rasterize edilir. | Bazı bağlanma bölgeleri yapısal bağlam taşır. | Orta, DSSP dependency. |

### C. Makale genişletmesi için güçlü ama maliyetli adaylar

| Aday feature | Nasıl üretilebilir? | Neden faydalı olabilir? | Risk |
|---|---|---|---|
| Evolutionary conservation | MSA + conservation score, residue skorlarını grid'e rasterize etme. | ConCavity gibi yöntemlerde güçlü sinyal. | Büyük veri için pahalı ve workflow karmaşık. |
| ESM-2 residue embeddings | Protein sequence embeddingleri residue/atom/grid seviyesine aktarılır. | Modern benchmarklarda embedding tabanlı yöntemler güçlü. | Feature boyutu yüksek; model ve storage maliyeti artar. |
| ESM-IF/structure embeddings | 3B yapıdan inverse-folding embeddingleri. | Yapısal bağlamı yoğun biçimde taşır. | Büyük compute, farklı model ailesi. |
| Molecular interaction fields | Farklı probe atomlarıyla enerji gridleri. | SiteHound/EASYMIFs çizgisine yakın, biyofiziksel olarak anlamlı. | Parametre seçimi ve hesaplama maliyeti. |
| Water/solvation features | Su yerleşimi, desolvation, buried unsatisfied polar hesapları. | Binding energetics için güçlü olabilir. | En pahalı ve en karmaşık aday. |
| Ensemble/flexibility features | B-factor, normal mode, MD veya AlphaFold ensemble. | Apo/holo farklarını yakalayabilir. | Tez takvimine göre çok pahalı. |

## Tez İçin Önerilen Feature Stratejisi

### 1. Ana iddia

APBS elektrostatik potansiyel alanı, protein-ligand bağlanma bölgesi tahmininde tek başına ve atomik/kimyasal özniteliklerle birlikte anlamlı sinyal taşır.

Bu iddiayı desteklemek için şu karşılaştırmalar yeterince temizdir:

- `shape_only`
- `apbs_only`
- `shape + apbs`
- `selected_chem`
- `apbs + selected_chem`
- `apbs + shape + selected_chem`
- `apbs + shape + selected_chem + dist_to_surface/hydrophobicity`

### 2. APBS özel ablation

APBS için ayrı olarak şunlar raporlanmalı:

- v1 ligand-proximal APBS ile v2 full-protein APBS farkı
- clip5/clip10/clip20/full_signed150/clip150_signed farkı
- pozitif ve negatif APBS bileşenlerinin ayrı kanal olarak etkisi
- 36, 72:120 ve 161 grid boyutlarının APBS etkisine katkısı

### 3. Leakage kontrolü

`dist_to_ligand` ve `ligand` eğitim girdisi olmayacak. Bunlar sadece:

- görselleştirme,
- sanity check,
- DCC/DCA/DVO hesaplama doğrulaması,
- label üretim kontrolü

için tutulacak.

### 4. Dış veri seti kontrolü

APBS katkısı sadece scPDB içinde gösterilirse zayıf kalabilir. BU48, COACH420 ve PDBBind dış testleri şu sorulara cevap verir:

- APBS-only sinyal dataset dışına taşınıyor mu?
- APBS kimyasal feature'larla birleşince daha stabil oluyor mu?
- 36'lık küçük kutu ile 161'lik APBS çözünürlüğü arasındaki fark genelleniyor mu?
- Top-1, Top-3 ve Top-(n+2) protokolleri altında APBS'nin etkisi devam ediyor mu?

## Tez Savunması İçin Net Cümleler

- "Bu çalışma sadece mevcut atomik feature setine APBS eklemekten ibaret değildir; APBS'nin kaynak seçimi, ölçeklendirmesi, işaretli/işaretsiz normalizasyonu ve tek başına öğrenilebilirliği sistematik olarak ölçülmektedir."
- "Ligand kaynaklı `dist_to_ligand` bilgisi model girdisi değildir; kör tahmin senaryosunda leakage oluşturacağı için auxiliary altında tutulmuştur."
- "Kalasanty ve PUResNet çizgisi atomik grid ve segmentation/residual CNN yaklaşımını kullanırken, bu çalışmanın ayırt edici yönü Poisson-Boltzmann temelli elektrostatik potansiyeli kontrollü ablation içinde incelemesidir."
- "DCC lokasyonu, DVO şekil örtüşmesini, pocket-F1 ise pocket-level başarıyı temsil eder; bu yüzden tek metrik yerine literature-style çoklu metrik raporlanmıştır."
- "APBS v1 eski pipeline uyumluluğu için, APBS v2 ise protein-only inference senaryosuna daha uygun ve temiz temsil için tutulmuştur."

## Kaynaklar

1. Stepniewska-Dziubinska M. M., Zielenkiewicz P., Siedlecki P. "Improving detection of protein-ligand binding sites with 3D segmentation." Scientific Reports, 2020. https://www.nature.com/articles/s41598-020-61860-z
2. Kandel J., Tayara H., Chong K. T. "PUResNet: prediction of protein-ligand binding sites using deep residual neural network." Journal of Cheminformatics, 2021. https://jcheminf.biomedcentral.com/articles/10.1186/s13321-021-00547-7
3. Kandel J. ve ark. "PUResNetV2.0: a deep learning model leveraging sparse representation for improved ligand binding site prediction." Journal of Cheminformatics, 2024. https://jcheminf.biomedcentral.com/articles/10.1186/s13321-024-00865-6
4. Krivak R., Hoksza D. "P2Rank: machine learning based tool for rapid and accurate prediction of ligand binding sites from protein structure." Journal of Cheminformatics, 2018. https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0285-8
5. Jimenez J. ve ark. "DeepSite: protein-binding site predictor using 3D-convolutional neural networks." Bioinformatics, 2017. https://academic.oup.com/bioinformatics/article/33/19/3036/3859178
6. Le Guilloux V., Schmidtke P., Tuffery P. "Fpocket: An open source platform for ligand pocket detection." BMC Bioinformatics, 2009. https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-168
7. Aggarwal R. ve ark. "DeepPocket: Ligand Binding Site Detection and Segmentation using 3D Convolutional Neural Networks." Journal of Chemical Information and Modeling, 2021/2022. https://pubs.acs.org/doi/10.1021/acs.jcim.1c00799
8. Capra J. A. ve ark. "Predicting Protein Ligand Binding Sites by Combining Evolutionary Sequence Conservation and 3D Structure." PLOS Computational Biology, 2009. https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000585
9. Mylonas S. K., Axenopoulos A. "DeepSurf: a surface-based deep learning approach for the prediction of ligand binding sites on proteins." Bioinformatics, 2021. https://academic.oup.com/bioinformatics/article/37/12/1681/6104838
10. Gainza P. ve ark. "Deciphering interaction fingerprints from protein molecular surfaces using geometric deep learning." Nature Methods, 2020. https://www.nature.com/articles/s41592-019-0666-6
11. Jiang ve ark. "A novel protein descriptor for the prediction of drug binding sites." BMC Bioinformatics, 2019. https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-3058-0
12. APBS documentation. "Visualization with PyMOL." https://apbs.readthedocs.io/en/stable/using/examples/visualization-pymol.html
13. "Comparative evaluation of methods for the prediction of protein-ligand binding sites." Journal of Cheminformatics, 2024. https://jcheminf.biomedcentral.com/articles/10.1186/s13321-024-00923-z

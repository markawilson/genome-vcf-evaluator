"""
Gene panels and hg38 coordinates for clinical genomic analysis.
Coordinates are (chrom, start, end) in 1-based positions matching VCF convention.
"""

from __future__ import annotations

# fmt: off

CANCER_GENES = {
    "BRCA1", "BRCA2", "MLH1", "MSH2", "MSH6", "PMS2",
    "APC", "TP53", "PALB2", "CHEK2", "ATM", "CDH1",
    "RB1", "VHL", "PTEN", "STK11", "MEN1", "RET",
    "SMAD4", "NF1", "NF2", "MUTYH",
}

CARDIOVASCULAR_GENES = {
    "LDLR", "APOB", "PCSK9",
    "MYH7", "MYBPC3", "TNNI3", "TNNT2", "TPM1", "MYL2", "MYL3", "ACTC1",
    "SCN5A", "KCNQ1", "KCNH2", "KCNE1", "KCNE2",
    "DSP", "PKP2", "DSG2", "DSC2",
    "FBN1", "FBN2", "COL3A1", "LMNA", "TTN",
}

LONGEVITY_GENES = {
    "APOE", "FOXO3", "CETP", "KL", "TERT", "TERC",
    "SIRT1", "SIRT3", "SIRT6", "MTOR", "IGF1", "IGF1R",
}

MITOCHONDRIAL_GENES = {
    # Nuclear-encoded mitochondrial maintenance genes
    "POLG", "POLG2", "SLC25A4", "TFAM", "SURF1",
    "SCO1", "SCO2", "TWNK", "DGUOK", "COX10",
    # mtDNA-encoded genes (annotated on chrM)
    "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6",
    "MT-CYB",
    "MT-CO1", "MT-CO2", "MT-CO3",
    "MT-ATP6", "MT-ATP8",
    "MT-TL1", "MT-TK", "MT-TE",
}

CONNECTIVE_TISSUE_GENES = {
    # EDS — Classical type
    "COL5A1", "COL5A2",
    # EDS — Classical / Osteogenesis imperfecta overlap
    "COL1A1", "COL1A2",
    # EDS — Vascular (COL3A1 also in cardiovascular panel)
    "COL3A1",
    # EDS — Hypermobile (rare Mendelian cause)
    "TNXB",
    # EDS — Kyphoscoliotic
    "PLOD1", "FKBP14",
    # EDS — Dermatosparaxis
    "ADAMTS2",
    # EDS — Brittle Cornea Syndrome
    "ZNF469", "PRDM5",
    # EDS — Spondylodysplastic
    "B4GALT7", "B3GALT6",
    # EDS — Musculocontractural
    "CHST14", "DSE",
    # Marfan / related (FBN1, FBN2 also in cardiovascular)
    "FBN1", "FBN2",
    # Loeys-Dietz syndrome
    "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "TGFB2", "TGFB3",
    # Osteogenesis imperfecta
    "IFITM5", "SERPINF1",
}

PHARMACOGENOMICS_GENES = {
    # CYP450 drug-metabolizing enzymes
    "CYP2D6", "CYP2C19", "CYP2C9", "CYP3A4", "CYP3A5", "CYP1A2",
    # Anticoagulation
    "VKORC1", "CYP4F2",
    # Chemotherapy toxicity
    "TPMT", "NUDT15", "DPYD", "UGT1A1",
    # Statin transport
    "SLCO1B1",
    # G6PD deficiency (drug-triggered hemolysis)
    "G6PD",
    # Acetylation / isoniazid
    "NAT2",
    # Cystic fibrosis / ivacaftor eligibility
    "CFTR",
}

NEUROLOGICAL_GENES = {
    # Parkinson's disease
    "LRRK2", "GBA", "SNCA", "PRKN", "PINK1", "PARK7",
    # Alzheimer's disease (APOE also in longevity panel)
    "APP", "PSEN1", "PSEN2",
    # ALS / FTD
    "SOD1", "FUS", "TARDBP", "C9orf72", "GRN",
    # Huntington's disease
    "HTT",
    # Frontotemporal dementia / tauopathy
    "MAPT",
    # Spinocerebellar ataxias
    "ATXN1", "ATXN2", "ATXN3",
    # Epilepsy channelopathies
    "SCN1A", "SCN2A", "KCNQ2",
}

DEPRESSION_ANXIETY_GENES = {
    # Serotonin system
    "SLC6A4",   # SERT — serotonin transporter; 5-HTTLPR affects SSRI response
    "HTR2A",    # 5-HT2A receptor — major antidepressant target
    "HTR1A",    # 5-HT1A autoreceptor — SSRI augmentation target
    "HTR2C",    # 5-HT2C receptor — weight/metabolic side-effect risk
    "TPH2",     # tryptophan hydroxylase-2 — rate-limiting brain serotonin synthesis
    # Dopamine/catecholamine system
    "COMT",     # catechol-O-methyltransferase — Val158Met (rs4680); dopamine catabolism
    "MAOA",     # MAO-A — monoamine oxidase; rs6323 affects monoamine levels
    "MAOB",     # MAO-B — selegiline target; dopamine metabolism
    "DRD2",     # D2 receptor — atypical antipsychotic adjunct target
    "DRD4",     # D4 receptor — novelty-seeking; ADHD comorbidity
    # Norepinephrine system
    "SLC6A2",   # NET — norepinephrine transporter; SNRI / TCA target
    # Neuroplasticity
    "BDNF",     # BDNF Val66Met (rs6265) — antidepressant response, memory
    "NTRK2",    # TrkB — BDNF receptor; antidepressant mechanism of action
    # Stress / HPA axis
    "FKBP5",    # FK506-binding protein-5 — glucocorticoid sensitivity; rs1360780
    "NR3C1",    # glucocorticoid receptor — HPA axis dysregulation in depression
    "CRHR1",    # CRH receptor-1 — stress-related depression
    # Glutamate / GABA
    "GRIN2B",   # NMDA receptor GluN2B — ketamine / esketamine target
    # Drug transport across blood-brain barrier
    "ABCB1",    # MDR1 / P-gp — brain drug efflux; affects CNS drug levels
    # Antidepressant pharmacogenomics (drug metabolism)
    "CYP2D6",   # metabolises: fluoxetine, paroxetine, venlafaxine, TCAs, aripiprazole
    "CYP2C19",  # metabolises: escitalopram, citalopram, sertraline, amitriptyline
    # Mood stabiliser / treatment-resistant markers
    "GNB3",     # G-protein β3 — C825T affects antidepressant response
    "ANK3",     # ankyrin-G — mood regulation; lithium response predictor
}

MALE_HORMONE_GENES = {
    # Androgen signalling
    "AR",       # androgen receptor — CAG repeat length ↓ = higher sensitivity
    "SRD5A2",   # 5α-reductase type 2 — testosterone→DHT; finasteride target
    "SRD5A1",   # 5α-reductase type 1 — alternate DHT pathway
    # Steroidogenesis (cholesterol → testosterone)
    "STAR",     # StAR protein — cholesterol transport (rate-limiting step)
    "CYP11A1",  # P450scc — cholesterol→pregnenolone
    "HSD3B1",   # 3β-HSD — pregnenolone→progesterone pathway
    "CYP17A1",  # 17α-hydroxylase/17,20-lyase — DHEA / androstenedione
    "HSD17B3",  # 17β-HSD3 — androstenedione→testosterone (testicular)
    "NR5A1",    # SF-1 — master regulator of steroidogenesis / gonadal development
    # Oestrogen conversion & binding
    "CYP19A1",  # aromatase — testosterone→estradiol; affects male bone & fertility
    "ESR1",     # oestrogen receptor α — male bone density, spermatogenesis
    "ESR2",     # oestrogen receptor β — prostate, testicular function
    # SHBG / free hormone
    "SHBG",     # sex hormone-binding globulin — determines free testosterone fraction
    # HPG axis
    "GNRH1",    # GnRH — hypothalamic driver of LH/FSH
    "KISS1R",   # kisspeptin receptor — upstream master regulator of GnRH pulse
    "LHCGR",    # LH/hCG receptor — Leydig cell stimulation
    "FSHR",     # FSH receptor — Sertoli cell / spermatogenesis
    # Testicular markers
    "INSL3",    # insulin-like factor 3 — Leydig cell function marker
    # Metabolic axis
    "IGF1",     # IGF-1 — anabolic axis; interacts with testosterone
    "CYP11B1",  # 11β-hydroxylase — cortisol synthesis; HPA/testosterone interaction
}

FEMALE_HORMONE_GENES = {
    # Oestrogen signalling
    "ESR1",     # oestrogen receptor α — breast, bone, cardiovascular
    "ESR2",     # oestrogen receptor β — ovary, brain, colon
    "CYP19A1",  # aromatase — major oestrogen source (especially post-menopause)
    "HSD17B1",  # 17β-HSD1 — oestrone→oestradiol in ovary/placenta
    "CYP1B1",   # 4-hydroxylation of oestradiol — genotoxic catechol oestrogen pathway
    "COMT",     # catechol-O-methyltransferase — oestrogen detoxification
    # Progesterone signalling
    "PGR",      # progesterone receptor — PROGINS variant affects cycle / HRT response
    # HPG axis
    "GNRH1",    # GnRH — hypothalamic pulse generator
    "FSHB",     # FSH β subunit — FSH level; menopause timing
    "FSHR",     # FSH receptor — ovarian stimulation response (key for IVF)
    "LHCGR",    # LH/hCG receptor — ovulation trigger
    # Ovarian reserve / AMH axis
    "AMH",      # anti-Müllerian hormone — ovarian reserve marker
    "AMHR2",    # AMH receptor-2 — folliculogenesis regulation
    "GDF9",     # growth differentiation factor-9 — folliculogenesis
    # Adrenal steroidogenesis / androgen excess
    "CYP17A1",  # 17α-hydroxylase — shared step; excess → PCOS-like phenotype
    "CYP11A1",  # P450scc — cholesterol→pregnenolone; steroidogenesis start
    "CYP21A2",  # 21-hydroxylase — deficiency → CAH / PCOS / androgen excess
    # Stress-hormone / cortisol interaction
    "HSD11B1",  # 11β-HSD1 — cortisone↔cortisol in adipose/liver; metabolic risk
    "NR3C1",    # glucocorticoid receptor — stress-cycle interaction
    # Prolactin
    "PRLR",     # prolactin receptor — cycle suppression, galactorrhoea
    # Thyroid (closely coupled with reproductive axis)
    "TSHR",     # TSH receptor — subclinical hypothyroidism common cause of cycle disruption
    # Thrombosis risk (critical for OCP / HRT safety counselling)
    "F5",       # Factor V Leiden (rs6025) — VTE risk ×4–8 with oestrogen
    "F2",       # prothrombin G20210A (rs1799963) — VTE risk with oestrogen
}

ALL_PANELS = {
    "cancer": CANCER_GENES,
    "cardiovascular": CARDIOVASCULAR_GENES,
    "longevity": LONGEVITY_GENES,
    "mitochondrial": MITOCHONDRIAL_GENES,
    "connective_tissue": CONNECTIVE_TISSUE_GENES,
    "pharmacogenomics": PHARMACOGENOMICS_GENES,
    "neurological": NEUROLOGICAL_GENES,
    "depression_anxiety": DEPRESSION_ANXIETY_GENES,
    "male_hormones": MALE_HORMONE_GENES,
    "female_hormones": FEMALE_HORMONE_GENES,
}

# ── Curated common-variant rsIDs per panel ────────────────────────────────────
# These rsIDs are always looked up in addition to the standard clinical filter.
# They are well-established common variants (often MAF > 1%) with documented
# clinical or physiological relevance for their panel that the strict rarity
# filter would otherwise exclude.
#
# Each rsID maps to a gene that IS part of the corresponding panel's gene set.
PANEL_KEY_RSIDS: dict[str, list[str]] = {

    "female_hormones": [
        # ── ESR1 (oestrogen receptor α) ──
        "rs2234693",   # ESR1 PvuII — menstrual cycle, bone density, mood, menopausal symptoms
        "rs9340799",   # ESR1 XbaI  — cardiovascular risk with OCP, HRT response
        "rs3020394",   # ESR1        — breast cancer risk modifier
        # ── ESR2 (oestrogen receptor β) ──
        "rs4986938",   # ESR2 RsaI  — ovarian function, bone
        "rs1256049",   # ESR2        — endometriosis risk
        # ── PGR (progesterone receptor) ──
        "rs1042838",   # PGR PROGINS (Alu insertion proxy) — cycle length, HRT response
        "rs10895068",  # PGR         — progesterone receptor expression
        # ── CYP19A1 (aromatase) ──
        "rs10046",     # CYP19A1     — aromatase activity, post-menopausal oestrogen
        "rs700518",    # CYP19A1     — aromatase activity
        "rs4646",      # CYP19A1     — endometriosis, breast cancer
        # ── CYP1B1 (oestrogen 4-hydroxylation) ──
        "rs1056836",   # CYP1B1 L432V — 4-OH oestradiol (genotoxic pathway)
        # ── COMT (oestrogen catabolism) ──
        "rs4680",      # COMT Val158Met — methylation of catechol oestrogens
        "rs4633",      # COMT           — linked marker
        # ── CYP17A1 ──
        "rs743572",    # CYP17A1 5'-UTR — androgen / oestrogen balance
        # ── FSHR (FSH receptor / IVF response) ──
        "rs6166",      # FSHR N680S     — ovarian stimulation response; key for IVF
        "rs1394205",   # FSHR           — FSH receptor expression
        # ── NR3C1 (glucocorticoid receptor) ──
        "rs41423247",  # NR3C1 N363S    — cortisol sensitivity, HPA axis
        # ── HSD11B1 (cortisol–cortisone interconversion) ──
        "rs846910",    # HSD11B1        — adipose cortisol, metabolic risk
        # ── F5 + F2 (thrombosis — critical for OCP/HRT safety) ──
        "rs6025",      # F5 Factor V Leiden  — VTE risk ×4 on oestrogen
        "rs1799963",   # F2 prothrombin G20210A — VTE risk ×3 on oestrogen
    ],

    "male_hormones": [
        # ── AR (androgen receptor) ──
        "rs6152",      # AR StuI (exon 1) — androgen sensitivity modifier
        # ── SRD5A2 (5α-reductase type 2) ──
        "rs523349",    # SRD5A2 V89L     — DHT production; finasteride context
        "rs9282858",   # SRD5A2 A49T     — 5α-reductase activity
        # ── CYP19A1 (aromatase) ──
        "rs10046",     # CYP19A1         — testosterone→oestradiol conversion
        "rs700518",    # CYP19A1         — aromatase activity
        # ── ESR1 (oestrogen signalling in males) ──
        "rs2234693",   # ESR1 PvuII      — male bone density, spermatogenesis
        "rs9340799",   # ESR1 XbaI       — male cardiovascular risk
        # ── ESR2 ──
        "rs4986938",   # ESR2 RsaI       — prostate, testicular function
        # ── SHBG (free testosterone) ──
        "rs1008805",   # SHBG            — SHBG levels, free testosterone fraction
        "rs6259",      # SHBG D327N      — SHBG binding affinity
        "rs6258",      # SHBG            — SHBG levels
        # ── CYP17A1 ──
        "rs743572",    # CYP17A1         — androgen/oestrogen balance
        # ── COMT ──
        "rs4680",      # COMT Val158Met  — catecholamine / dopamine catabolism
    ],

    "depression_anxiety": [
        # ── COMT (dopamine / mood) ──
        "rs4680",      # COMT Val158Met  — dopamine catabolism; stress response
        "rs4633",      # COMT            — linked marker
        "rs4818",      # COMT            — linked marker
        # ── BDNF (neuroplasticity / antidepressant response) ──
        "rs6265",      # BDNF Val66Met   — antidepressant response, memory, anxiety
        # ── HTR2A (serotonin 2A receptor) ──
        "rs6311",      # HTR2A -1438A/G  — SSRI response, suicide risk
        "rs6313",      # HTR2A           — linked marker
        "rs7997012",   # HTR2A           — antidepressant response
        # ── SLC6A4 (serotonin transporter) ──
        "rs4795541",   # SLC6A4          — 5-HTTLPR region proxy (SSRI response)
        "rs2020936",   # SLC6A4          — SERT expression
        # ── MAOA (monoamine oxidase A) ──
        "rs6323",      # MAOA            — monoamine catabolism; stress response
        "rs909525",    # MAOA            — MAOA uVNTR-linked; activity variant
        # ── TPH2 (serotonin synthesis) ──
        "rs4570625",   # TPH2            — brain serotonin synthesis rate
        "rs11178997",  # TPH2            — TPH2 expression
        # ── DRD2 (dopamine D2 receptor) ──
        "rs1800497",   # DRD2/ANKK1 TaqIA — dopamine reward; antipsychotic response
        # ── FKBP5 (glucocorticoid sensitivity / depression vulnerability) ──
        "rs1360780",   # FKBP5           — childhood trauma → depression interaction
        "rs3800373",   # FKBP5           — stress-sensitisation
        "rs9296158",   # FKBP5           — linked marker
        # ── NR3C1 (glucocorticoid receptor) ──
        "rs41423247",  # NR3C1 N363S     — cortisol hypersensitivity
        # ── GNB3 (G-protein β3; antidepressant response) ──
        "rs5443",      # GNB3 C825T      — antidepressant response; SSRI efficacy
        # ── CRHR1 (CRH receptor; stress response) ──
        "rs110402",    # CRHR1           — HPA axis reactivity
        "rs242924",    # CRHR1           — stress-related depression
        # ── ANK3 (mood regulation / lithium response) ──
        "rs10994336",  # ANK3            — bipolar disorder / mood regulation
        # ── ABCB1 (P-gp; CNS drug penetration) ──
        "rs1045642",   # ABCB1 C3435T    — P-gp expression; antidepressant CNS levels
        "rs2032583",   # ABCB1           — P-gp; SSRI brain penetration
        # ── CYP2D6 (antidepressant metabolism — key star allele SNPs) ──
        "rs3892097",   # CYP2D6 *4       — poor metaboliser (most common Caucasian PM allele)
        "rs35742686",  # CYP2D6 *3       — poor metaboliser (frameshift)
        "rs1065852",   # CYP2D6 *10      — reduced function (common in Asians)
        "rs16947",     # CYP2D6 *2 proxy — near-normal / ultra-rapid context
        # ── CYP2C19 (antidepressant metabolism) ──
        "rs4244285",   # CYP2C19 *2      — poor metaboliser (escitalopram, sertraline)
        "rs4986893",   # CYP2C19 *3      — poor metaboliser (East Asian)
        "rs12248560",  # CYP2C19 *17     — ultra-rapid metaboliser
    ],

    "pharmacogenomics": [
        # ── CYP2D6 star alleles ──
        "rs3892097",   # *4  — poor metaboliser (most common; opioids, TCAs, SSRIs)
        "rs35742686",  # *3  — poor metaboliser (frameshift)
        "rs1065852",   # *10 — reduced function (Asian populations)
        "rs16947",     # *2 proxy
        "rs28371706",  # *6  — poor metaboliser
        "rs5030655",   # *8  — poor metaboliser
        # ── CYP2C19 star alleles ──
        "rs4244285",   # *2  — poor metaboliser (clopidogrel, PPIs, SSRIs)
        "rs4986893",   # *3  — poor metaboliser (East Asian)
        "rs12248560",  # *17 — ultra-rapid metaboliser
        # ── CYP2C9 ──
        "rs1799853",   # *2  — reduced warfarin metabolism
        "rs1057910",   # *3  — reduced warfarin metabolism
        # ── VKORC1 (warfarin sensitivity) ──
        "rs9923231",   # VKORC1 -1639G>A — warfarin dose requirement
        # ── TPMT (thiopurine toxicity) ──
        "rs1800462",   # TPMT *2
        "rs1800460",   # TPMT *3A component
        "rs1142345",   # TPMT *3C
        # ── DPYD (5-fluorouracil toxicity) ──
        "rs3918290",   # DPYD *2A — severe 5-FU toxicity
        "rs55886062",  # DPYD HapB3
        "rs67376798",  # DPYD c.2846A>T
        # ── SLCO1B1 (statin myopathy) ──
        "rs4149056",   # SLCO1B1 *5 — statin-induced myopathy risk
        # ── UGT1A1 (irinotecan toxicity) ──
        "rs887829",    # UGT1A1 *6 proxy
        "rs8175347",   # UGT1A1 *28 (TA repeat) — irinotecan toxicity
        # ── CYP3A5 (tacrolimus metabolism) ──
        "rs776746",    # CYP3A5 *3 — reduced/no CYP3A5 expression (non-expressors)
        # ── CYP4F2 (vitamin K / warfarin) ──
        "rs2108622",   # CYP4F2 V433M — warfarin dose modifier
        # ── NAT2 (isoniazid / caffeine metabolism) ──
        "rs1041983",   # NAT2 — slow acetylator
        "rs1799929",   # NAT2 — slow acetylator
        "rs1799930",   # NAT2 — slow acetylator
    ],

    "longevity": [
        # ── APOE (Alzheimer / cardiovascular / longevity) ──
        "rs429358",    # APOE ε4 allele — Alzheimer risk, reduced longevity
        "rs7412",      # APOE ε2 allele — cardioprotective, longevity-associated
        # ── FOXO3 ──
        "rs2764264",   # FOXO3 — longevity-associated in multiple centenarian studies
        "rs13217795",  # FOXO3 — longevity
        "rs4946936",   # FOXO3 — longevity
        # ── CETP (HDL cholesterol / longevity) ──
        "rs708272",    # CETP TaqIB — HDL levels; longevity in Ashkenazi
        "rs5882",      # CETP        — HDL metabolism
        # ── TERT (telomere length) ──
        "rs10936599",  # TERT — telomere length GWAS
        "rs2736100",   # TERT — telomere length; lung cancer risk
        # ── KL (Klotho / ageing) ──
        "rs9536314",   # KL F352V — longevity-associated heterozygous advantage
        "rs9527025",   # KL — longevity
    ],

    # ── Cancer hereditary risk ─────────────────────────────────────────────────
    # Common / founder moderate-penetrance variants that pass the GWAS/ClinVar
    # significance threshold but have MAF too high for the standard clinical
    # filter (AF ≥ 1%).  High-penetrance BRCA1/2 / Lynch variants are rare
    # enough to clear the filter on their own; the rsIDs below target the
    # moderate-penetrance tier that the filter would otherwise miss.
    "cancer": [
        # ── APC ──
        "rs1801155",   # APC I1307K — colorectal/adenoma risk; ~6% Ashkenazi Jewish
        # ── CHEK2 ──
        "rs17879961",  # CHEK2 I157T — moderate breast/colorectal risk; ~1-2% European
        "rs555607708", # CHEK2 1100delC — European founder, ~0.5-1.4% Northern European
        # ── MUTYH (MAP — MUTYH-associated polyposis) ──
        "rs34612342",  # MUTYH Y165C — MAP founder variant (biallelic → CRC)
        "rs36053993",  # MUTYH G382D — MAP founder variant (biallelic → CRC)
        # ── ATM ──
        "rs1800054",   # ATM S49C — common missense; breast cancer risk modifier
        "rs28904921",  # ATM P1054R — moderate penetrance breast/pancreatic
        # ── PALB2 ──
        "rs152451",    # PALB2 — common intronic tag for GWAS breast-cancer locus
        # ── NBN (Nijmegen breakage syndrome gene) ──
        "rs1805794",   # NBN E185Q — moderate breast/prostate risk (~1% European)
        # ── BRCA1/2 common low-penetrance modifiers ──
        "rs8176318",   # BRCA1 region — common modifier in CIMBA studies
        "rs11200014",  # BRCA2 region — GWAS breast-cancer susceptibility
        # ── TERT / telomere-associated cancer risk ──
        "rs2736100",   # TERT intron — lung cancer + glioma + multiple myeloma GWAS
        # ── TP53 (Li-Fraumeni modifier) ──
        "rs1042522",   # TP53 R72P (Pro72Arg) — common; modifies cancer risk & apoptosis
    ],

    # ── Cardiovascular ─────────────────────────────────────────────────────────
    # Complements the rare LQTS/Brugada/cardiomyopathy variants that clear the
    # standard filter.  These are common variants with well-established effects
    # on lipids, QTc, or arrhythmia penetrance.
    "cardiovascular": [
        # ── PCSK9 ──
        "rs11591147",  # PCSK9 R46L — LOF; cardioprotective (↓LDL ~15%, ~2.6% European)
        "rs562556",    # PCSK9 I474V — common variant affecting LDL-C
        # ── APOB ──
        "rs1367117",   # APOB — common variant; LDL-C GWAS locus
        # ── LPA ──
        "rs10455872",  # LPA KIV-2 tag — elevated Lp(a), CAD risk (~8% European)
        "rs3798220",   # LPA I4399M — elevated Lp(a) and MI risk (~2% European)
        # ── SCN5A (Brugada / conduction modifier) ──
        "rs1805124",   # SCN5A H558R — common modifier of Brugada/LQT penetrance (~20%)
        # ── KCNQ1 (QTc interval) ──
        "rs1057128",   # KCNQ1 — QTc interval GWAS modifier
        # ── KCNH2 / HERG modifier ──
        "rs36210421",  # KCNH2 — common QTc modifier
        # ── APOE (cardiovascular / lipids) ──
        "rs429358",    # APOE ε4 — raised LDL-C; coronary disease risk
        "rs7412",      # APOE ε2 — lower LDL-C; triglyceride modulation
        # ── CYP2C19 (clopidogrel metabolism) ──
        "rs4244285",   # CYP2C19 *2 — clopidogrel poor metaboliser (loss of platelet inhibition)
        "rs12248560",  # CYP2C19 *17 — ultra-rapid metaboliser (increased bleeding risk)
        # ── VKORC1 (warfarin / anticoagulation) ──
        "rs9923231",   # VKORC1 -1639G>A — warfarin dose requirement
        # ── F5 / F2 (thrombosis modifying cardiovascular events) ──
        "rs6025",      # F5 Leiden — VTE; paradoxically affects MI risk in some studies
        "rs1799963",   # F2 prothrombin — VTE / stroke risk
    ],

    # ── Connective tissue ──────────────────────────────────────────────────────
    # True pathogenic FBN1/COL3A1/COL5A1 variants are rare and clear the
    # standard filter.  The rsIDs below target common variants that modulate
    # penetrance or bone/tendon phenotype in connective-tissue conditions.
    "connective_tissue": [
        # ── COL1A1 (EDS / osteoporosis) ──
        "rs1800012",   # COL1A1 Sp1 site — modifies BMD; increases fracture/EDS severity
        # ── COL1A2 ──
        "rs42524",     # COL1A2 — BMD / osteoporosis modifier
        # ── COL5A1 (hypermobile EDS / tendon laxity) ──
        "rs12722",     # COL5A1 3'UTR — tendon stiffness; ACL injury risk
        "rs3196378",   # COL5A1 — tendon / ligament phenotype modifier
        # ── COL3A1 ──
        "rs1800255",   # COL3A1 A698A — associated with aortic aneurysm
        # ── FBN1 (Marfan modifier) ──
        "rs2118181",   # FBN1 — common intronic variant; aortic root diameter modifier
        # ── TNXB (Tenascin-X; hEDS overlap) ──
        "rs2073941",   # TNXB — common variant; hypermobility / EDS modifier
        # ── MMP3 (matrix metalloproteinase; joint laxity) ──
        "rs679620",    # MMP3 — joint laxity / spine degeneration
        "rs3025058",   # MMP3 5A/6A promoter — joint / disc disease
        # ── PLOD1 (lysyl hydroxylase; kyphoscoliotic EDS) ──
        "rs2071436",   # PLOD1 — collagen cross-linking modifier
        # ── ADAMTS10 / ADAMTS17 (Weill-Marchesani) ──
        "rs2280026",   # ADAMTS10 — common variant in fibrillin pathway
        # ── VDR (vitamin D receptor; bone density) ──
        "rs1544410",   # VDR BsmI — BMD / fracture risk modifier
        "rs2228570",   # VDR FokI — BMD modifier; EDS symptom modulator
    ],

    # ── Mitochondrial & nuclear-encoded mitochondrial disease ──────────────────
    # Mitochondrial DNA variants: rsIDs are assigned to chrM positions in
    # dbSNP and should be detectable if the VCF includes chrM/MT calls.
    # Nuclear POLG / TWNK / ANT1 pathogenic variants are typically rare enough
    # to clear the standard filter, but the recurrent alleles below can appear
    # at higher frequency in certain populations.
    "mitochondrial": [
        # ── LHON (Leber hereditary optic neuropathy) — mtDNA primary mutations ──
        "rs28358600",  # MT-ND4  m.11778G>A — most common LHON primary mutation
        "rs28358596",  # MT-ND1  m.3460G>A  — LHON primary mutation
        "rs28358234",  # MT-ND6  m.14484T>C — LHON (better visual prognosis)
        # ── MELAS ──
        "rs199474657", # MT-TL1  m.3243A>G  — MELAS / maternally inherited diabetes + deafness
        # ── MERRF ──
        "rs199474659", # MT-TK   m.8344A>G  — MERRF syndrome
        # ── NARP / Leigh ──
        "rs28357679",  # MT-ATP6 m.8993T>G  — NARP / Leigh syndrome
        "rs28357680",  # MT-ATP6 m.8993T>C  — milder NARP allele
        # ── POLG (nuclear; recurrent pathogenic alleles) ──
        "rs113994096", # POLG A467T — most common recessive POLG pathogenic variant
        "rs111033691", # POLG W748S — second most common recessive allele
        "rs113994097", # POLG G848S — recurrent allele
        # ── TWNK / C10orf2 (Twinkle helicase — mtDNA depletion) ──
        "rs3743916",   # TWNK — common variant in mtDNA copy-number GWAS
        # ── Haplogroup / population-differentiating mtSNPs ──
        "rs2853826",   # MT-CYB m.15326A>G — haplogroup H defining SNP
        "rs28357376",  # MT-ND2 m.5178C>A  — haplogroup D (longevity-associated in Japanese)
    ],

    # ── Neurological ───────────────────────────────────────────────────────────
    # Early-onset Mendelian neurological variants (PSEN1/HTT/etc.) are rare
    # enough to clear the standard filter.  Below are common / founder
    # moderate-penetrance variants relevant to the late-onset / risk-modifier tier.
    "neurological": [
        # ── GBA (Parkinson's / Gaucher) ──
        "rs76763715",  # GBA N370S — most common AJ Gaucher/Parkinson's variant
        "rs421016",    # GBA L444P — Parkinson's risk & Gaucher type 1/3
        "rs75548401",  # GBA E326K — moderate Parkinson's risk modifier
        "rs2230288",   # GBA T369M — Parkinson's risk modifier
        # ── LRRK2 (Parkinson's) ──
        "rs34637584",  # LRRK2 G2019S — Parkinson's; Ashkenazi/N.African founder
        "rs33939927",  # LRRK2 R1441C — Parkinson's European founder
        "rs35870237",  # LRRK2 I2020T — Parkinson's Japanese founder
        # ── SNCA (Parkinson's / DLB) ──
        "rs356219",    # SNCA 3'UTR — Parkinson's GWAS locus (~25% MAF)
        "rs2583988",   # SNCA — Parkinson's risk
        # ── MAPT (tauopathy / FTD / PSP) ──
        "rs17649553",  # MAPT H1 haplotype tag — FTD, PSP, and CBD risk
        "rs8070723",   # MAPT H1 — inversion haplotype tag (corticobasal degeneration)
        # ── APOE (Alzheimer's / DLB) ──
        "rs429358",    # APOE ε4 — major Alzheimer's risk allele
        "rs7412",      # APOE ε2 — protective against Alzheimer's
        # ── APP (Alzheimer's) ──
        "rs63750671",  # APP A673T — protective variant (Icelandic founder; ~0.5% MAF)
        # ── CLU / BIN1 / CR1 (Alzheimer's GWAS loci) ──
        "rs11136000",  # CLU — Alzheimer's susceptibility GWAS
        "rs744373",    # BIN1 — Alzheimer's susceptibility GWAS
        "rs3818361",   # CR1 — Alzheimer's susceptibility GWAS
        # ── C9orf72 (ALS / FTD) ──
        "rs3849942",   # C9orf72 — ALS/FTD locus tag SNP
        # ── HTT modifier ──
        "rs7685686",   # HTT CAG-repeat length modifier locus tag
        # ── TARDBP (ALS) ──
        "rs80356730",  # TARDBP G298S — ALS-associated
        # ── SMN1/SMN2 (SMA copy number tag) ──
        "rs143838139", # SMN1/SMN2 region — SMA copy-number proxy tag
    ],
}

# ── Ancestry-informative markers ───────────────────────────────────────────────
# Genes whose regions contain well-established ancestry-informative markers.
# Queried in unfiltered mode (no clinical significance filter) so common
# population-differentiating alleles are not excluded.
ANCESTRY_GENES = {
    # Pigmentation — skin
    "SLC24A5",   # A111T (rs1426654) — major European skin-lightening variant
    "SLC45A2",   # L374F (rs16891982) — European/East-Asian pigmentation
    "TYR",       # S192Y — melanin synthesis
    "TYRP1",     # melanin polymerisation
    # Pigmentation — eyes and hair
    "OCA2",      # P gene; HERC2 intron regulates OCA2 expression
    "HERC2",     # rs12913832 — blue-eye determinant in Europeans
    "MC1R",      # red hair / pale skin variants
    "TPCN2",     # hair colour modifier
    # Morphology / adaptation
    "EDAR",      # V370A (rs3827760) — East-Asian hair/tooth/sweat-gland signature
    "EDARADD",   # EDAR pathway
    # Carbohydrate metabolism / diet history
    "LCT",       # lactase gene body
    "MCM6",      # contains LCT regulatory element (rs4988235) for lactase persistence
    # Alcohol metabolism
    "ALDH2",     # E504K (rs671) — alcohol flushing, East Asian
    "ADH1B",     # R47H (rs1229984) — alcohol metabolism
    # Immune / malaria adaptation
    "ACKR1",     # DARC/Duffy antigen; Duffy-null near-fixed in sub-Saharan Africans
    # Kidney disease / selection (West African)
    "APOL1",
    # Fatty acid metabolism (reflects ancestral diet — marine vs terrestrial)
    "FADS1",
    "FADS2",
    # Stature / skeletal
    "GDF5",
    # Full mitochondrial genome — scanned unfiltered for haplogroup assignment
    "MT_GENOME",
}

# rsIDs of key ancestry-informative SNPs, looked up directly (bypasses all filters).
# Each rsID is annotated with the population signal it tags.
ANCESTRY_RSIDS = [
    # ── Pigmentation — skin ──────────────────────────────────────────────────
    "rs1426654",    # SLC24A5 A111T   — near-fixed in Europeans (~99%), rare in Africans/E. Asians
    "rs16891982",   # SLC45A2 L374F   — high in Europeans, moderate in E. Asians
    "rs1042602",    # TYR S192Y       — pigmentation modifier
    # ── Pigmentation — eyes ─────────────────────────────────────────────────
    "rs12913832",   # HERC2/OCA2      — major blue-eye locus (Northern European)
    "rs1800407",    # OCA2 R419Q      — green/hazel eye modifier
    # ── Pigmentation — hair ─────────────────────────────────────────────────
    "rs1805007",    # MC1R R151C      — red hair / fair skin (Northern European)
    "rs1805008",    # MC1R R160W      — red hair
    "rs1805009",    # MC1R D294H      — red hair
    # ── Morphology ──────────────────────────────────────────────────────────
    "rs3827760",    # EDAR V370A      — East Asian (thick hair, shovel incisors, more sweat glands)
    # ── Diet / carbohydrate metabolism ──────────────────────────────────────
    "rs4988235",    # MCM6 (LCT)      — lactase persistence, high in N. Europeans & some Africans
    "rs41525747",   # LCT             — lactase persistence, East African pastoralist lineages
    "rs41380347",   # LCT             — lactase persistence, East African lineages
    # ── Alcohol metabolism ───────────────────────────────────────────────────
    "rs671",        # ALDH2 E504K     — alcohol-flush reaction; ~30-40% of East Asians carry it
    "rs1229984",    # ADH1B R47H      — rapid ethanol oxidation; high in East Asians & some S. Asians
    # ── Malaria resistance / immune ─────────────────────────────────────────
    "rs2814778",    # ACKR1/DARC -46T>C — Duffy-null allele; near-fixed in W/C Africans (vivax malaria protection)
    # ── West African specific ────────────────────────────────────────────────
    "rs73885319",   # APOL1 G1a       — West African; associated with kidney disease risk
    "rs60910145",   # APOL1 G1b       — West African
    # ── Fatty acid desaturase ────────────────────────────────────────────────
    "rs174546",     # FADS1           — high-activity allele enriched in Africans; reflects dietary history
    # ── Skeletal / stature ───────────────────────────────────────────────────
    "rs143384",     # GDF5            — height/joint morphology; differs across continents
    # ── APOE — ancestry + disease risk ──────────────────────────────────────
    "rs7412",       # APOE ε2 allele  — varies by ancestry; ε4 elevated in W. Africans
    "rs429358",     # APOE ε4 allele
    # ── Mitochondrial haplogroup anchors — back up the full chrM region scan ──
    # rCRS (the reference sequence) is haplogroup H2a2a1, so H individuals show
    # few chrM variants vs reference; non-H show branch-defining variants below.
    "rs2853826",    # MT-CYB m.15326A>G — haplogroup H marker (Western European)
    "rs28357376",   # MT-ND2 m.5178C>A  — haplogroup D marker (East Asian)
]

# hg38 coordinates (chrom, start, end) — 1-based, inclusive
# Used for coordinate-based gene lookup when VEP/SnpEff annotations are absent.
GENE_COORDINATES: dict[str, tuple[str, int, int]] = {
    # --- Cancer ---
    "BRCA1":   ("chr17",  43_044_295,  43_125_483),
    "BRCA2":   ("chr13",  32_315_474,  32_400_266),
    "MLH1":    ("chr3",   36_993_325,  37_050_846),
    "MSH2":    ("chr2",   47_403_067,  47_559_557),
    "MSH6":    ("chr2",   47_695_572,  47_810_110),
    "PMS2":    ("chr7",   5_970_955,    6_009_130),
    "APC":     ("chr5",  112_707_536, 112_846_239),
    "TP53":    ("chr17",   7_661_779,   7_687_550),
    "PALB2":   ("chr16",  23_603_160,  23_641_310),
    "CHEK2":   ("chr22",  28_687_743,  28_741_956),
    "ATM":     ("chr11", 108_222_484, 108_369_102),
    "CDH1":    ("chr16",  68_737_254,  68_835_516),
    "RB1":     ("chr13",  48_303_747,  48_481_890),
    "VHL":     ("chr3",   10_141_778,  10_153_669),
    "PTEN":    ("chr10",  87_862_625,  87_971_930),
    "STK11":   ("chr19",   1_177_558,   1_228_434),
    "MEN1":    ("chr11",  64_803_516,  64_823_836),
    "RET":     ("chr10",  43_077_029,  43_130_351),
    "SMAD4":   ("chr18",  51_028_394,  51_085_044),
    "NF1":     ("chr17",  31_094_013,  31_377_677),
    "NF2":     ("chr22",  29_603_060,  29_696_491),
    "MUTYH":   ("chr1",   45_329_163,  45_340_893),
    # --- Cardiovascular ---
    "LDLR":    ("chr19",  11_089_862,  11_133_820),
    "APOB":    ("chr2",   21_001_429,  21_044_073),
    "PCSK9":   ("chr1",   55_039_548,  55_064_852),
    "MYH7":    ("chr14",  23_412_572,  23_436_004),
    "MYBPC3":  ("chr11",  47_331_466,  47_374_285),
    "TNNI3":   ("chr19",  55_154_299,  55_163_839),
    "TNNT2":   ("chr1",  201_328_498, 201_371_640),
    "TPM1":    ("chr15",  63_335_034,  63_388_818),
    "MYL2":    ("chr12",  11_060_948,  11_073_498),
    "MYL3":    ("chr3",   46_905_782,  46_913_765),
    "ACTC1":   ("chr15",  35_087_117,  35_097_902),
    "SCN5A":   ("chr3",   38_548_061,  38_649_687),
    "KCNQ1":   ("chr11",   2_447_580,   2_891_637),
    "KCNH2":   ("chr7",  150_930_196, 150_965_832),
    "KCNE1":   ("chr21",  34_394_756,  34_414_818),
    "KCNE2":   ("chr21",  34_385_709,  34_392_299),
    "DSP":     ("chr6",   7_541_801,   7_586_512),
    "PKP2":    ("chr12",  32_779_921,  32_901_863),
    "DSG2":    ("chr18",  29_078_008,  29_131_760),
    "DSC2":    ("chr18",  29_055_794,  29_078_107),
    "FBN1":    ("chr15",  48_408_313,  48_645_709),
    "FBN2":    ("chr5",  127_614_285, 127_875_396),
    "COL3A1":  ("chr2",  188_974_372, 189_011_484),
    "LMNA":    ("chr1",  156_052_823, 156_109_881),
    "TTN":     ("chr2",  178_525_989, 178_807_423),
    # --- Longevity ---
    "APOE":    ("chr19",  44_905_754,  44_909_393),
    "FOXO3":   ("chr6",  108_876_554, 109_010_349),
    "CETP":    ("chr16",  56_961_136,  56_988_436),
    "KL":      ("chr13",  33_587_910,  33_630_069),
    "TERT":    ("chr5",    1_253_167,   1_295_047),
    "TERC":    ("chr3",  169_764_367, 169_765_679),
    "SIRT1":   ("chr10",  67_884_874,  67_918_390),
    "SIRT3":   ("chr11",   4_918_839,   4_939_586),
    "SIRT6":   ("chr19",  4_184_834,   4_196_215),
    "MTOR":    ("chr1",  11_106_535,  11_262_551),
    "IGF1":    ("chr12",  102_395_834, 102_481_941),
    "IGF1R":   ("chr15",  98_648_961,  99_059_852),
    # --- Mitochondrial (nuclear) ---
    "POLG":    ("chr15",  89_858_594,  89_921_077),
    "POLG2":   ("chr17",  62_470_732,  62_483_519),
    "SLC25A4": ("chr4",   185_556_619, 185_562_020),
    "TFAM":    ("chr10",  58_385_274,  58_408_488),
    "SURF1":   ("chr9",   133_351_914, 133_356_956),
    "SCO1":    ("chr17",  11_712_889,  11_722_867),
    "SCO2":    ("chr22",  50_529_959,  50_533_330),
    "TWNK":    ("chr10",  100_990_646, 101_018_241),
    "DGUOK":   ("chr2",   74_121_596,  74_151_498),
    "COX10":   ("chr17",  13_989_979,  14_085_148),
    # --- Connective Tissue / EDS ---
    "COL5A1":  ("chr9",  134_638_362, 134_875_570),
    "COL5A2":  ("chr2",  189_869_459, 190_048_047),
    "COL1A1":  ("chr17",  50_184_967,  50_201_632),
    "COL1A2":  ("chr7",   94_023_977,  94_082_381),
    "TNXB":    ("chr6",   32_039_810,  32_112_134),
    "ADAMTS2": ("chr5",  178_972_997, 179_128_908),
    "PLOD1":   ("chr1",  211_967_488, 212_008_528),
    "FKBP14":  ("chr7",   92_011_779,  92_025_554),
    "ZNF469":  ("chr16",  88_441_801,  88_477_121),
    "PRDM5":   ("chr4",  121_619_261, 121_792_083),
    "B4GALT7": ("chr5",  177_396_547, 177_420_466),
    "B3GALT6": ("chr1",  156_758_444, 156_769_038),
    "CHST14":  ("chr15",  40_700_024,  40_714_598),
    "DSE":     ("chr6",   35_765_642,  35_797_612),
    "TGFBR1":  ("chr9",   99_104_916,  99_186_903),
    "TGFBR2":  ("chr3",   30_648_340,  30_745_054),
    "SMAD2":   ("chr18",  47_808_995,  47_932_556),
    "SMAD3":   ("chr15",  67_063_087,  67_195_068),
    "TGFB2":   ("chr1",  218_346_189, 218_510_011),
    "TGFB3":   ("chr14",  75_958_626,  76_018_083),
    "IFITM5":  ("chr11",     311_015,     312_100),
    "SERPINF1":("chr17",   1_677_503,   1_699_938),
    # --- Pharmacogenomics ---
    "CYP2D6":  ("chr22",  42_522_501,  42_526_883),
    "CYP2C19": ("chr10",  94_762_681,  94_853_000),
    "CYP2C9":  ("chr10",  94_938_657,  95_020_500),
    "CYP3A4":  ("chr7",   99_354_585,  99_465_071),
    "CYP3A5":  ("chr7",   99_245_044,  99_277_621),
    "CYP1A2":  ("chr15",  74_749_498,  74_761_452),
    "CYP4F2":  ("chr19",  15_990_431,  16_009_065),
    "VKORC1":  ("chr16",  31_091_769,  31_096_024),
    "TPMT":    ("chr6",   18_128_544,  18_155_000),
    "NUDT15":  ("chr13",  48_026_944,  48_050_428),
    "DPYD":    ("chr1",   97_543_300,  98_386_615),
    "UGT1A1":  ("chr2",  233_750_129, 233_773_086),
    "SLCO1B1": ("chr12",  21_283_661,  21_395_730),
    "G6PD":    ("chrX",  154_531_392, 154_547_697),
    "NAT2":    ("chr8",   18_391_245,  18_401_875),
    "CFTR":    ("chr7",  117_480_025, 117_668_665),
    # --- Neurological ---
    "LRRK2":   ("chr12",  40_618_813,  40_763_087),
    "GBA":     ("chr1",  155_204_238, 155_214_653),
    "SNCA":    ("chr4",   89_724_099,  89_838_315),
    "PRKN":    ("chr6",  161_347_417, 162_727_497),
    "PINK1":   ("chr1",   20_959_948,  20_978_722),
    "PARK7":   ("chr1",    7_961_654,   7_985_505),
    "APP":     ("chr21",  25_880_550,  26_170_600),
    "PSEN1":   ("chr14",  73_136_416,  73_223_691),
    "PSEN2":   ("chr1",  226_870_612, 226_901_271),
    "SOD1":    ("chr21",  31_659_617,  31_669_877),
    "FUS":     ("chr16",  31_180_213,  31_194_951),
    "TARDBP":  ("chr1",   11_012_654,  11_030_308),
    "C9orf72": ("chr9",   27_546_542,  27_573_863),
    "HTT":     ("chr4",    3_074_680,   3_243_960),
    "MAPT":    ("chr17",  45_894_382,  46_028_417),
    "GRN":     ("chr17",  44_345_302,  44_353_106),
    "ATXN1":   ("chr6",   16_299_402,  16_749_253),
    "ATXN2":   ("chr12", 111_452_718, 111_599_676),
    "ATXN3":   ("chr14",  92_071_012,  92_103_748),
    "SCN1A":   ("chr2",  165_984_570, 166_149_171),
    "SCN2A":   ("chr2",  165_387_882, 165_643_256),
    "KCNQ2":   ("chr20",  62_045_958,  62_131_370),
    # --- Depression / Anxiety / Antidepressant Pharmacogenomics ---
    "SLC6A4":  ("chr17",  28_526_524,  28_620_888),  # SERT serotonin transporter
    "HTR2A":   ("chr13",  46_896_812,  46_998_716),  # 5-HT2A receptor
    "HTR1A":   ("chr5",   63_255_756,  63_258_377),  # 5-HT1A autoreceptor
    "HTR2C":   ("chrX",  113_266_271, 113_560_815),  # 5-HT2C receptor
    "TPH2":    ("chr12",  72_084_246,  72_188_899),  # tryptophan hydroxylase-2
    "COMT":    ("chr22",  19_929_261,  19_959_928),  # catechol-O-methyltransferase
    "MAOA":    ("chrX",   43_654_907,  43_746_824),  # monoamine oxidase A
    "MAOB":    ("chrX",   43_755_721,  43_836_984),  # monoamine oxidase B
    "DRD2":    ("chr11", 113_280_384, 113_419_337),  # dopamine D2 receptor
    "DRD4":    ("chr11",     636_790,     640_370),  # dopamine D4 receptor
    "SLC6A2":  ("chr16",  55_696_671,  55_755_849),  # norepinephrine transporter
    "BDNF":    ("chr11",  27_654_893,  27_722_058),  # brain-derived neurotrophic factor
    "NTRK2":   ("chr9",   84_666_797,  85_023_636),  # TrkB / BDNF receptor
    "FKBP5":   ("chr6",   35_545_834,  35_770_478),  # FK506-binding protein-5
    "NR3C1":   ("chr5",  142_657_495, 142_783_635),  # glucocorticoid receptor
    "CRHR1":   ("chr17",  43_702_487,  43_760_388),  # CRH receptor-1
    "GRIN2B":  ("chr12",  13_536_799,  14_113_530),  # NMDA receptor GluN2B
    "ABCB1":   ("chr7",   87_132_682,  87_343_698),  # MDR1 / P-gp drug efflux pump
    "GNB3":    ("chr12",   6_848_817,   6_861_932),  # G-protein β3 subunit
    "ANK3":    ("chr10",  61_501_259,  62_434_571),  # ankyrin-G / mood regulation
    # --- Hormone Genes (shared Male & Female) ---
    "ESR1":    ("chr6",  151_656_690, 152_129_619),  # oestrogen receptor α
    "ESR2":    ("chr14",  64_227_492,  64_320_025),  # oestrogen receptor β
    "CYP19A1": ("chr15",  51_208_066,  51_338_316),  # aromatase
    "CYP17A1": ("chr10", 104_588_030, 104_601_498),  # 17α-hydroxylase/17,20-lyase
    "CYP11A1": ("chr15",  74_635_459,  74_682_283),  # P450scc cholesterol side-chain cleavage
    "GNRH1":   ("chr8",   25_115_136,  25_121_783),  # gonadotropin-releasing hormone
    "LHCGR":   ("chr2",   48_687_823,  48_879_895),  # LH/hCG receptor
    "FSHR":    ("chr2",   49_101_867,  49_358_428),  # FSH receptor
    "NR5A1":   ("chr9",  124_479_729, 124_534_477),  # SF-1 steroidogenesis master regulator
    "SHBG":    ("chr17",   7_462_060,   7_477_948),  # sex hormone-binding globulin
    # --- Male Hormone specific ---
    "AR":      ("chrX",   67_544_020,  67_730_619),  # androgen receptor
    "SRD5A2":  ("chr2",   31_734_792,  31_802_745),  # 5α-reductase type 2
    "SRD5A1":  ("chr5",    6_640_362,   6_681_958),  # 5α-reductase type 1
    "STAR":    ("chr8",   38_513_034,  38_527_395),  # steroidogenic acute regulatory protein
    "HSD3B1":  ("chr1",  119_458_677, 119_472_561),  # 3β-hydroxysteroid dehydrogenase 1
    "HSD17B3": ("chr9",  100_765_043, 100_855_855),  # 17β-HSD3 (testicular testosterone)
    "KISS1R":  ("chr19",   1_033_697,   1_052_403),  # kisspeptin receptor
    "INSL3":   ("chr19",  47_526_793,  47_529_003),  # insulin-like factor 3 (Leydig)
    "CYP11B1": ("chr8",  143_991_519, 144_001_150),  # 11β-hydroxylase (cortisol synthesis)
    # --- Female Hormone specific ---
    "PGR":     ("chr11", 100_900_750, 101_047_000),  # progesterone receptor
    "HSD17B1": ("chr17",  37_003_432,  37_010_553),  # 17β-HSD1 (ovarian oestradiol)
    "CYP1B1":  ("chr2",   38_294_000,  38_318_000),  # CYP1B1 oestrogen 4-hydroxylation
    "FSHB":    ("chr11",  30_133_924,  30_141_095),  # FSH β-subunit
    "AMH":     ("chr19",   2_249_325,   2_252_773),  # anti-Müllerian hormone
    "AMHR2":   ("chr12",  53_432_949,  53_455_785),  # AMH receptor-2
    "GDF9":    ("chr5",  132_860_523, 132_866_571),  # growth differentiation factor-9
    "CYP21A2": ("chr6",   31_975_382,  31_983_391),  # 21-hydroxylase (CAH)
    "HSD11B1": ("chr1",  209_694_064, 209_754_744),  # 11β-HSD1 cortisone↔cortisol
    "PRLR":    ("chr5",   35_055_012,  35_238_776),  # prolactin receptor
    "TSHR":    ("chr14",  81_395_560,  81_629_428),  # TSH receptor
    "F5":      ("chr1",  169_511_955, 169_586_090),  # coagulation Factor V (Leiden)
    "F2":      ("chr11",  46_709_978,  46_739_906),  # prothrombin (Factor II)
    # --- Ancestry-informative genes ---
    "SLC24A5": ("chr15",  48_413_169,  48_434_589),
    "SLC45A2": ("chr5",   33_944_721,  34_033_844),
    "TYR":     ("chr11",  88_911_281,  89_031_847),
    "TYRP1":   ("chr9",   12_693_953,  12_710_232),
    "OCA2":    ("chr15",  27_622_465,  28_028_897),
    "HERC2":   ("chr15",  28_356_368,  29_000_525),
    "MC1R":    ("chr16",  89_919_487,  89_921_977),
    "TPCN2":   ("chr11",  68_790_906,  68_862_758),
    "EDAR":    ("chr2",  108_879_983, 108_989_888),
    "EDARADD": ("chr1",  236_523_397, 236_606_960),
    "LCT":     ("chr2",  135_787_851, 135_837_186),
    "MCM6":    ("chr2",  135_837_187, 135_890_922),
    "ALDH2":   ("chr12",  111_766_917, 111_817_529),
    "ADH1B":   ("chr4",   99_318_018,  99_363_870),
    "ACKR1":   ("chr1",  159_173_802, 159_176_421),
    "APOL1":   ("chr22",  36_253_070,  36_267_530),
    "FADS1":   ("chr11",  61_797_697,  61_818_102),
    "FADS2":   ("chr11",  61_751_390,  61_799_488),
    "GDF5":    ("chr20",  34_012_145,  34_022_396),
    # --- Mitochondrial (chrM) ---
    "MT-ND1":  ("chrM",     3_307,    4_262),
    "MT-ND2":  ("chrM",     4_470,    5_511),
    "MT-ND3":  ("chrM",    10_059,   10_404),
    "MT-ND4":  ("chrM",    10_760,   12_137),
    "MT-ND4L": ("chrM",    10_470,   10_766),
    "MT-ND5":  ("chrM",    12_337,   14_148),
    "MT-ND6":  ("chrM",    14_149,   14_673),
    "MT-CYB":  ("chrM",    14_747,   15_887),
    "MT-CO1":  ("chrM",     5_904,    7_445),
    "MT-CO2":  ("chrM",     7_586,    8_269),
    "MT-CO3":  ("chrM",     9_207,    9_990),
    "MT-ATP6": ("chrM",     8_527,    9_207),
    "MT-ATP8": ("chrM",     8_366,    8_572),
    "MT-TL1":  ("chrM",     3_230,    3_304),
    "MT-TK":   ("chrM",     8_295,    8_364),
    "MT-TE":   ("chrM",    14_674,   14_742),
    # Full mitochondrial genome — used for ancestry haplogroup scanning
    "MT_GENOME": ("chrM",      1,   16_569),
}

# Build reverse lookup: (chrom, pos) → gene name for unannotated VCFs
def build_interval_index() -> dict[str, list[tuple[int, int, str]]]:
    """Return dict of chrom → sorted list of (start, end, gene)."""
    index: dict[str, list[tuple[int, int, str]]] = {}
    for gene, (chrom, start, end) in GENE_COORDINATES.items():
        index.setdefault(chrom, []).append((start, end, gene))
    # Also index without "chr" prefix for VCFs that omit it
    for gene, (chrom, start, end) in GENE_COORDINATES.items():
        bare = chrom.lstrip("chr") or chrom  # "chr1" → "1", "chrM" → "M"
        if bare != chrom:
            index.setdefault(bare, []).append((start, end, gene))
    for chrom in index:
        index[chrom].sort()
    return index


INTERVAL_INDEX: dict[str, list[tuple[int, int, str]]] = build_interval_index()


def lookup_gene_by_position(chrom: str, pos: int) -> str | None:
    """Return the gene name if pos falls within a known gene region, else None."""
    intervals = INTERVAL_INDEX.get(chrom, [])
    for start, end, gene in intervals:
        if start <= pos <= end:
            return gene
        if start > pos:
            break  # sorted — no point continuing
    return None

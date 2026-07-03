# TCR Clone Assignment

Pipeline em R para atribuição de clonótipos de TCR (T Cell Receptor) a partir do output do **Cell Ranger VDJ** (10x Genomics), com tratamento apropriado de cadeias α/β duplas por barcode celular.

---

## Contexto Biológico

### Rearranjo V(D)J e diversidade de TCR

Os receptores de células T (TCRs) são proteínas heterodiméricas compostas por uma cadeia α (alfa) e uma cadeia β (beta), cada uma gerada por recombinação V(D)J durante o desenvolvimento tímico. Cada T cell madura possui, em teoria, **uma cadeia α funcional e uma cadeia β funcional**, determinando sua especificidade antigênica.

### O problema das cadeias duplas (dual-TCR)

Em single-cell VDJ sequencing, frequentemente observamos **mais de uma cadeia α e/ou β produtiva** para o mesmo cell barcode. Isso pode ter duas causas:

#### 1. Inclusão alélica (fenótipo biológico real)

Durante a seleção tímica, o rearranjo do **locus TCRα** ocorre em ambos os cromossomos. Diferentemente do TCRβ — que tem silenciamento alélico rigoroso —, o TCRα frequentemente escapa do *feedback* negativo antes que um alelo seja silenciado. Resultado:

- **~25-35% das T cells maduras** carregam **duas cadeias α produtivas** (dual-TCRα)
- **~3%** carregam **duas cadeias β produtivas** (dual-TCRβ, raro — ratio α:β ≈ 10:1)

Essas células dual-TCR são biologicamente relevantes: podem ter afinidade alterada por antígenos, participar em respostas cruzadas, e estão implicadas em autoimunidade e tolerância.

**Referências (verificadas via PubMed):**
- Arstila et al. (1999) *Science* 286:958-961 — estimativa fundamental da diversidade TCRαβ humana
- Tuovinen, Salminen & Arstila (2006) *Blood* 108:1601 (PMID: 16926292) — 10-30% de T cells maduras expressam duas cadeias TCR
- Dupic et al. (2019) *PLoS Comput Biol* 15(3):e1006874 (PMID: 30830899) — modelo quantitativo: ~35% das cells expressam ambas as cadeias α
- Schuldt et al. (2015) *J Immunol* (PMID: 26693713) — dual-TCRα supera dual-TCRβ em ~10:1
- Schuldt et al. (2017) *J Immunol* 199(1):33-38 (PMID: 28539428) — dual-TCRα como risco autoimune via limitação de Treg
- Schuldt & Binstadt (2019) *Front Immunol* (PMID: 30670579) — revisão abrangente sobre dual-TCR

#### 2. Artefatos técnicos

| Artefato | Causa | Diagnóstico |
|----------|-------|-------------|
| **Doublets/multiplets** | 2+ células capturadas no mesmo GEM (gel bead-in-emulsion) | Alta proporção de dual-α (>50%) |
| **Ambient RNA** | mRNA livre contaminando o barcode | Cadeia secundária com UMI muito baixo |
| **Index hopping** | Troca de índice durante sequenciação | Padrão aleatório, não reprodutível |

### Diferenciando biológico de artefato

Este pipeline usa três critérios para classificar dual-TCR cells como **confiáveis** vs. **artefato**:

1. **UMI ratio**: a cadeia secundária deve ter ≥10% do UMI da cadeia primária (ajustável)
2. **UMI absoluto**: a cadeia secundária deve ter ≥5 UMIs (ajustável)
3. **Frequência esperada**: se >50% das cells têm dual-α, isso indica doublet contamination

---

## O Pipeline

### Visão geral

```
all_contig_annotations.csv (Cell Ranger VDJ output)
         │
         ▼
   ┌─────────────────┐
   │  1. Filtragem    │  productive + high_confidence + is_cell + TRA/TRB
   └────────┬────────┘
            ▼
   ┌─────────────────┐
   │  2. Ranking      │  Para cada barcode, rankear cadeias por UMI
   │     por UMI      │  (primária vs. secundária)
   └────────┬────────┘
            ▼
   ┌─────────────────┐
   │  3. Dual-TCR     │  Filtrar cadeias secundárias por qualidade
   │     filtering    │  (UMI ratio ≥ 10% + UMI absoluto ≥ 5)
   └────────┬────────┘
            ▼
   ┌─────────────────┐
   │  4. Clone        │  Definir clones por CDR3 (TRB ou TRA+TRB)
   │     assignment   │  Atribuir clone IDs ordenados por frequência
   └────────┬────────┘
            ▼
   ┌─────────────────┐
   │  5. QC + Plots   │  Shannon diversity, Pielou evenness, distribuição
   │                  │  de clone size, dual-TCR rates
   └─────────────────┘
```

### Definição de clone

O pipeline suporta três modos de definição de clone:

| Modo | Descrição | Uso recomendado |
|------|-----------|-----------------|
| `TRB` *(default)* | Clone = CDR3 aa + V gene + J gene da TRB | Padrão da área; TRB é o determinante principal de especificidade |
| `TRA_TRB` | Clone = combinação CDR3 de TRA + TRB | Mais stringente; útil quando TRB sozinho não discrimina |
| `TRB_cdr3_only` | Clone = CDR3 aa da TRB apenas | Mais permissivo; para análises exploratórias |

### Funções principais

| Função | Descrição |
|--------|-----------|
| `load_contigs()` | Carrega e filtra `all_contig_annotations.csv` |
| `rank_chains()` | Rankeia cadeias por UMI para cada barcode |
| `filter_dual_tcr()` | Filtra dual-TCR cells por critérios de qualidade |
| `build_clones()` | Constrói tabela de clones e atribui clone IDs |
| `qc_summary()` | Métricas de QC: dual-TCR rate, Shannon diversity, Pielou |
| `to_screpertoire()` | Formata contigs para integração com `scRepertoire` |
| `run_tcr_pipeline()` | Executa o pipeline completo em um comando |

### Métricas de QC calculadas

- **Dual-α rate**: % de cells com 2+ cadeias α produtivas (esperado: 25-35%)
- **Dual-β rate**: % de cells com 2+ cadeias β produtivas (esperado: ~3%)
- **Shannon diversity index**: $H' = -\sum p_i \ln(p_i)$
- **Pielou evenness**: $J' = H' / \ln(S)$, onde S = nº de clones
- **Top 10 clone share**: % de cells pertencentes aos 10 maiores clones
- **Singleton rate**: % de clones com apenas 1 célula

---

## Instalação

### Pré-requisitos

```r
install.packages(c("tidyverse", "igraph"))
```

### Opcional: scRepertoire (para análise downstream)

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("scRepertoire")
```

---

## Uso

### Exemplo básico (uma amostra)

```r
source("tcr_clone_assignment.R")

results <- run_tcr_pipeline(
  all_contig_csv = "vdj_t/outs/all_contig_annotations.csv",
  output_dir     = "tcr_results",
  clone_def      = "TRB",     # "TRB", "TRA_TRB", ou "TRB_cdr3_only"
  min_umi        = 3,
  min_read       = 10,
  min_ratio      = 0.1,       # cadeia secundária ≥ 10% do UMI da primária
  min_abs_umi    = 5,         # cadeia secundária ≥ 5 UMIs absolutos
  make_plots     = TRUE
)
```

### Múltiplas amostras (batch)

```r
samples <- list(
  sample1 = "sample1/vdj_t/outs/all_contig_annotations.csv",
  sample2 = "sample2/vdj_t/outs/all_contig_annotations.csv"
)

all_results <- lapply(names(samples), function(s) {
  run_tcr_pipeline(samples[[s]], output_dir = paste0("tcr_results/", s))
})
names(all_results) <- names(samples)
```

### Integração com scRepertoire

```r
library(scRepertoire)

# Usar o output formatado pelo pipeline
contigs <- read.csv("tcr_results/contigs_screpertoire_format.csv")
combined <- combineTCR(contigs, samples = NULL, ID = NULL)

# Análises downstream
clonalHomeostasis(combined)
clonalProportion(combined)
clonalOverlap(combined, method = "morisita")
```

---

## Outputs

```
tcr_results/
├── tcr_cells_with_clones.csv       ← Tabela principal: barcode × clonotype
├── clonotype_table.csv             ← Frequência de cada clone
├── dual_chain_filtered.csv         ← Dual-TCR cells aprovadas pelos critérios
├── contigs_screpertoire_format.csv ← Para scRepertoire::combineTCR()
├── clone_size_distribution.pdf     ← Top 50 clones (log scale)
├── chain_status_pie.pdf            ← Proporção single/dual alpha/beta
└── dual_umi_ratio.pdf              ← Histograma UMI ratio (sec/prim)
```

---

## Parâmetros ajustáveis

| Parâmetro | Default | Descrição | Quando ajustar |
|-----------|---------|-----------|----------------|
| `clone_def` | `"TRB"` | Definição de clone | `"TRA_TRB"` para mais stringência; `"TRB_cdr3_only"` para explorar |
| `min_umi` | `3` | UMI mínimo por contig | Aumentar em bibliotecas de alta complexidade |
| `min_read` | `10` | Reads mínimos por contig | Aumentar se suspeitar de contigs espúrios |
| `min_ratio` | `0.1` | UMI ratio sec/prim mínimo | `0.05` para capturar mais dual-TCRs; `0.2` para ser mais rigoroso |
| `min_abs_umi` | `5` | UMI absoluto mínimo na secundária | Aumentar para reduzir falsos positivos |

---

## Estrutura do input

O pipeline espera o arquivo `all_contig_annotations.csv` gerado pelo Cell Ranger VDJ, com as seguintes colunas-chave:

| Coluna | Descrição |
|--------|-----------|
| `barcode` | Cell barcode |
| `is_cell` | Se o contig foi associado a uma célula (vs. background) |
| `high_confidence` | Confiança da assemblagem do contig |
| `chain` | Tipo de cadeia (TRA, TRB, TRG, TRD, IGH, IGL, IGK) |
| `v_gene`, `d_gene`, `j_gene`, `c_gene` | Genes V(D)J e constante |
| `cdr3` | Sequência aminoacidica do CDR3 |
| `cdr3_nt` | Sequência nucleotídica do CDR3 |
| `productive` | Se o rearranjo é produtivo (in-frame, sem stop) |
| `umis` / `umi_count` | Número de UMIs |
| `reads` / `read_count` | Número de reads |

---

## Autores

- **Thiago Y. Oliveira** — FMRP-USP, Bioinformática
- **Lilian Nogueira** — Colaboradora

## Licença

MIT

# MBG — Match-Between-Glycans

**MBG** is a lightweight, modular tool that expands glycopeptide identification in mass spectrometry-based glycoproteomics by inferring additional glycoforms at the MS1 level. It is fully integrated into [FragPipe](https://fragpipe.nesvilab.org/) and designed to work seamlessly after an initial MSFragger-Glyco database search.

---

## Overview

Protein glycosylation is one of the most important and complex post-translational modifications, regulating protein function, stability, and localization. Because glycosylation is a non-template-driven biosynthetic process, a single glycosite typically carries many co-occurring glycoforms that differ by the sequential addition or removal of monosaccharide units. This dilutes the signal of individual glycopeptides, leading low-abundance glycoforms to be missed in standard DDA database searches.

MBG addresses this by exploiting a key property of reverse-phase LC separation: **glycopeptides sharing the same peptide backbone elute at nearly the same retention time (RT), regardless of their glycan**. After an initial MSFragger-Glyco search provides a set of high-confidence glycopeptide identifications, MBG:

1. **Generates candidate glycoforms** — for each identified glycopeptide, it predicts neighboring glycoforms differing by one or more monosaccharide units (e.g., +Hex, +HexNAc, +NeuAc, +Fuc) or user-defined glycan modifications (e.g., NH₄⁺ or Fe³⁺ adducts).
2. **Searches for MS1 evidence** — using IonQuant, it looks for precursor signals at the expected m/z and within a narrow RT/IM window anchored to the parent glycopeptide's observed RT/IM plus a learned monosaccharide-specific shift.
3. **Scores and filters candidates** — a linear discriminant analysis (LDA) model with 7 features (RT/IM shift, mass error, precursor intensity, Y0/Y1 ion relative intensities, isotope envelope KL divergence, glycan shift frequency) separates targets from decoys (decoys use a +11 Da mass offset). FDR is controlled at the precursor level.

The result is an expanded PSM table containing both the original MSFragger-Glyco identifications and the new MBG-inferred glycopeptides, providing a more complete quantitative profile of glycosylation at each glycosite.

---

## Key Features

- **MS1-based inference** — identifies glycopeptides that lack MS2 spectra or have low-quality MS2, using only precursor-level evidence.
- **Learned RT/IM shifts** — monosaccharide-specific RT and ion-mobility (IM) shifts are estimated from the data itself, making MBG robust across different LC gradients and instrument platforms.
- **Adduct and modification recovery** — can recover glycoforms bearing adducts (NH₄⁺, Fe³⁺, Na⁺) or modifications (e.g., phosphorylation for M6P glycans) without expanding the original database search space.
- **Target-decoy FDR control** — rigorous statistical filtering at a user-defined FDR threshold.
- **FragPipe integration** — runs as a one-click step inside the FragPipe glycoproteomics workflow, compatible with label-free and TMT-labeled experiments, DDA and PASEF data.

---

## Workflow

```
Raw files
   │
   ▼
MSFragger-Glyco (database search)
   │
   ▼
PTM-Shepherd + Philosopher (PSM validation, FDR filtering)
   │
   ▼
IonQuant (quantification of identified glycopeptides)
   │
   ▼
MBG  ◄─── this tool
   │   1. Group PSMs by glycosite
   │   2. Estimate per-monosaccharide RT/IM shifts
   │   3. Generate candidate glycoforms (+/- monosaccharides, adducts)
   │   4. Search MS1 via IonQuant within RT/IM windows
   │   5. Score with LDA; apply FDR filter
   │
   ▼
Expanded PSM table (original + inferred glycopeptides)
   │
   ▼
Downstream analysis / FragPipe-Analyst
```

---

## Performance (from the manuscript)

| Dataset | Instrument | Increase in glycopeptide IDs |
|---|---|---|
| Fission yeast | Orbitrap Fusion | +7.5% at 1% FDR; +23.6% at 5% FDR |
| Human plasma (PASEF) | timsTOF HT | +14.6% |
| GBM (CPTAC, TMT-11) | Orbitrap Fusion Lumos | 740 exclusive MBG IDs (4.6% of total) |
| Mouse liver (adducts) | Orbitrap Fusion | +1,234 glycopeptides incl. NH₄⁺ and Fe³⁺ adducts |

Entrapment analysis on the yeast dataset showed an estimated false inference rate of 0.63%.

---

## Requirements

- **Java 11** or later
- **FragPipe** (for integrated use) — MBG is bundled as a tool within FragPipe
- **IonQuant** (bundled as a dependency) — used internally for MS1 feature detection
- **BatMassIO** (bundled) — used for reading raw spectrum files (Thermo `.raw`, Bruker `.d`, mzML, etc.)

### Building from source

```bash
./gradlew jar
```

The output JAR will be at `build/libs/MBG-<version>.jar`.
---
### Installation

**Via FragPipe (recommended):**

MBG is bundled with FragPipe. Simply download and install [FragPipe](https://fragpipe.nesvilab.org/) — no separate installation of MBG is needed.

**Standalone:**

1. Install Java 11 or later.
2. Download the latest `mbg-<version>.jar` from the [Releases](https://github.com/Nesvilab/MBG/releases) page, or build from source with `./gradlew jar`.

**Typical install time:** Less than 5 minutes (as part of FragPipe installation on a standard desktop computer). Building from source takes approximately 1–2 minutes.

---

## Usage

MBG is typically run automatically within FragPipe. For standalone use:

```bash
java -jar mbg-<version>.jar --match [args]
```

### Required arguments

| Argument | Description |
|---|---|
| `--psm <path>` | Path to the input PSM file (FragPipe `psm.tsv` format) |

### Optional arguments

| Argument | Default | Description |
|---|---|---|
| `--manifest <path>` | — | Path to the FragPipe `.fp-manifest` file mapping raw file names to full paths |
| `--residuedb <path>` | built-in | Path to a custom glycan residue definitions file |
| `--glycanmoddb <path>` | built-in | Path to a custom glycan modification definitions file (for adducts, etc.) |
| `--maxq <float>` | `0.01` | Maximum glycan q-value to accept as a high-confidence input glycoPSM |
| `--minpsms <int>` | `2` | Minimum number of PSMs required for a glycan to be used for RT/IM shift estimation |
| `--minglycans <int>` | `2` | Minimum number of distinct glycans observed per peptide to enable inference |
| `--rttol <float>` | `0.4` | RT tolerance (minutes) for matching inferred glycoforms |
| `--imtol <float>` | `0.05` | Ion mobility tolerance (V·s·cm⁻²) for PASEF data |
| `--mztol <int>` | `5` | Mass tolerance (ppm) for MS1 matching |
| `--fdr <float>` | `0.05` | FDR threshold for accepting inferred glycopeptides |
| `--nopasef <bool>` | `false` | Set to `true` for non-PASEF (Orbitrap) data |
| `--runtmt <bool>` | `false` | Set to `true` for TMT-labeled experiments |
| `--toaddresiduals <str>` | — | Comma-separated list of monosaccharides or glycan compositions to add as candidate shifts (e.g., `HexNAc,Fuc,HexNAc(1)Hex(1)`) |
| `--expanddb <int>` | `0` | Number of additional rounds of iterative inference on newly found peaks |
| `--maxskips <int>` | `0` | Number of missed peaks allowed during iterative expansion |
| `--numthreads <int>` | `4` | Number of threads for parallel processing |
| `--allowchimeric <bool>` | `false` | Allow chimeric spectra when searching for supporting MS2 |

### Example

```bash
java -jar mbg-0.3.7.jar --match \
    --psm /data/experiment/psm.tsv \
    --manifest /data/experiment/fragpipe.fp-manifest \
    --maxq 0.01 \
    --minpsms 1 \
    --rttol 0.4 \
    --fdr 0.05 \
    --nopasef true
```
---

## Demo: Fission Yeast Dataset

This demo walks through running MBG on the fission yeast glycoproteomics dataset (PXD005565), a well-characterized benchmark used in the manuscript.

### 1. Download the demo data

Download the fission yeast dataset from ProteomeXchange:
- **Accession:** [PXD005565](https://www.ebi.ac.uk/pride/archive/projects/PXD005565)
- Download the `.raw` files (3 files).
- (optional) Convert to mzML using [MSConvert](https://proteowizard.sourceforge.io/).

### 2. Run FragPipe with MSFragger-Glyco and MBG

1. Open FragPipe and load the mzML files.
2. Select the **glyco-N-HCD** workflow.
3. Click **Database** panel to set the protein database to the *S. pombe* proteome (download from UniProt, Proteome ID UP000002485, add decoys and contaminants). 
4. Click **Glyco** panel set the glycan database as **Mouse_N-glycans-1670-pGlyco**. Check **Run Glycoform Inference** and adjust MBG setting through **MS1 Glycoform Inference** block as the following:
   - Max Glycan q-value: 0.01
   - Min PSMs: 1
   - Min Glycans: 2
   - MBG FDR: 0.01
   - Max Inference Step: 2
   - Max Inference Skips: 0
   - Allow inferring chimeric spectra: checked
   - Residuals to add: Hex(1),NH4(1)

5. Click **Quant(MS1)** panel, in **Feature dectection and peak tracing** block, set m/z tolerance to 10, RT tolerance to 0.4

6. Run Fragpipe

### 3. Run MBG out of Fragpipe

Run Fragpipe as the above parameters without checking **Run Glycoform Inference**, after the FragPipe search completes, run MBG on the output:

```bash
java -jar mbg-0.3.7.jar --match \
    --psm /path/to/yeast_output/psm.tsv \
    --manifest /path/to/yeast_output/fragpipe.fp-manifest \
    --toaddresiduals Hex(1),NH4(1) \
    --maxq 0.01 \
    --minpsms 1 \
    --minglycans 2 \
    --rttol 0.4 \
    --mztol 10 \
    --fdr 0.01 \
    --expanddb 2 \
    --nopasef true
```


### 4. Expected output

After MBG completes, you should see:

- **Updated `psm.tsv`** — the original PSMs plus MBG-inferred glycopeptides. At 5% FDR, expect approximately a **7.5%** (or 23.6% with 5% MBG FDR) increase in unique glycopeptide identifications.
- **`rt_shifts.csv`** — retention time shifts per glycan transition (e.g., HexNAc(2)Hex(8) → HexNAc(2)Hex(9)). Most high-mannose transitions should fall within ~0.4 minutes.
- MBG-inferred entries are flagged in the PSM table and can be distinguished from MSFragger-Glyco identifications for downstream analysis.

### 5. Expected run time

- **FragPipe search (MSFragger-Glyco + PTM-Shepherd + IonQuant):** approximately 15–20 minutes for 3 raw files on a standard desktop.
- **MBG step only:** approximately 2–5 minutes on the same hardware.


---

## Output

MBG appends inferred glycopeptide PSMs to the input PSM table and writes:

- **Updated `psm.tsv`** — original PSMs plus MBG-inferred entries (marked with an `MBG` flag in the `Hyperscore` or source columns)
- **`rt_shifts.csv`** / **`im_shifts.csv`** — per-monosaccharide RT and IM shift distributions (useful for QC and plotting)
- **Skyline modifications file** — definitions for any novel glycan compositions inferred by MBG, for import into Skyline

---

## Citation

> Shen J, Polasky DA, Jager S, Yu F, Heck AJR, Reiding KR, Nesvizhskii AI.
> **Expanding Glycopeptide Identification with Match-Between-Glycans in FragPipe.**
> *Manuscript in preparation*, 2026.
> Data: PRIDE PXD074575.

---

## License

See [LICENSE](LICENSE) for details.

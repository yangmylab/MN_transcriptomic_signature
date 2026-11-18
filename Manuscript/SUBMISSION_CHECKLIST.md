# Submission checklist — Paper 1 (*EBioMedicine*)

**Status:** prepared 2026-11-17; every item below is either *agent-complete* or an explicit
**user sign-off**. Nothing in this list may be auto-filled with invented detail.

---

## 1. Files to upload

| # | File | State | Notes |
|---|---|---|---|
| 1 | `Manuscript/Paper1_manuscript.md` → `.docx` | agent-complete | EBioMedicine structure; numbers from the audited pipeline |
| 2 | `Manuscript/Cover_letter.md` → `.docx` | agent-complete | 519 words; author list and reviewers are placeholders |
| 3 | `output/figures/Figure1–6.pdf` | agent-complete | vector, 0 non-embedded fonts |
| 4 | `output/figures/supplementary/FigureS1–S3.pdf` | agent-complete | vector, 0 non-embedded fonts |
| 5 | `Manuscript/references.tsv` | agent-complete | 44 entries, every DOI machine-verified |
| 6 | `Manuscript/IF_validation_plan.md` | internal | not uploaded; drives the P0 wet-lab work |
| 7 | `Manuscript/DATA_REQUESTS.md` | internal | not uploaded; drives the B2/B7 request |

## 2. USER sign-offs required before submission

| # | Item | Why it blocks | Current placeholder |
|---|---|---|---|
| U1 | **Author list and order** | The manuscript cannot be submitted with a placeholder author block; CRediT, affiliations and the reference set can all shift with it | `[further authors to be confirmed before submission]` |
| U2 | **Affiliations and ORCIDs** | Required by the submission system for every author | only the corresponding author's two affiliations are recorded |
| U3 | **CRediT contributions** | EBioMedicine requires a contribution statement per author | Meng only |
| U4 | **Funding statement** | Currently "No specific funding"; if any grant supported this, the grant number must appear in the Summary's *Funding* line **and** the Acknowledgements | "No specific funding" |
| U5 | **Declaration of interests** | Must be true for every author, not only the corresponding author | "The author declares no competing interests" |
| U6 | **Corresponding author contact** | Must be the address the journal will use | `m.lingzhang@hotmail.com`; `lzmeng@gxams.org.cn` |
| U7 | **Three suggested reviewers** | Optional but strongly recommended; must not have co-authored recently | none |
| U8 | **Data availability wording** | Current text promises deposition "at acceptance" — confirm you will deposit, or switch to on-request only | deposition statement |
| U9 | **Ethics statement** | Confirm no institutional review was required for public de-identified data at your institution | stated as not required |
| U10 | **Word / figure caps** | Unconfirmed (publisher pages blocked, P-50). Check the submission system's live limits when you create the submission | Summary 371 words; 6 main figures + 3 supplementary |

## 3. Agent items still open (do on request)

- [ ] Re-derive the separability gate under cross-platform harmonisation (05a/05b) and propagate
      the resulting numbers to Methods, Results, Discussion and the Summary.
- [ ] Confirm the EBioMedicine caps from the live submission system once you are logged in.
- [ ] Renumber citations if U1 changes the reference set.
- [ ] Rebuild every `.docx` after any text change:
      `pandoc <in>.md --reference-doc=reference/tnr_fangsong_redcode.docx -o <out>.docx`.

## 4. Verification gates already satisfied

- Six verification scripts exit 0 (51 / 86 / 27 / 42 / 117 / 103 assertions).
- All nine figure plates: vector PDF, 0 non-embedded fonts, 300 dpi JPEG previews.
- Manuscript numbers are script-generated; a drift fails verification before it reaches the text.

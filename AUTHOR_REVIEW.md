# Author review before submission

This is a completed technical author-review package, not a submitted or accepted paper. The request to prepare a manuscript was authorized on 29 October 2026. The author confirmed the saved name/affiliation/email and stated: **No specific funding and no conflicts of interest.**

## Material journal-policy issue

The official Optica policy checked on 29 October 2026 states that AI tools may enhance/edit original text but may not generate manuscript text. AI-generated code and figures must be documented, including tool information and how it was used. See:

https://opg.optica.org/content/author/portal/item/review-general-policies

The current narrative was drafted with OpenAI Codex assistance; it must not be represented as compliant original author prose merely after proofreading. The author needs to prepare the final narrative in their own words, with accurate acknowledgment of computational assistance, and should seek editorial clarification if the permitted boundary remains unclear. This requirement was reported to the author during preparation. No submission has occurred.

OpenAI Codex assisted the mathematical investigation, implementation, computational checks, figure-generation scripts, and drafting. The figures are deterministic MATLAB plots of frozen numerical records, not generated illustrations or invented measurements. MATLAB scripts and exact certificates are retained so the plots and claims can be checked without a generative model. The local Codex CLI version observed during preparation was `0.165.0-alpha.1.2`; an exact model build identifier covering all historical research sessions was not retained and must not be invented. A final disclosure should accurately reflect the tools and versions that can be established from the author's records.

## Scientific review priorities

1. Review the problem's practical relevance: the small fixed library, normalized XYZ metric, slope class, and continuous relative-SPD assumptions. These are the main external-validity limitations.
2. Review the proofs, particularly the uncancelled residual norm in model uncertainty and actual-white normalization in the shared-SPD bound. The explicit calibrated predictor includes target white-response information.
3. Retain the complete lower-protocol coverage. The main contribution is the certified count/precision/regularity result and robust quantitative regions, not a new convex framework or general recovery theorem.
4. Keep unresolved frontier and tolerance intervals and sufficient uncertainty budgets accurately described. More optimization is not required to support the claims presently made.
5. Complete or appropriately qualify the remaining closest full-text comparisons. Worthey (1988/1994), Trussell (1991), Yang--Stark (1994), Belforte--Tay (1993), Pronzato--Walter (1990), and recent Forsythe--Funt conference papers remain incomplete full-text comparisons in the research record. The manuscript makes no first-in-literature claim and cites Worthey only for the broad historical topic. Abstracts, metadata, code, or secondary descriptions are not full-text audits.
6. Assess whether a dataset or physical demonstration is desirable for the chosen journal. Such additions would have a different role from the existing universal mathematical guarantees and are not silently included here.

## Release and submission

The associated GitHub repository is private, consistent with the standing research release restrictions. Provide reviewer access or authorize an appropriate public release before using a public-data availability statement. The current statement explicitly reports private author-review status. A stable archived data release and DOI can be added later; none has been invented.

No journal submission, public release, paid service, automated restart, additional optimization campaign, or sub-agent delegation was performed. Publication likelihood cannot be assigned from numerical verification alone.

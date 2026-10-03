---
name: lab-protocols
description: Look up wet-lab standard operating procedures from Addgene and Thermo Fisher — viral vector production/titration, cloning, transfection, ELISA, flow cytometry, cell culture, staining, and general lab practice. Use when a task asks how to perform a molecular or cell biology bench protocol.
---

# Lab Protocols

Find the applicable protocol at its publisher and cite the version used.
This skill provides discovery guidance and links; it does not bundle protocol
full text.

## Sources

- [Addgene protocols](https://www.addgene.org/protocols/): plasmid preparation,
  cloning, transfection, viral vectors and antibody applications.
- [Thermo Fisher protocols](https://www.thermofisher.com/us/en/home/references/protocols.html):
  cell culture, staining, ELISA, flow cytometry and reagent-specific methods.
  Match the procedure to the product/catalog number and its current manual.
- [protocols.io](https://www.protocols.io/): discover published protocols and
  select the relevant version; use the API as described below when requested.

## Usage

1. Search the publisher's site for the topic, organism, assay or reagent.
   Open the original protocol or manual rather than relying on a search snippet.
2. Record its title, publisher, canonical URL or DOI, revision/version when
   available, and access date. State when no revision is supplied.
3. Retain the published conditions and distinguish them from any requested
   adaptation. Cite the specific protocol, not only the collection homepage.
   If the source cannot be accessed, report that limitation instead of inventing
   steps or treating the result as an empty search.

Summarize and link to the source. Publisher content remains subject to its own
terms; BioV's MIT licence does not cover it. Do not add protocol full-text
collections to the plugin or Python distribution.

These are reference procedures, not executable code — adapt volumes and
concentrations to the user's actual experiment only when asked.

For procedures outside these collections (including oligo annealing, Golden
Gate assembly and chondrogenic aggregate assays), obtain the applicable
manufacturer or published method and its revision; do not reuse the old
Biomni hard-coded recipe merely because its function name matches. Record
reagent/enzyme identity, input assumptions and the source for quantitative
conditions. The presence of a reference protocol does not mean an experiment
was performed or its yield/purity was measured.

For protocols.io, use its [documented API](https://apidoc.protocols.io/) with
the user's authorized token to search public protocols and retrieve a chosen
protocol/version. Keep the native JSON, ID, DOI/URL, version and access date;
do not embed credentials in outputs or treat an API error as an empty search.
Use a browser for discovery when no API access was requested or available.
No BioV search wrapper or credential fallback is required.

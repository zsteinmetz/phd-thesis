# Agricultural Plastic Covers—Source of Plastic Debris in Soil?

Submitted dissertation thesis for the partial fulfillment of the requirements
for a Doctor of Natural Sciences

*By Zacharias Steinmetz*

## Summary

The use of agricultural plastic covers has become common practice for
its agronomic benefits such as improving yields and crop quality,
managing harvest times better, and increasing pesticide and water use
efficiency. However, plastic covers are suspected of partially breaking
down into smaller debris and thereby contributing to soil pollution with
microplastics. A better understanding of the sources and fate of plastic
debris in terrestrial systems has so far been hindered by the lack of
adequate analytical techniques for the mass-based and polymer-selective
quantification of plastic debris in soil. The aim of this dissertation
was thus to assess, develop, and validate thermoanalytical methods for
the mass-based quantification of relevant polymers in and around
agricultural fields previously covered with fleeces, perforated foils,
and plastic mulches. Thermogravimetry/mass spectrometry (TGA/MS) enabled
direct plastic analyses of 50 mg of soil without any sample preparation.
With polyethylene terephthalate (PET) as a preliminary model, the method
limit of detection (LOD) was 0.7 g kg<sup>−1</sup>. But the missing
chromatographic separation complicated the quantification of polymer mixtures.
Therefore, a pyrolysis-gas chromatography/mass spectrometry (Py-GC/MS) method
was developed that additionally exploited the selective solubility of
polymers in specific solvents prior to analysis. By dissolving polyethylene
(PE), polypropylene (PP), and polystyrene (PS) in a mixture of
1,2,4-trichlorobenzene and *p*-xylene after density separation,
up to 50 g soil became amenable to routine plastic analysis. Method
LODs were 0.3–2.2 mg kg<sup>−1</sup>, and the recovery of 20 mg kg<sup>−1</sup>
PE, PP, and PS from a reference loamy sand was 86–105%. In the reference
silty clay, however, poor PS recoveries, potentially induced by the additional
separation step, suggested a semi-quantitative evaluation of PS. Yet, the
new solvent-based Py-GC/MS method enabled a first exploratory screening of
plastic-covered soil. It revealed PE, PP, and PS contents above LOD in
six fields (6% of all samples). In three fields, PE levels of
3–35 mg kg<sup>−1</sup> were associated with the use of 40 μm thin perforated
foils. By contrast, 50 μm PE films were not shown to induce plastic
levels above LOD. PP and PS contents of 5–19 mg kg<sup>−1</sup>
were restricted to single observations on four sites and potentially
originated from littering. The results suggest that the short-term use
of thicker and more durable plastic covers should be preferred to limit
plastic emissions and accumulation in soil. By providing mass-based
information on the distribution of the three most common plastics in
agricultural soil, this work may facilitate comparisons with modeling
and effect data and thus contribute to a better risk assessment and
regulation of plastics. However, the fate of plastic debris in the
terrestrial environment remains incompletely understood and needs to be
scrutinized in future, more systematic research. This should include the
study of aging processes, the interaction of plastics with other organic
and inorganic compounds, and the environmental impact of biodegradable
plastics and nanoplastics.

## Included Files and Folders

<dl>
  <code>README.md</code>
  <dd>
    Readme file
  </dd>
  <code>thesis.pdf</code>
  <dd>
    Dissertation as PDF file
  </dd>
  <code>thesis.tex</code>
  <dd>
    Dissertation as editable LaTeX source file including all thesis contents
    from the following subdirectories: <code>frontbackmatter/</code>,
    <code>chapters/</code>, <code>appendices/</code>
  </dd>
  <code>proposal.pdf</code>
  <dd>
    Thesis proposal as PDF file
  </dd>
  <code>proposal.tex</code>
  <dd>
    Thesis proposal as editable LaTeX source file
  </dd>
  <code>references.bib</code>
  <dd>
    Complete bibliographical data of the references cited
  </dd>
  <code>figures/</code>
  <dd>
    Figures
  </dd>
    <code>tables/</code>
  <dd>
    Tables
  </dd>
    <code>supplements/</code>
  <dd>
    Supplemental (raw) data and code
  </dd>
    <code>cv/</code>
  <dd>
    Curriculum vitae
  </dd>
</dl>

## Build

Building PDF files from LaTeX sources requires a LaTeX distribution
installed on your computer, such as TeX Live or MiKTeX. See
[The LaTeX Project](https://www.latex-project.org/get/) for a brief overview of
LaTeX distributions and installation instructions.

This dissertation was built on Manjaro Linux with TeX Live 2021 using `latexmk`.
The following code chunk produces the PDF:

```shell
latexmk -f -pdf thesis.tex
```

After building, you may want to clean your environment with `latexmk -c`.

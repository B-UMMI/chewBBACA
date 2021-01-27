import React from "react";
import ReactMarkdown from "react-markdown";
import gfm from "remark-gfm";
import classes from "./Header.css";

// Material-UI ExpansionPanel components
import Accordion from "@material-ui/core/Accordion";
import AccordionDetails from "@material-ui/core/AccordionDetails";

export default function Header() {
  const markdown = `
  # Schema Evaluator

  Provides summary charts that allows users to explore:
  - The diversity (number of alleles) at each locus (**Panel A**);
  - The variation of allele mode sizes per locus (**Panel B**);
  - Summary statistics (minimum allele size in blue, minimum allele size in orange and median allele size in green) for each locus (**Panel C**);
  - The loci size distribution (**Panel D**);
  - The presence of alleles that are not CDSs (when evaluating schemas called by other algorithms) (**Panel E**).

  Users are able to select an **individual locus** to be analysed, by clicking on:
  - each **point (locus)** of the **Locus Statistics** and **Locus Size Variation charts**;
  - the **locus name** on the **Locus** column of the **CDS Analysis** table or (if it exists) the **Loci with high variability** table.

  By selecting a locus the following will be displayed:
  - 2 charts (histogram and scatterplot) containing an **analysis of the allele sizes**;
  - a table with **summary statistics** of the alleles;
  - a **Neighbor Joining tree** built by clustal based on the mafft alignment;
  - a **multiple sequence alignment** of the alleles produced by mafft.
  `;

  return (
    <div>
      <Accordion defaultExpanded>
        <AccordionDetails>
          <div
            className={classes.mainPaper}
            style={{ width: "100%", height: "100%" }}
          >
            <ReactMarkdown plugins={[gfm]} children={markdown} />
          </div>
        </AccordionDetails>
      </Accordion>
    </div>
  );
}

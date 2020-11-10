import React from "react";
import ReactMarkdown from "react-markdown";
import gfm from "remark-gfm";
import classes from "./Header.css";

// Material-UI components
import Button from "@material-ui/core/Button";
import Typography from "@material-ui/core/Typography";
import { createMuiTheme, MuiThemeProvider } from "@material-ui/core/styles";

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

  After the panels, a **dropdown menu** allows users to select an **individual locus** to analysed.
  By selecting a locus the following will be displayed:
  - 2 charts (histogram and scatterplot) containing an **analysis of the allele sizes**;
  - a table with **summary statistics** of the alleles;
  - a **multiple sequence alignment** of the alleles.
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

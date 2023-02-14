import AccordionMUI from './AccordionMUI';


const Header = () => {
  const markdown = `
  # Schema Evaluator

  Provides summary charts that allows users to explore:
  - The diversity (number of alleles) at each locus (**Panel A**);
  - The variation of allele mode sizes per locus (**Panel B**);
  - Summary statistics (minimum allele size in blue, maximum allele size in orange and median allele size in green) for each locus (**Panel C**);
  - The loci size distribution (**Panel D**);
  - The presence of alleles that do not comply with the parameters used to define the schema (this is particularly relevant when evaluating schemas created by other algorithms and that are adapted for use with chewBBACA) (**Panel E**).

  Users are able to select an **individual locus** to be analysed, by clicking on:
  - each **point (locus)** of the **Locus Statistics** and **Locus Size Variation charts**;
  - the **locus name** on the **Locus** column of the **CDS Analysis** table, the **Loci with high variability** table or the **Loci with only 1 allele** table.

  By selecting a locus the following will be displayed:
  - 2 charts (histogram and scatterplot) containing an **analysis of the allele sizes**;
  - a table with **summary statistics** of the alleles;
  - a **Neighbor Joining tree** based on the mafft alignment;
  - a **multiple sequence alignment** of the alleles produced by mafft.
  `;

  return (
    <div>
      <AccordionMUI 
        summaryText='Report Description' 
        detailsText={markdown} 
        expanded={false} >
      </AccordionMUI>
    </div>
  );
}

export default Header;

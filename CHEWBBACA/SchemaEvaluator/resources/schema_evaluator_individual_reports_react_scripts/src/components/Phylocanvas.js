import React from "react";

// Phylocanvas imports
import Phylocanvas from "phylocanvas";

// Phylocanvas plugins imports
import contextMenu from "phylocanvas-plugin-context-menu";
import scalebar from "phylocanvas-plugin-scalebar";
// import history from 'phylocanvas-plugin-history';

// Register Phylocanvas plugins
Phylocanvas.plugin(contextMenu);
Phylocanvas.plugin(scalebar);
// Phylocanvas.plugin(history);

export class PhylogeneticTree extends React.Component {
    constructor(props) {
      super(props);
  
      this.tree = null;
    }
  
    setTreeAttributes = () => {
      // Set the tree type to rectangular
      this.tree.setTreeType(this.props.treeType);
      this.tree.disableZoom = !this.props.zoom;
      this.tree.setTextSize(15);
      this.tree.resizeToContainer();
      this.tree.draw();
    };
  
    componentDidMount = () => {
      this.tree = Phylocanvas.createTree(this.node, { nodeAlign: true });
      let newickTree = this.props.newickString;
  
      // Load the newick tree from props
      this.tree.load(newickTree);
      this.setTreeAttributes();
      // tree.load("((_R_CC0067_NODE_1_length_10197_cov_734.723715_pilon:0.0006430160208146113,(CC0061_k77_1_flag_1_multi_4641.2458_len_10267_pilon:0.00953411150529446,_R_CC0116_NODE_1_length_10196_cov_675.686135_pilon:0.004676216709352007)56:0.003523710625083518)100:0.5306604729360941,Spike_NODE_2_length_10199_cov_229.021834_pilon:0.29987078286758695,_R_Spike_NODE_1_length_10319_cov_2021.808436_pilon:0.2754770121464598);");
    };
  
    render() {
      if (this.node) {
        this.setTreeAttributes();
      }
  
      return (
        <div style={{ width: "100%", zIndex: "1" }}>
          <div
            style={{ height: "700px" }}
            ref={(node) => (this.node = node)}
          ></div>
        </div>
      );
    }
  }
  
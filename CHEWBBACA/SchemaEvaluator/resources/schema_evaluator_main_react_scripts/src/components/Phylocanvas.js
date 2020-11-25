import PropTypes from "prop-types";
import React from "react";
import Phylocanvas from "phylocanvas";
import { treeTypes } from "phylocanvas";
import contextMenu from "phylocanvas-plugin-context-menu";
import scalebar from 'phylocanvas-plugin-scalebar';
import _keys from "lodash.keys";

// Register Phylocanvas plugins
Phylocanvas.plugin(contextMenu);
Phylocanvas.plugin(scalebar);

export default class PhylocanvasComponent extends React.Component {
  static propTypes = {
    className: PropTypes.string,
    data: PropTypes.string,
    style: PropTypes.object,
    treeType: PropTypes.oneOf(_keys(treeTypes)),
  };

  componentDidMount() {
    this.tree = Phylocanvas.createTree(this.refs.phyloCanvasDiv);
    this.componentDidUpdate({});
  }

  componentDidUpdate(prevProps) {
    const props = this.props;
    if (prevProps.data !== props.data) {
      this.tree.load(props.data);
    }
    if (prevProps.treeType !== props.treeType) {
      this.tree.setTreeType(props.treeType);
    }
  }

  render() {
    const { className, style } = this.props;
    return <div ref="phyloCanvasDiv" style={style} className={className} />;
  }
}

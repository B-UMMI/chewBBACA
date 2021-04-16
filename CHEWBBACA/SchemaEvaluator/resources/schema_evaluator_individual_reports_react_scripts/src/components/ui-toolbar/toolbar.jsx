import React from 'react';
import PropTypes from 'prop-types';

// import {fitToViewer} from '../features/zoom';
import IconCursor from './icon-cursor.jsx';
import IconPan from './icon-pan.jsx';
import IconZoomIn from './icon-zoom-in.jsx';
import IconZoomOut from './icon-zoom-out.jsx';
// import IconFit from './icon-fit';
import ToolbarButton from './toolbar-button.jsx';

export default function Toolbar({tool, value, onChangeValue, onChangeTool, activeToolColor, position}) {

  let handleChangeTool = (event, tool) => {
    onChangeTool(tool);
    event.stopPropagation();
    event.preventDefault();
  };

  // let handleFit = event => {
  //   onChangeValue(fitToViewer(value, SVGAlignX, SVGAlignY));
  //   event.stopPropagation();
  //   event.preventDefault();
  // };

  let isHorizontal = ['top', 'bottom'].indexOf(position) >= 0;

  let style = {
    //position
    position: "absolute",
    transform: ['top', 'bottom'].indexOf(position) >= 0 ? "translate(-50%, 0px)" : "none",
    top: ['left', 'right', 'top'].indexOf(position) >= 0 ? "5px" : "unset",
    left: ['top', 'bottom'].indexOf(position) >= 0 ? "50%" : ('left' === position ? "5px" : "unset"),
    right: ['right'].indexOf(position) >= 0 ? "5px" : "unset",
    bottom: ['bottom'].indexOf(position) >= 0 ? "5px" : "unset",

    //inner styling
    backgroundColor: "rgba(19, 20, 22, 0.90)",
    borderRadius: "2px",
    display: "flex",
    flexDirection: isHorizontal ? "row" : "column",
    padding: isHorizontal ? "1px 2px" : "2px 1px"
  };

  return (
    <div style={style} role="toolbar">
      <ToolbarButton
        toolbarPosition={position}
        active={tool === 'none'}
        activeColor={activeToolColor}
        name="unselect-tools"
        title="Selection"
        onClick={ event => handleChangeTool(event, 'none') }>
        <IconCursor/>
      </ToolbarButton>

      <ToolbarButton
        toolbarPosition={position}
        active={tool === 'pan'}
        activeColor={activeToolColor}
        name="select-tool-pan"
        title="Pan"
        onClick={ event => handleChangeTool(event, 'pan') }>
        <IconPan/>
      </ToolbarButton>

      <ToolbarButton
        toolbarPosition={position}
        active={tool === 'zoom-in'}
        activeColor={activeToolColor}
        name="select-tool-zoom-in"
        title="Zoom in"
        onClick={ event => handleChangeTool(event, 'zoom-in') }>
        <IconZoomIn/>
      </ToolbarButton>

      <ToolbarButton
        toolbarPosition={position}
        active={tool === 'zoom-out'}
        activeColor={activeToolColor}
        name="select-tool-zoom-out"
        title="Zoom out"
        onClick={ event => handleChangeTool(event, 'zoom-out') }>
        <IconZoomOut/>
      </ToolbarButton>

      {/* <ToolbarButton
        toolbarPosition={position}
        active={false}
        activeColor={activeToolColor}
        name="fit-to-viewer"
        title="Fit to viewer"
        onClick={ event => handleFit(event) }>
        <IconFit/>
      </ToolbarButton> */}
    </div>
  )
}

Toolbar.propTypes = {
  tool: PropTypes.string.isRequired,
  onChangeTool: PropTypes.func.isRequired,
  value: PropTypes.object.isRequired,
  onChangeValue: PropTypes.func.isRequired,

  //customizations
  position: PropTypes.oneOf(['top', 'right', 'bottom', 'left']),
  // SVGAlignX: PropTypes.oneOf([ALIGN_CENTER, ALIGN_LEFT, ALIGN_RIGHT]),
  // SVGAlignY: PropTypes.oneOf([ALIGN_CENTER, ALIGN_TOP, ALIGN_BOTTOM]),
  activeToolColor: PropTypes.string
};

Toolbar.defaultProps = {
  position: 'right',
  // SVGAlignX: ALIGN_LEFT,
  // SVGAlignY: ALIGN_TOP,
  activeToolColor: '#1CA6FC'
};

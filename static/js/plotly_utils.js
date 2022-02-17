function getPlotlySize() {
  return 300;
};

function getPlotlyMargin() {
  return {l: 40, r: 40, b: 40, t: 40, pad: 0};
};

function getbackgroundColor() {
  return "#dee2e6";
};

function getPlotlyConfig(){
  var REMOVE_MODEBARS = [
      //'resetScale2d', 'resetCameraDefault3d',
      //'zoom2d', 'pan2d', 'select2d', 'lasso2d',
      'zoomIn2d', 'zoomOut2d', 'autoScale2d',
      'zoom3d', 'pan3d', 'orbitRotation', 'tableRotation', 'handleDrag3d',
      'resetCameraLastSave3d', 'hoverClosest3d', 'hoverClosestCartesian', 'hoverCompareCartesian',
      'hoverClosestGl2d', 'hoverClosestPie', 'toggleHover', 'resetViews',
      'sendDataToCloud', 'toggleSpikelines', 'resetViewMapbox', //'toImage',
  ];
  return {
    displaylogo: false,
    modeBarButtonsToRemove: REMOVE_MODEBARS,
    responsive: true,
    displayModeBar: true,
    toImageButtonOptions: {
      format: 'png',
      filename: "downloaded-plot",
      scale: 4,
    },
  };
};
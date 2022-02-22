
function plotHistogram() {
  const bins = dom.get('bin-select').value;
  const similarities = JSON.parse(dom.get("file-upload").dataset.similarities);
  const styleMap = {color: "#198754", size: 3, label: "inside similarity range"};
  const traces = [{
    x: similarities.filter(x => x !== 1),
    type: 'histogram',
    marker: {color: styleMap.color},
    //label: v.label,
    //hovertext: "hii",
    xbins: {start: 0, size: 1/bins, end: 1},
  }];
  const layout = {
      barmode: "overlay",
      hovermode: 'closest',
      margin: getPlotlyMargin(),
      width: getPlotlySize() - 50,
      height: getPlotlySize() - 50,
      xaxis: {title: 'Similarity score', showgrid: false},
      yaxis: {title: 'Counts', showgrid: false},
      showlegend: false,
  };
  dom.hideLoading('histogram-plot');
  Plotly.newPlot('histogram-plot', traces, layout, getPlotlyConfig());
};



function plotPCA() {

  const smiles0 = dom.get("file-upload").dataset.selected;
  const pca = JSON.parse(dom.get("file-upload").dataset.pca);
  const smilesList = JSON.parse(dom.get("file-upload").dataset.smilesList);
  const similarities = JSON.parse(dom.get("file-upload").dataset.similarities);
  // get index of currently selected molecule
  const selectedIdx = smilesList.indexOf(smiles0);

  // get text to show on hover
  var hoverText = [];
  for (let i=0; i < pca.length; i++){
    hoverText.push(`<b>${smilesList[i]}</b><br>Similarity: <b>${Math.round(similarities[i]*100)/100}</b>`);
  };

  const traces = [
    {
      x: pca.map(x => x[0]),
      y: pca.map(x => x[1]),
      type: 'scattergl',
      mode: 'markers',
      //text: smilesList,
      text: hoverText,
      marker: {
          color: similarities,
          colorscale: "Jet",
          reversescale: true,
          cmax: 0.4,
          cmin: 0,
          size: 4,
          opacity: 1,
          line: {color: 'white', width: 0.5, opacity: 1},
          colorbar: {thickness: 10},
      },
    },
    {
      x: [pca[selectedIdx][0]],
      y: [pca[selectedIdx][1]],
      text: `Selected molecule:<br><b>${smiles0}</b>`,
      type: 'scatter',
      mode: 'markers',
      marker: {
          color: 'rgba(100,0,0,0)',
          size: 10,
          opacity: 1,
          line: {color: 'gray', width: 1, opacity: 1},
      },
    },
    {
      x: [pca[selectedIdx][0]],
      y: [pca[selectedIdx][1]],
      text: `Selected molecule:<br><b>${smiles0}</b>`,
      type: 'scatter',
      mode: 'markers',
      marker: {
          color: 'rgba(100,0,0,0)',
          size: 30,
          opacity: 1,
          line: {color: 'gray', width: 1, opacity: 1},
      },
    },
    {
      x: [pca[selectedIdx][0]],
      y: [pca[selectedIdx][1]],
      text: `Selected molecule:<br><b>${smiles0}</b>`,
      type: 'scatter',
      mode: 'markers',
      marker: {
          color: 'rgba(100,0,0,0)',
          size: 50,
          opacity: 1,
          line: {color: 'gray', width: 1, opacity: 1},
      },
    },
  ];
  const layout = {
      hovermode: 'closest',
      margin: getPlotlyMargin(),
      width: getPlotlySize() + 100,
      height: getPlotlySize(),
      xaxis: {title: 'PCA-1'},
      yaxis: {title: 'PCA-2'},
      showlegend: false,
  };

  dom.hideLoading('pca-plot');
  Plotly.newPlot('pca-plot', traces, layout, getPlotlyConfig());
  dom.get('pca-plot').dataset.initialized = "true"

  /*
  myPlot.on('plotly_click', function(data){
    var pts = '';
    for(var i=0; i < data.points.length; i++){
        pts = 'x = '+data.points[i].x +'\ny = '+
            data.points[i].y.toPrecision(4) + '\n\n';
    }
    alert('Closest point clicked:\n\n'+pts);
  });
  */

  dom.get('pca-plot').on('plotly_click', function(plotData){
      for(var i=0; i < plotData.points.length; i++){
          var newIdx = plotData.points[i].pointNumber;
          var newCurveNum = plotData.points[i].curveNumber;
      };       
      if (newCurveNum === 0){
        const newMolecule = smilesList[newIdx];
        selectMolecule(newMolecule);
      };
  });
};

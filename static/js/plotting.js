
function plotScatter() {
  var data = JSON.parse(document.getElementById("file-upload").dataset.data);
  var xcol = document.getElementById('x-axis-select').value;
  var ycol = document.getElementById('y-axis-select').value;
  var xvals = data.map(x => x[xcol]);
  var yvals = data.map(x => x[ycol]);

  var plot = [{
      x: xvals,
      y: yvals,
      type: 'scatter',
      mode: 'markers',
      //text: colorVals.map(x => `VALUE: ${x}`),
      marker: {
          size: 8,
          color: 'darkorange',
          //colorscale: 'Jet',
          symbol: 'circle',
          line: {'width': 1, 'color': 'white'},
          //color: colorVals,
      },
  }];
  var layout = {  
      xaxis: {title: xcol, automargin: true},
      yaxis: {title: ycol, automargin: true},
      hovermode:'closest',
      //hoverlabel: {font: {size: 10}},
      width: getPlotlySize()*1.5,
      height: getPlotlySize(),
      margin: getPlotlyMargin(),
      showlegend: false,
      plot_bgcolor: getbackgroundColor(),
      paper_bgcolor: getbackgroundColor(),
  };
  //if ([...new Set(colorVals)].length > 1) {
  //    plot[0]['marker']['colorbar'] = {len: 0.6, thickness: 10, thicknessmode: 'pixels',};
  //};
  Plotly.newPlot('scatter-plot-div', plot, layout, getPlotlyConfig());

};









function plotHistogram() {

  var data = JSON.parse(document.getElementById("file-upload").dataset.data);
  var columnToPlot = document.getElementById('histogram-select').value;
  var vals = data.map(x => x[columnToPlot]);

  var plot = [{x: vals, type: 'histogram', marker: {color: 'darkorange'},}];

  var layout = {  
    xaxis: {title: columnToPlot, automargin: true},
    yaxis: {title: "Counts", automargin: true},
    hovermode:'closest',
    //hoverlabel: {font: {size: 10}},
    width: getPlotlySize(),
    height: getPlotlySize(),
    margin: getPlotlyMargin(),
    showlegend: false,
    plot_bgcolor: getbackgroundColor(),
    paper_bgcolor: getbackgroundColor(),
  };

  Plotly.newPlot('histogram-plot-div', plot, layout, getPlotlyConfig());

  vals = vals.map(x => parseFloat(x));
  var tableData = [
    ['Min', d3.min(vals)],
    ['Max', d3.max(vals)],
    ['Min index', d3.minIndex(vals)],
    ['Max index', d3.maxIndex(vals)],
    ['Sum', d3.sum(vals)],
    ['Mean', d3.mean(vals)],
    ['Median', d3.median(vals)],
    ['Mode', d3.mode(vals)],
    ['Variance', d3.variance(vals)],
    ['Std. deviation', d3.deviation(vals)],
  ];
  var tableString = `<table class="table table-sm table-borderless" style="font-size: 0.8rem;"><tbody>`;
  for (let r of tableData){
    tableString += `<tr><th>${r[0]}</th><td>${r[1]}</td><tr>`
  };
  tableString += `</tbody></table>`;
  document.getElementById('histogram-stats-div').innerHTML = tableString;
};





// from a dataset in the form of an array of
// objects, create a map of all the column correlations

function getCorrelationMap(data){
  var columns = Object.keys(data[0]);
  r2Map = [];
  var x;
  var y;
  var r2;
  for (let c of columns){
    mapRow = [];
    for (let c2 of columns){
      x = data.map(x => parseFloat(x[c]));
      y = data.map(x => parseFloat(x[c2]));
      r2 = getR2(x, y);
      mapRow.push(r2);
    };
    r2Map.push(mapRow);
  };
  return [r2Map, columns];
};





function plotCorrelationMap() {

  var data = JSON.parse(document.getElementById("file-upload").dataset.data);
  var [r2Map, columns] = getCorrelationMap(data);

  var plot = [{
      z: r2Map,
      x: columns,
      y: columns,
      type: 'heatmap',
      hoverongaps: false,
      colorscale: 'Portland',
    }];
  
  var layout = {
    xaxis: {automargin: true},
    yaxis: {scaleanchor: "x", automargin: true},
    hovermode:'closest',
    width: getPlotlySize()*1.5,
    height: getPlotlySize()*1.5,
    margin: getPlotlyMargin(),
    showlegend: false,
    plot_bgcolor: getbackgroundColor(),
    paper_bgcolor: getbackgroundColor(),
  };

  Plotly.newPlot('correlation-map-div', plot, layout, getPlotlyConfig());
};












document.getElementById('secret-btn').addEventListener('click', () => {



  const [win, windowBody] = createWindow(title="My new plot (this title is editable)");
  
  const myId  = getRandomId();
  const plotId = `${myId}-plot`;

  const plotDiv = document.createElement('div');
  plotDiv.id = plotId;

  windowBody.appendChild(plotDiv);
  plotResponsiveScatter(plotId);

  function resizePlot() {
    if (document.getElementById(plotId)){
      Plotly.relayout(plotId, {
        'height': windowBody.clientHeight - 10, // these offset values must be exactly twice the size of the window content padding!
        'width': windowBody.clientWidth - 10,
        'xaxis.autorange': true,
        'yaxis.autorange': true,
      });
    };
  };
   
   new ResizeObserver(resizePlot).observe(win);


}, false);




function plotResponsiveScatter(plotId) {
  var traces = [{
    x: [-2, 2, 3, 4],
    y: [-0.2, 0.9, 1.5, 3.6],
    mode: 'lines+markers',
    type: 'scatter',
    marker: {
        size: 12,
        color: 'darkorange',
        symbol: 'circle',
        line: {'width': 1, 'color': 'darkorange'},
    },
  }];
  var layout = {  
    xaxis: {title: "X-axis", automargin: true},
    yaxis: {title: "Y-axis", automargin: true},
    autosize: true,
    hovermode:'closest',
    margin: {l: 60, r: 20, b: 100, t: 0, pad: 0},
    showlegend: false,
    plot_bgcolor: "white",
    paper_bgcolor: "white",
  };
  Plotly.newPlot(plotId, traces, layout, getPlotlyConfig());
};



// Create a movable, stretchable window.
// Inputs are title (string), and content (DOM element)
function createWindow(title="My window title") {
  // create window
  const win = document.createElement('div');
  Object.assign(win.style, {
    position: 'absolute',
    top: '300px',
    left: '30px',
    backgroundColor: 'white',
    boxShadow: '3px 3px 10px 3px rgba(0,0,0,0.5)',
    //zIndex: 100,
    resize: 'both',
    minWidth: '250px',
    minHeight: '50px',
    borderRadius: "5px",
    overflow: "hidden",
  });

  // add window to main
  document.getElementsByTagName('main')[0].appendChild(win);

  // create window top title bar
  const titleBar = document.createElement('div');
  titleBar.classList.add('d-flex', 'justify-content-between');
  win.appendChild(titleBar);
  Object.assign(titleBar.style, {
    padding: '5px',
    cursor: 'move',
    //marginBottom: '5px',
    borderBottom: '1px solid #ccc',
    backgroundColor: '#dee2e6',
  });
  // create icon on title bar
  const titleBarIcon = document.createElement('span');
  titleBarIcon.innerHTML = `<i class="bi bi-arrows-move"></i>`;
  titleBar.appendChild(titleBarIcon);
  // create window title on title bar
  const titleDiv = document.createElement('span');
  titleDiv.contentEditable = "true";
  Object.assign(titleDiv.style, {
    fontWeight: "bold",
    outline: "0px solid transparent",
    cursor: "text",
  });
  titleDiv.innerHTML = title;
  titleBar.appendChild(titleDiv);


  // create div for buttons
  const btnDiv = document.createElement('div');
  titleBar.appendChild(btnDiv);

  // create dock window button on title bar
  const dockBtn = document.createElement('button');
  dockBtn.type = 'button';
  dockBtn.classList.add("btn", "btn-sm", "bg-transparent", "text-dark");
  dockBtn.innerHTML = "<i class='bi bi-arrows-expand'></i>";
  dockBtn.dataset.expanded = "true";
  btnDiv.appendChild(dockBtn);
  //Object.assign(dockBtn.style, {cursor: 'pointer',});
  dockBtn.addEventListener('click', () => {
    if (dockBtn.dataset.expanded === "true"){
      // if window is currently expanded, minimize it
      dockBtn.dataset.height = win.clientHeight;
      //dockBtn.dataset.width = win.clientWidth;
      Object.assign(win.style, {
        height: '50px',
        //width: '400px',
      });
      dockBtn.dataset.expanded = "false";
    } else {
      // if window is currently minimized, expand it
      Object.assign(win.style, {
        height: `${dockBtn.dataset.height}px`,
        //width: `${dockBtn.dataset.width}px`,
      });
      dockBtn.dataset.expanded = "true";
    };
  }, false);

  // create close window button on title bar
  const closeBtn = document.createElement('button');
  closeBtn.type = 'button';
  closeBtn.classList.add("btn", "btn-sm", "bg-transparent", "text-dark")
  closeBtn.innerHTML = "<i class='bi bi-x-lg'></i>";
  btnDiv.appendChild(closeBtn);
  //Object.assign(closeBtn.style, {cursor: 'pointer'});
  closeBtn.addEventListener('click', () => {win.remove()}, false);

  // create main window content
  const windowBody = document.createElement('div');
  win.appendChild(windowBody);
  Object.assign(windowBody.style, {
    padding: '5px',
    height: '100%',
  });

  // add draggable logic to window title bar
  titleBar.addEventListener('mousedown', function(evt) {
    // record where the window started
    var real = window.getComputedStyle(win);
    var winX = parseFloat(real.left);
    var winY = parseFloat(real.top);
    // record where the mouse started
    var mX = evt.clientX;
    var mY = evt.clientY;
    // drag the window until the mouse button comes up
    document.body.addEventListener('mousemove', drag, false);
    document.body.addEventListener('mouseup', function() {
      document.body.removeEventListener('mousemove', drag, false);
    }, false);
    function drag(evt){
      // add difference between where the mouse is now
      // versus where it was last to the original positions
      win.style.left = winX + evt.clientX-mX + 'px';
      win.style.top  = winY + evt.clientY-mY + 'px';
    };
  }, false);

  return [win, windowBody];
};
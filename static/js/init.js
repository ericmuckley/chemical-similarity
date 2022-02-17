

function showDataLoading(){
  document.getElementById('loading-div').style.display = "block";
  document.getElementById('plot-row').style.display = "none";
  document.getElementById('table-div').innerHTML = "";
  document.getElementById('data-stats-div').innerHTML = "";
};


document.getElementById('use-example-csv-btn').addEventListener("click", function() {
  showDataLoading();
  var dataUrl = "https://raw.githubusercontent.com/ericmuckley/datasets/master/wikipedia-molecule-properties.csv";
  d3.csv(dataUrl).then(function(data) {
    load_data(data);
  });
}, false);

// when file is selected
document.getElementById("file-upload").addEventListener("change", function() {
  showDataLoading();

  // configure file reader action once file is finished being read
  var reader = new FileReader();
  reader.onload = function(e) {
    var contents = e.target.result;
    var data = d3.csvParse(contents);
    load_data(data);
  };    
  // read the selected file
  var file = this.files[0];
  //if (typeof file === 'Blob'){
  reader.readAsText(file);
  //};
  
}, false);  


function load_data(data) {

  data = data.slice(0, 10000);


  document.getElementById("file-upload").dataset.data = JSON.stringify(data);

  document.getElementById('data-stats-div').innerHTML = `
    Dataset contains ${data.length} rows, ${Object.keys(data[0]).length} columns
  `;

  var [table, thead, tbody] = createFilledTable(
    data=data,
    classes=['table', 'table-sm', 'table-bordered'],
    divId='table-div',
  );
  Object.assign(table.style, {
    'font-size': '0.75rem',
  });



  // build histogram
  populateSelect('histogram-select', Object.keys(data[0]));
  plotHistogram();
  [`histogram-select`].forEach(x => {
      var el = document.getElementById(x);
      el.addEventListener('change', event => {
        plotHistogram();
      });
  });
  // build scatter plot
  var cols = Object.keys(data[0]);
  populateSelect('x-axis-select', cols, selectedIndex=0);
  populateSelect('y-axis-select', cols, selectedIndex=1);
  plotScatter();
  [`x-axis-select`, `y-axis-select`].forEach(x => {
      var el = document.getElementById(x);
      el.addEventListener('change', event => {
        plotScatter();
      });
  });



  // build correlation heatmap
  plotCorrelationMap();



  /*

  // show column details
  var tableHeaders = ['Column Name', 'Data type', 'Missing values'];
  var [table, thead, tbody] = createFilledTable();
  var row = tbody.insertRow();
  var rowString = ``;
  for (let header of tableHeaders){
    rowString += `<th>${header}</th>`;
  };
  row.innerHTML = rowString;
  for (let c of cols){
    var row = tbody.insertRow();
   row.innerHTML = `<th>${c}</th><td>${typeof parseFloat(data[0][c])}</td><td>${3}</td>`;
  };
  document.getElementById('column-details-div').innerHTML = ""
  document.getElementById('column-details-div').appendChild(table);
  */






  // show the plots
  document.getElementById('loading-div').style.display = "none";
  document.getElementById('plot-row').style.display = "block";
};


// when help button is clicked
document.getElementById("help-btn").addEventListener("click", function() {
  document.getElementById('help-alert').style.display = 'block';
  window.scrollTo({top: 0, behavior: 'smooth'});
}, false);
document.getElementById("close-help-alert-btn").addEventListener("click", function() {
  document.getElementById('help-alert').style.display = 'none';
}, false);
// when user clicks on example file link
dom.get('use-example-csv').addEventListener("click", function() {
  showLoading('loading-div');
  var dataUrl = "https://raw.githubusercontent.com/ericmuckley/datasets/master/qm9-small.csv";
  d3.csv(dataUrl).then(function(data) {
    load_data(data);
  });
}, false);

// when file is selected from upload field
dom.get("file-upload").addEventListener("change", function() {
  showLoading('loading-div');
  // configure file reader action once file is finished being read
  var reader = new FileReader();
  reader.onload = function(e) {
    var contents = e.target.result;
    var data = d3.csvParse(contents);
    load_data(data);
  };    
  // read the selected file
  var file = this.files[0];
  reader.readAsText(file);
}, false);  


function load_data(data) {
  if (data.length > 10000){
    alert("Thats a lot of data. We will use only the first 10,000 rows to ensure performance.")
    data = data.slice(0, 10000);
  };
  dom.get("file-upload").dataset.data = JSON.stringify(data);
  console.log(data);

  getSimilarities()


  hideLoading('loading-div');
};








function showLoading(id){
  dom.make('div', {
    parent: id,
    classes: 'text-center m-2',
    children: [
      ['div', {classes: "spinner-border text-success", styles:{width: '5rem', height:"5rem"}}]
    ],
  });
};

function hideLoading(id){
  dom.get(id).innerHTML = "";
};

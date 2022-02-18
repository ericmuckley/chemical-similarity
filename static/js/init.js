// when user clicks on example file link
dom.get('use-example-csv').addEventListener("click", function() {
  var dataUrl = "https://raw.githubusercontent.com/ericmuckley/datasets/master/qm9-small.csv";
  d3.csv(dataUrl).then(function(data) {
    load_data(data);
  });
}, false);

// when file is selected from upload field
dom.get("file-upload").addEventListener("change", function() {
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
  cleanPage();

  if (data.length > 2500){
    alert("That's a lot of data. We will use only the first 10,000 rows to ensure performance.")
    data = data.slice(0, 2500);
  };
  // embed smiles list in frontend element
  const smilesList = data.map(x => x.smiles);
  dom.get("file-upload").dataset.smilesList = JSON.stringify(smilesList);

  // use first molecule in list as default selected molecule
  selectMolecule(smilesList[0]);

  dom.get("details-div").style.display = "block";
  dom.get('bin-select').addEventListener('change', plotHistogram);

};



function selectMolecule(s){
  // save new molecule smiles string
  dom.get("file-upload").dataset.selected = s;

  // draw the selected molecule
  postData(url="/draw_smiles", data={"smiles": s}).then(x => {
    dom.get('selected-mol-header').innerHTML = s;
    var img = svgToImg(x.svg);
    img.width = 200;
    dom.hideLoading('selected-mol-div');
    dom.get('selected-mol-div').appendChild(img);
  });

  // get the molecule similarities
  const similarities = getSimilarities();
};


function getPCA() {
  const smilesList = JSON.parse(dom.get("file-upload").dataset.smilesList);
  postData(url="/get_pca", data={"smiles_list": smilesList}).then(x => {
    const pca = x.pca;
    dom.get("file-upload").dataset.pca = JSON.stringify(pca);
    plotPCA();
  });
};


function getSimilarities() {
  const smiles0 = dom.get("file-upload").dataset.selected;
  const smilesList = JSON.parse(dom.get("file-upload").dataset.smilesList);
  postData(url="/get_similarities", data={"smiles0": smiles0, "smiles_list": smilesList}).then(x => {
    const similarities = x.similarities;
    dom.get("file-upload").dataset.similarities = JSON.stringify(similarities);
    plotHistogram();
    getPCA();
    getSiblingMolecules();
  });
};



function getSiblingMolecules(){

  //const smiles0 = dom.get("file-upload").dataset.selected;
  //const pca = JSON.parse(dom.get("file-upload").dataset.pca);
  const smilesList = JSON.parse(dom.get("file-upload").dataset.smilesList);
  const similarities = JSON.parse(dom.get("file-upload").dataset.similarities);

  // sort all molecules by theiir similarities
  const info = [];
  for (let i=0; i < smilesList.length; i++){
    info.push({smiles: smilesList[i], similarity: similarities[i], idx: i});
  };
  info.sort(function(a, b) {
    return ((a.similarity < b.similarity) ? -1 : ((a.similarity == b.similarity) ? 0 : 1));
  });
  const Nsims = 4;
  const lowSims = info.slice(0, Nsims+1);
  const highSims = info.slice(-Nsims-2, -1);

  console.log(lowSims);
  console.log(highSims);

};





function cleanPage(){
  //dom.get('selected-mol-div').innerHTML = '';
  const divIds = [
    "pca-plot",
    "histogram-plot",
    "selected-mol-div",
  ]
  divIds.forEach(function(id) {
    dom.get(id).innerHTML = "";
    dom.showLoading(id);
  });
};


function svgToImg(svg){
  // convert an SVG string to an image object
  const blob = new Blob([svg], {type: 'image/svg+xml'});
  const url = URL.createObjectURL(blob);
  const img = dom.make('img', {src: url});
  img.addEventListener('load', () => URL.revokeObjectURL(url), {once: true});
  return img;
};


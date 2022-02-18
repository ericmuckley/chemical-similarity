  // show all tooltips
  var tooltipTriggerList = [].slice.call(document.querySelectorAll('[data-bs-toggle="tooltip"]'))
  var tooltipList = tooltipTriggerList.map(function (tooltipTriggerEl) {
    return new bootstrap.Tooltip(tooltipTriggerEl)
  })


// when user clicks on example file link
dom.get('use-example-csv').addEventListener("click", function() {
  var dataUrl = "https://raw.githubusercontent.com/ericmuckley/datasets/master/qm9-small.csv";
  d3.csv(dataUrl).then(function(data) {
    load_data(data);
  });
}, false);

// when file is selected from upload field
dom.get("file-upload").addEventListener("input", function() {
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
  clearPage();

  // check if smiles header was the wrong case
  if (! Object.keys(data[0]).includes('smiles')){
    for (let k of Object.keys(data[0])){
      if (k.toLowerCase() === 'smiles'){
        for (let d0 of data){
          d0.smiles = d0[k];
        };
      };
    };
  };

  const smilesList = data.map(x => x.smiles);

  if (smilesList.length < 2){
    alert("File does not contain a 'smiles' column header with valid smiles strings in rows.")
  } else {
    if (data.length > 2500){
      alert("That's a lot of data. We will use only the first 10,000 rows to ensure performance.")
      data = data.slice(0, 2500);
    };
    // embed smiles list in frontend element
    dom.get("file-upload").dataset.smilesList = JSON.stringify(smilesList);
    // use first molecule in list as default selected molecule
    selectMolecule(smilesList[0]);
    dom.get('bin-select').addEventListener('change', plotHistogram);
  };
};



function selectMolecule(s){

  ["similar-div", "dissimilar-div"].forEach(function(id) {
    dom.get(id).innerHTML = "";
  })

  // save new molecule smiles string
  dom.get("file-upload").dataset.selected = s;

  // draw the selected molecule
  postData(url="/draw_smiles", data={"smiles": s}).then(x => {
    if (x.msg === "success"){
      dom.get('selected-mol-header').innerHTML = s;
      var img = svgToImg(x.svg);
      img.width = 200;
      dom.hideLoading('selected-mol-div');
      dom.get('selected-mol-div').appendChild(img);
      getSimilarities();      
    } else {
      clearPage();
      alert("File does not contain valid smiles strings");
    };
  });


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
  dom.get("details-div").style.display = "block";
  const smiles0 = dom.get("file-upload").dataset.selected;
  const smilesList = JSON.parse(dom.get("file-upload").dataset.smilesList);
  postData(url="/get_similarities", data={"smiles0": smiles0, "smiles_list": smilesList}).then(x => {
    if (x.msg === "fail"){
      clearPage();
      alert('Some SMILES strings are not valid');
    } else {
      const similarities = x.similarities;
      dom.get("file-upload").dataset.similarities = JSON.stringify(similarities);
      plotHistogram();
      getSiblingMolecules();
      if (dom.get('pca-plot').dataset.initialized === "true"){
        plotPCA();
      } else {
        getPCA();
      };
    };
   
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


  var Nsims = Math.floor(similarities.length/10) + 1;
  Nsims = (Nsims > 12) ? 12 : Nsims;
  // get high and low similarity samples
  const lowSims = info.slice(0, Nsims+1);
  const highSims = info.slice(-Nsims-2, -1);

  // regions of molecules: similar or dissimilar
  const regions = [
    {data: highSims, div: 'similar-div', color: '#198754', hoverclass: "link-div-similar", code: "success"},
    {data: lowSims, div: 'dissimilar-div', color: '#dc3545', hoverclass: "link-div-dissimilar", code: "danger"},
  ];

  // iterate over each region
  for (let reg of regions){
    // iterate over each molecule in the region
    for (let m of reg.data){
      // get molecule image
      postData(url="/draw_smiles", data={"smiles": m.smiles}).then(x => {
        const div = dom.make('div', {
          classes:`text-center m-2 cdiv p-3 ${reg.hoverclass}`,
          parent: reg.div,
          children: [["p", {classes: `text-center text-${reg.code}`, innerHTML: m.smiles}]],
        });
        div.style.border = `1px solid ${reg.color}`;
        div.addEventListener('click', () => {
          selectMolecule(m.smiles);
        });
        // add image to page
        //svg = x.svg.replaceAll("fill:#FFFFFF", "fill:#198754");
        var img = svgToImg(x.svg);
        img.width = 60;
        div.appendChild(img);
      });
    };
  };

};





function clearPage(){
  //dom.get('selected-mol-div').innerHTML = '';
  dom.get('details-div').style.display = 'none';
  dom.get('pca-plot').dataset.initialized = "false"
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


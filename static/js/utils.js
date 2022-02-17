
// get Pearson r^2 correlation coefficient
function getR2(x, y){
  let sumX = 0,
    sumY = 0,
    sumXY = 0,
    sumX2 = 0,
    sumY2 = 0;
  const minLength = x.length = y.length = Math.min(x.length, y.length),
    reduce = (xi, idx) => {
      const yi = y[idx];
      sumX += xi;
      sumY += yi;
      sumXY += xi * yi;
      sumX2 += xi * xi;
      sumY2 += yi * yi;
    }
  x.forEach(reduce);
  var R = (minLength * sumXY - sumX * sumY) / Math.sqrt((minLength * sumX2 - sumX * sumX) * (minLength * sumY2 - sumY * sumY));
  return Math.pow(R, 2);
};


// create a select input using an array of options
function createSelect(options, label='', id='', parentId='', classList=[]){
  var select = document.createElement("select");
  select.classList.add('form-select');
  for (let c of classList){
    select.classList.add(c);
  };
  if (id !== ''){
      select.id = id
  };
  if (label !== '' && parentId !== ''){
      const lab = document.createElement('label');
      lab.classList.add('form-label', 'mb-0');
      lab.innerHTML = label;
      document.getElementById(parentId).appendChild(lab);
  };
  if (parentId !== ''){
      document.getElementById(parentId).appendChild(select)
  };
  for (let op of options){
      var option = document.createElement("option");
      //option.value = strToId(op);
      option.value = op;
      option.text = op;
      select.appendChild(option);
  };
  return select;
};


function clearSelectOptions(selectId){
  var select = document.getElementById(selectId);
  while (select.options.length) {                
    select.remove(0);
  };
};


function populateSelect(selectId, options, selectedIndex=0){
  clearSelectOptions(selectId)
  var select = document.getElementById(selectId);
  for (let op of options){
    var option = document.createElement("option");
    option.value = op;
    option.text = op;
    select.appendChild(option);
  };
  select.value = options[selectedIndex];
};


// Create a bootstrap row with multiple columns.
// Returns an array containing each column.
function createRowCols(parentdiv, ncols=2) {
  var row = document.createElement('div');
  row.classList.add('row');
  for (let i=0; i<ncols; i++) {
      var col = document.createElement('div');
      col.classList.add('col');
      row.appendChild(col);
  };
  parentdiv.appendChild(row);
  return row
};


// Create and populate a table from an array of data.
// The *data* argument should be an array of objects.
function createFilledTable(data=null, classes=['table', 'table-sm', 'table-borderless'], divId=null){
  // create tabel
  const table = document.createElement("table");
  var thead = document.createElement('thead');
  var tbody = document.createElement('tbody');
  table.appendChild(thead);
  table.appendChild(tbody);
  for (var c of classes){
      table.classList.add(c);
  };
  if (data !== null){
      // add table headers
      var headers = Object.keys(data[0]);
      var headrow = thead.insertRow();
      for (var h of headers) {
          var th = document.createElement('th');
          //th.style.textAlign = 'center';
          //th.classList.add('text-nowrap');
          th.innerHTML = h;
          headrow.appendChild(th);
      };
      // add table body content
      for (let record of data){
          var row = tbody.insertRow();
          for (let h of headers){
              var cell = row.insertCell();
              cell.innerHTML = record[h];
          };
      };
  };
  if (divId !== null) {
      document.getElementById(divId).appendChild(table);
  };
  return [table, thead, tbody];
};




// add a spinner to a parent div
function addSpinner(parentId, color='primary', size=null){
  const spinner = document.createElement('div');
  spinner.classList.add('mx-2', 'spinner', 'spinner-border', `text-${color}`);
  if (size !== null){
      spinner.classList.add(`spinner-border-${size}`);
  };
  document.getElementById(parentId).appendChild(spinner);
};
// remove spinner from a parent div
function removeSpinner(parentId){
  const spinners = document.getElementById(parentId).getElementsByClassName("spinner");
  for(const s of spinners){
      s.remove()
  };
};


// from a hex color code and an alpha level, return a new color code
function addAlpha(color, alpha) {
  var _alpha = Math.round(Math.min(Math.max(alpha || 1, 0), 1) * 255);
  return color + _alpha.toString(16).toUpperCase();
};

// generate a random ID of strings
function getRandomId(n=12) {
    var chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghiklmnopqrstuvwxyz'.split('');
    if (! n) {
        n = Math.floor(Math.random() * chars.length);
    };
    var id = '';
    for (var i = 0; i < n; i++) {
        id += chars[Math.floor(Math.random() * chars.length)];
    };
    return id;
};

// check if a string is a number
function stringIsNumeric(str) {
  return !isNaN(str) && !isNaN(parseFloat(str));
};






function createWindow(title="My window title", content="Window content", loc=[30, 30]) {

  console.log(title);
  // create window
  var win = document.createElement('div');
  Object.assign(win.style, {
    position: 'absolute',
    top: `${loc[0]}px`,
    left: `${loc[1]}px`,
    backgroundColor: 'white',
    boxShadow: '3px 3px 10px 3px rgba(0,0,0,0.5)',
    //zIndex: 100,
    resize: 'both',
    overflow: 'auto',
    minWidth: '250px',
    minHeight: '50px',
  });

  // add window to main
  document.getElementsByTagName('main')[0].appendChild(win);

  // create window top title bar
  var titleBar = document.createElement('div');
  titleBar.classList.add('d-flex', 'justify-content-between');
  win.appendChild(titleBar);
  Object.assign(titleBar.style, {
    padding: '5px',
    cursor: 'move',
    marginBottom: '5px',
    borderBottom: '1px solid #ccc',
    backgroundColor: '#dee2e6',
  });
  // create icon on title bar
  var titleBarIcon = document.createElement('span');
  titleBarIcon.innerHTML = `<i class="bi bi-arrows-move"></i>`;
  titleBar.appendChild(titleBarIcon);
  // create window title on title bar
  var titleDiv = document.createElement('span');
  titleDiv.innerHTML = `<strong>${title}</strong>`;
  titleBar.appendChild(titleDiv);
  // create close window button on title bar
  var closeBtn = document.createElement('button');
  closeBtn.type = 'button';
  closeBtn.classList.add("btn-close")
  titleBar.appendChild(closeBtn);
  Object.assign(closeBtn.style, {
    cursor: 'pointer',
  });
  closeBtn.addEventListener('click', () => {
    //win.style.display = "none";
    closeBtn.parentNode.parentNode.remove();
  }, false);

  // create main window content
  var windowContent = document.createElement('div');
  win.appendChild(windowContent);
  Object.assign(windowContent.style, {
    padding: '5px',
  });
  windowContent.innerHTML = content;

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
};



// generate random string to use as a DOM ID
function getRandomId(n=12) {
  const chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghiklmnopqrstuvwxyz'.split('');
  var id = '';
  for (let i = 0; i < n; i++) {
      id += chars[Math.floor(Math.random() * chars.length)];
  };
  return id;
};
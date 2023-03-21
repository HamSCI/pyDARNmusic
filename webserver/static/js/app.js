// source: https://www.codinglabweb.com/2021/06/dropdown-sidebar-menu-html-css.html
let arrow = document.querySelectorAll(".arrow");
  for (var i = 0; i < arrow.length; i++) {
    arrow[i].addEventListener("click", (e)=>{
   let arrowParent = e.target.parentElement.parentElement;//selecting main parent of arrow
   arrowParent.classList.toggle("showMenu");
    });
  }
  let sidebar = document.querySelector(".sidebar");
  let sidebarBtn = document.querySelector(".bx-menu");
  console.log(sidebarBtn);
  sidebarBtn.addEventListener("click", ()=>{
    sidebar.classList.toggle("close");
  });

  //show loading screen
var pageContent = document.getElementById("classify-content");
var loadingScreen = document.getElementById("loading_screen")
var btnClassification = document.getElementById("classification-btn")

btnClassification.onclick = function () {
  loadingScreen.classList.remove("dis-non");
}

  // hemisphere selection buttons
var table1 = document.getElementById("north");
var table2 = document.getElementById("south");

var btnTab1 = document.getElementById("northernhemisphere");
var btnTab2 = document.getElementById("southernhemisphere");

btnTab1.onclick = function() {
  table1.style.display = "table";
  table2.style.display = "none";
  btnTab2.style.backgroundColor = "#6c757d";
  btnTab1.style.backgroundColor = "#041c32";
}

btnTab2.onclick = function() {
  table1.style.display = "none";
  table2.style.display = "table";
  btnTab1.style.backgroundColor = "#6c757d";
  btnTab2.style.backgroundColor = "#041c32";
}
/* eslint-disable @blueprintjs/classes-constants -- do not import files */
const dataset =
  typeof document !== "undefined" &&
  document.querySelector("#obsolete-browsers").dataset;

const datasetRegex = dataset.regex;

const regex = new RegExp(datasetRegex.slice(1, -1));

if (dataset && !regex.test(navigator.userAgent)) {
  const root = document.getElementById("root");
  root.remove();
  const portals = document.getElementsByClassName("bp3-portal");
  for (let i = 0; i < portals.length; i += 1) {
    portals[i].remove();
  }

  const element = document.createElement("div");
  element.innerHTML = dataset.template;
  document.body.appendChild(element);
}
/* eslint-enable @blueprintjs/classes-constants -- do not import files */

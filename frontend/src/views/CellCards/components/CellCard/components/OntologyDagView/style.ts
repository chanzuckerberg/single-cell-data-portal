import styled from "@emotion/styled";

export const Wrapper = styled.div`
  .legendText {
    font: 10px sans-serif;
  }
  .axis {
    font: 10px sans-serif;
  }
  .axis path,
  .axis line {
    fill: none;
    stroke: #000;
    shape-rendering: crispEdges;
  }

  .outlined {
    fill: none;
    stroke-width: 4;
  }

  .node {
    cursor: pointer;
  }

  .node circle {
    fill: #fff;
    stroke-width: 1.5px;
  }

  .node text {
    font: 11px sans-serif;
  }

  .overlay {
    background-color: #eee;
  }

  .link {
    fill: none;
    stroke: #ccc;
    stroke-width: 1.5px;
  }

  .search {
    width: 30%;
  }

  .button {
    background-color: #0073ff;
    border: none;
    color: white;
    padding: 10px 15px;
    text-align: center;
    text-decoration: none;
    display: inline-block;
    font-size: 10px;
    margin: 5px 10px;
    cursor: pointer;
    border-radius: 8px;
  }
  .hover {
    background: lightsteelblue;
    font: 12px sans-serif;
    position: absolute;
    text-align: center;
    width: 80px;
    /*height: 54px;*/
    padding: 2px;
    font: 12px sans-serif;
    border: 0px;
    border-radius: 8px;
    pointer-events: none;
    opacity: 80%;
  }
  .toggle-label {
    font: 12px sans-serif;
  }
  h1 {
    font-family: Copperplate;
    text-align: center;
    color: black;
  }
`;

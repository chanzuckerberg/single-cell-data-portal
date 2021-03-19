import styled from "styled-components";

export const Layout = styled.div`
  width: 100%;
  background-color: white;
`;

export const CommonStyle = styled.div`
  max-width: 700px;
  margin: 0px auto;
  padding: 50px;
  line-height: 1.5;
  font-family: "Roboto", sans-serif;

  h4,
  h5 {
    display: inline;
    font-size: 16px;
  }

  .caps {
    text-transform: uppercase;
  }
`;

export const CZILogo = styled.img`
  height: 40px;
  padding-right: 12px;
  border-right-color: black;
  border-right-style: solid;
  border-right-width: 2px;
`;

export const CellxgeneLogo = styled.img`
  height: 35px;
`;

export const PrivacyStyle = styled.div`
  h1,
  h2,
  h3 {
    margin: 30px 0 10px 0;
  }

  h1,
  h2,
  h3,
  h4,
  h5,
  h6,
  p,
  ul,
  ol {
    margin: 0;
    margin-block-start: 0;
    margin-block-end: 0;
  }

  p,
  li {
    margin-top: 10px;
  }

  ol > li,
  ul > li,
  ol > p,
  ul > p {
    padding-left: 8px;
  }

  .section1 > li,
  .section2 > li,
  .section3 > li,
  .section4 > li,
  .section5 > li,
  .section6 > li,
  .section7 > li,
  .section8 > li,
  .section9 > li,
  .section10 > li {
    counter-increment: increase-by;
  }

  ol.section1,
  ol.section2,
  ol.section3,
  ol.section4,
  ol.section5,
  ol.section6,
  ol.section7,
  ol.section8,
  ol.section9,
  ol.section10 {
    list-style-type: none;
  }

  .section1 > li:before,
  .section2 > li:before,
  .section3 > li:before,
  .section4 > li:before,
  .section5 > li:before,
  .section8 > li:before,
  .section9 > li:before,
  .section10 > li:before {
    padding-right: 8px;
  }

  .section2 > li:before {
    content: "2." counter(increase-by, decimal) " ";
  }
  .section3 > li:before {
    content: "3." counter(increase-by, decimal) " ";
  }
  .section4 > li:before {
    content: "4." counter(increase-by, decimal) " ";
  }
  .section5 > li:before {
    content: "5." counter(increase-by, decimal) " ";
  }
  .section6 > li:before {
    content: "6." counter(increase-by, decimal) " ";
  }
  .section7 > li:before {
    content: "7." counter(increase-by, decimal) " ";
  }
  .section8 > li:before {
    content: "8." counter(increase-by, decimal) " ";
  }
  .section9 > li:before {
    content: "9." counter(increase-by, decimal) " ";
  }
  .section10 > li:before {
    content: "10." counter(increase-by, decimal) " ";
  }
`;

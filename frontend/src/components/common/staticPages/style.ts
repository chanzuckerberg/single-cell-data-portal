import styled from "@emotion/styled";

export const Layout = styled.div`
  width: 100%;
  background-color: white;
`;

export const CommonStyle = styled.div`
  max-width: 700px;
  margin: 0px auto;
  padding: 50px;
  line-height: 1.5;
  font-family: "Inter", sans-serif;

  h4,
  h5 {
    display: inline;
    font-size: 16px;
  }

  .caps {
    text-transform: uppercase;
  }
`;

const AgreementDocumentStyle = styled.div`
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
`;

export const TOSStyle = styled(AgreementDocumentStyle)`
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
  .section9 > li {
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
  ol.section9 {
    list-style-type: none;
  }

  .section1 > li:before,
  .section2 > li:before,
  .section3 > li:before,
  .section4 > li:before,
  .section5 > li:before,
  .section6 > li:before,
  .section7 > li:before,
  .section8 > li:before,
  .section9 > li:before {
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
`;

export const PrivacyStyle = styled(AgreementDocumentStyle)`
  // "Last updated" text.
  h1 + p {
    margin-bottom: 20px;
  }

  h2,
  h3,
  h4 {
    margin: 10px 0;
  }

  main > h4 {
    display: block;
  }

  // Inline heading and paragraph.
  .inline {
    margin-top: 10px;
  }

  .inline > * {
    display: inline;
  }

  ul {
    padding-left: 16px;
  }

  ul.text-list {
    padding-left: 8px;
  }

  ul > li {
    list-style-type: square;
    padding: 0;
  }

  ul.text-list > li {
    list-style-type: none;
  }
`;

export const LandingStyle = styled.div`
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

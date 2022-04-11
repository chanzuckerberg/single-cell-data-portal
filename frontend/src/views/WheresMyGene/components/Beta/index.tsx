import { StyledCallout, SubmitIssue } from "./style";

export default function Beta(): JSX.Element {
  return (
    <StyledCallout intent="info" dismissed={false}>
      This feature is in beta. We would appreciate your feedback, please fill
      out a{" "}
      <SubmitIssue
        href="https://docs.google.com/forms/d/e/1FAIpQLSde_zIFZPQD2p0ovaX3Pb7lDOajWJCmOeuX4wQ8Z8Ab5NXUjw/viewform"
        target="_blank"
        rel="noopener"
      >
        quick survey
      </SubmitIssue>
      .
    </StyledCallout>
  );
}

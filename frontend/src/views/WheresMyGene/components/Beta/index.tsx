import { StyledCallout, SubmitIssue } from "./style";

export default function Beta(): JSX.Element {
  return (
    <StyledCallout intent="info" dismissed={false}>
      This feature is in beta. If you have any suggestions or feedback, please{" "}
      <SubmitIssue
        href="https://github.com/chanzuckerberg/single-cell-data-portal/issues/new/choose"
        target="_blank"
        rel="noopener"
      >
        submit an issue on our GitHub
      </SubmitIssue>
      .
    </StyledCallout>
  );
}

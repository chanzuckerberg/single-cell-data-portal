import { noop } from "src/common/constants/utils";
import { StyledNotification, SubmitIssue } from "./style";

export default function Beta({
  className,
}: {
  className?: string;
}): JSX.Element {
  const position = {
    bottom: "10px",
    position: "absolute",
    right: "20px",
    width: "360px",
    zIndex: "99",
  } as React.CSSProperties;

  return (
    <StyledNotification
      id="beta-notification"
      intent="info"
      autoDismiss={false}
      dismissDirection="left"
      onClose={noop}
      style={position}
      className={"elevated " + className}
    >
      We would appreciate your feedback, please fill out a{" "}
      <SubmitIssue
        href="https://airtable.com/shrLwepDSEX1HI6bo"
        target="_blank"
        rel="noopener"
      >
        quick survey
      </SubmitIssue>
      .
    </StyledNotification>
  );
}

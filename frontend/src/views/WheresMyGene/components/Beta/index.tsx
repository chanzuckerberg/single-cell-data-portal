import { noop } from "src/common/constants/utils";
import { StyledNotification, SubmitIssue } from "./style";

export default function Beta({
  className,
}: {
  className?: string;
}): JSX.Element {
  const position = {
    position: "absolute",
    bottom: "10px",
    right: "20px",
    zIndex: "99",
  } as React.CSSProperties;

  return (
    <StyledNotification
      intent="info"
      autoDismiss={false}
      dismissDirection="left"
      onClose={noop}
      style={position}
      className={"elevated " + className}
    >
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
    </StyledNotification>
  );
}

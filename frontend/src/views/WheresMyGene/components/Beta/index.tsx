import { StyledNotification, SubmitIssue } from "./style";

export default function Beta(): JSX.Element {
  const position = {
    "position": "absolute",
    "bottom": "60px",
    "right": "15px",
    "z-index": "1"
  } as React.CSSProperties;

  return (
    <StyledNotification 
      intent="info" 
      autoDismiss={false} 
      dismissDirection="left" 
      onClose={function(){}} 
      style={position}
      className="elevated">
      This feature is in beta. If you have any suggestions or feeback, please{" "}
      <SubmitIssue
        href="https://docs.google.com/forms/d/e/1FAIpQLSde_zIFZPQD2p0ovaX3Pb7lDOajWJCmOeuX4wQ8Z8Ab5NXUjw/viewform"
        target="_blank"
        rel="noopener"
      >
        submit an issue on our GitHub
      </SubmitIssue>
      .
    </StyledNotification>
  );
}

import React from "react";
import {
  StyledCallout,
  CalloutTitle,
  CalloutTextWrapper,
  StyledButton,
} from "./style";

const ByodCallout: React.FC = () => {
  const handleLearnMore = () => {
    // TODO: Implement BYOD info dialog
    console.log("Open BYOD info dialog");
  };

  return (
    <StyledCallout intent="info">
      <CalloutTextWrapper>
        <CalloutTitle>
          Want to use this model to explore your own data?
        </CalloutTitle>{" "}
        Upload and analyze your data, add annotations, and explore embeddings on
        CZI&apos;s AI Workspace.
      </CalloutTextWrapper>
      <StyledButton
        sdsType="primary"
        sdsStyle="minimal"
        onClick={handleLearnMore}
      >
        Learn more about AI Workspace
      </StyledButton>
    </StyledCallout>
  );
};

export default ByodCallout;

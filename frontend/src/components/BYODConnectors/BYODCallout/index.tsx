import React from "react";
import { StyledButton } from "../style";
import { useBYODModal } from "src/contexts/BYODModalContext";
import SparkleIcon from "src/common/images/sparkle-icon.svg";
import { CalloutTextWrapper, CalloutTitle, StyledCallout } from "./style";

const BYODCallout: React.FC = () => {
  const { openModal } = useBYODModal();

  return (
    <StyledCallout intent="info" icon={<SparkleIcon />}>
      <CalloutTextWrapper>
        <CalloutTitle>
          Want to use this model to explore your own data?
        </CalloutTitle>{" "}
        Upload and analyze your data, add annotations, and explore embeddings on
        CZI&apos;s AI Workspace.
      </CalloutTextWrapper>
      <StyledButton sdsType="primary" sdsStyle="minimal" onClick={openModal}>
        Learn more about AI Workspace
      </StyledButton>
    </StyledCallout>
  );
};

export default BYODCallout;

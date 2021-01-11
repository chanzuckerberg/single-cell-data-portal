import {
  Icon,
  Intent,
  PopoverInteractionKind,
  Tooltip as TooltipRaw,
} from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import React, { FC } from "react";
import { UPLOAD_STATUS, VALIDATION_STATUS } from "src/common/entities";
import { RED } from "src/components/common/theme";
import { StyledAnchor } from "./style";

const tooltipContent = (error: VALIDATION_STATUS | UPLOAD_STATUS) => {
  return error === VALIDATION_STATUS.INVALID ? (
    <span>
      You must validate your dataset locally before uploading. We provide a
      local CLI script to do this.{" "}
      <b>
        <StyledAnchor
          href="https://github.com/chanzuckerberg/cellxgene/blob/main/dev_docs/schema_guide.md"
          target="_blank"
        >
          Learn More
        </StyledAnchor>
      </b>
    </span>
  ) : (
    <span>There was a problem uploading your file. Please try again.</span>
  );
};

interface Props {
  error?: VALIDATION_STATUS | UPLOAD_STATUS;
}

const Tooltip: FC<Props> = ({ error }) => {
  if (!error) return null;

  return (
    <TooltipRaw
      intent={Intent.DANGER}
      interactionKind={PopoverInteractionKind.HOVER}
      hoverCloseDelay={500}
      content={tooltipContent(error)}
    >
      <Icon icon={IconNames.ISSUE} iconSize={16} color={RED.C} />
    </TooltipRaw>
  );
};

export default Tooltip;

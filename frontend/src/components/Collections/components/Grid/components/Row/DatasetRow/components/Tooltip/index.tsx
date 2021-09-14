import {
  Icon,
  Intent,
  PopoverInteractionKind,
  Tooltip as TooltipRaw,
} from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { FC } from "react";
import {
  CONVERSION_STATUS,
  UPLOAD_STATUS,
  VALIDATION_STATUS,
} from "src/common/entities";
import { ORANGE, RED } from "src/components/common/theme";
import { FAILED_RETURN_TYPE, FailReturn } from "../../utils";
import { StyledAnchor } from "./style";

interface Content {
  color: string;
  intent: Intent;
  content: JSX.Element;
}

const ERROR_TO_CONTENT: { [key: string]: Content } = {
  [FAILED_RETURN_TYPE.VALIDATION + VALIDATION_STATUS.INVALID]: {
    color: RED.C,
    content: (
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
    ),
    intent: Intent.DANGER,
  },
  [FAILED_RETURN_TYPE.CONVERSION + CONVERSION_STATUS.FAILED]: {
    color: ORANGE.C,
    content: (
      <span>
        The dataset was uploaded successfully, but one or more conversions from
        .h5ad failed. We&apos;ll attempt to fix this manually and follow-up in
        an email after we&apos;ve investigated.
      </span>
    ),
    intent: Intent.WARNING,
  },
  [FAILED_RETURN_TYPE.UPLOAD + UPLOAD_STATUS.FAILED]: {
    color: RED.C,
    content: <span>There was a problem uploading your file.</span>,
    intent: Intent.DANGER,
  },
};

const DEFAULT_CONTENT: Content = {
  color: RED.C,
  content: <span>There was an unexpected problem.</span>,
  intent: Intent.DANGER,
};

interface Props {
  error: FailReturn["error"];
  type: FailReturn["type"];
}

const Tooltip: FC<Props> = ({ error, type }) => {
  if (!error) return null;

  const content = ERROR_TO_CONTENT[type + error] || DEFAULT_CONTENT;

  return (
    <TooltipRaw
      intent={content.intent}
      interactionKind={PopoverInteractionKind.HOVER}
      hoverCloseDelay={500}
      content={content.content}
    >
      <Icon icon={IconNames.ISSUE} iconSize={16} color={content.color} />
    </TooltipRaw>
  );
};

export default Tooltip;
